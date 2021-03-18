library(REDCapR)
library(readxl)
library(writexl)

# Working directory
setwd("~/Projects/SerocoViD")

# ------------------------------- REDCap data ------------------------------- #

# API
uri   <- "https://redcap.unisante.ch/api/"

# Tokens
tokens <- list(
  corona_immunitas = "F1BFC49C8D21C355AE9946277DD52EDC",
  personal_data    = "2553F65D05A7FB4FB8AAACE5BDE0784C",
  research_data    = "C296C589D5FB1A89DBFC11D14EFA4365"
)

# Import
tmp_file <- "/tmp/serocovid_data.rda"
if (file.exists(tmp_file)) {
  load(tmp_file)
} else {
  serocovid_data <- lapply(tokens, function(token) {
    redcap_read(redcap_uri = uri, token = token)$data
  })
  save(serocovid_data, file = tmp_file)
}
rm(tmp_file)

# ------------------- Date of birth / Stratum - Survey 3 -------------------- #

# Date of birth
dob <- read_xlsx("data-fso/COVID19_VD_V2_REDCap_et_Publipostage.xlsx",
                 sheet = "Feuil1", guess_max = 10^4)
dob <- as.data.frame(dob[c("HID", "dateOfBirth")])
names(dob)[names(dob) == "HID"] <- "hid"
dob$dateOfBirth <- as.Date(dob$dateOfBirth)
if (any(duplicated(dob$hid))) stop("duplicated hid")
if (any(is.na(dob$dateOfBirth))) stop("missing date of birth")

# Stratum
smpl <- read.csv("data-fso/COVID19_VD_V2_Total.csv", sep = ";")
smpl$dateOfBirth <- as.Date(smpl$dateOfBirth)
tmp <- unique(smpl[!is.na(smpl$strate), c("dateOfBirth", "strate")])
tmp <- do.call(rbind, lapply(sort(unique(tmp$strate)), function(i) {
  z <- tmp[tmp$strate == i, "dateOfBirth"]
  data.frame(dateOfBirth = seq(min(z), max(z), by = 1), stratum = i)
}))
if (!all(aggregate(stratum ~ dateOfBirth, tmp, function(x)
           length(unique(x)))[, 2] == 1)) {
  stop("ill-defined strata")
}

# Merge: get a stratum for each hid
dob <- merge(dob, tmp, by = "dateOfBirth", all.x = TRUE)
if (any(is.na(dob$stratum))) stop("missing stratum")
if (any(duplicated(dob$hid))) stop("duplicated hid")

# ------------------------ Research data - Survey 3 ------------------------- #

# Personal data - type of participant and date of birth
vars <- c("uc_s_participant_hid", "uc_s_type_participant", "full_dob")
pdta <- unique(serocovid_data$personal_data[vars])
if (any(duplicated(pdta$uc_s_participant_hid))) {
  stop("Duplicated hid")
}

# Research data
vars <- c("hid", "uc_labo_coviggl_v2", "uc_labo_covigal_v2",
          "uc_ser_date_hour_v2")
dta <- serocovid_data$corona_immunitas[vars]
dta <- merge(dta, pdta, by.x = "hid", by.y = "uc_s_participant_hid",
             all.x = TRUE)
if (any(is.na(dta$uc_s_type_participant))) {
  stop("missing type of participant")
}
rm(vars, pdta)

# Select observation of survey 3
dta <- dta[dta$uc_s_type_participant == 13, ]

# Remove fake participants
dta <- dta[!grepl("^123123", dta$hid), ]

# -------------------------------- Serology --------------------------------- #

# Serology
if (FALSE) {
  with(dta, table(uc_labo_coviggl_v2, uc_labo_covigal_v2, useNA = "ifany"))
}
for (k in 1:2) {
  u <- c("uc_labo_coviggl_v2", "uc_labo_covigal_v2")[k]
  v <- c("serol_igg", "serol_iga")[k]
  dta[[v]] <- NA
  i <- !is.na(dta[[u]]) & grepl("(Négatif|Ind(e|é)terminé)", dta[[u]])
  dta[i, v] <- 0
  i <- !is.na(dta[[u]]) & grepl("Positif", dta[[u]])
  dta[i, v] <- 1
  if (all(is.na(dta[[v]]) | dta[[v]] %in% c("0", "1"))) {
    dta[[v]] <- as.numeric(dta[[v]])
  } else {
    stop("error in serology recoding")
  }
}
dta$serol_any <- pmin(dta$serol_igg + dta$serol_iga, 1)
rm(k, i, u, v)

# Select the observation for which serology is available
dta <- dta[!is.na(dta$serol_igg) | !is.na(dta$serol_iga), ]
with(dta, addmargins(table(serol_igg, serol_iga, useNA = "ifany")))

# Add stratum to research data
dta <- merge(dta, dob, by = "hid", all.x = TRUE)
if (any(is.na(dta$stratum))) stop("missing stratum")

# ------------------------------- Estimation -------------------------------- #

# Population totals
pop3 <- read.table(header = TRUE, text = "
  stratum              label      N
        1  '15 ans - 19 ans'  42140
        2  '20 ans - 39 ans' 214498
        3  '40 ans - 64 ans' 270755
        4  '65 ans - 74 ans'  67462
        5      '75 ans et +'  63605
")

# Estimation (function)
# Variance: Cochran, W. G., Sampling techniques John Wiley & Sons, 1977,
#           page 144, formula 5A.75
# domain: any variable
Mean <- function(data = data, variable = "serol", stratum = "stratum",
                 domain = NULL, pop = pop) {
  # select observed values
  data <- data[!is.na(data[[variable]]) & !is.na(data[[stratum]]), ]
  # formula variable ~ stratum
  f <- as.formula(paste(variable, "~", stratum))
  # stratum sizes
  strata <- aggregate(f, data, length)
  names(strata) <- c("stratum", "n")
  strata <- merge(pop, strata, by = "stratum")
  # convert N and n to numeric to avoid the 'integer overflow' problem
  strata$N <- as.numeric(strata$N)
  strata$n <- as.numeric(strata$n)
  # no domain
  if (is.null(domain)) {
    data$.all <- 1
    domain <- ".all"
  }
  # convert domains to factors
  for (d in domain) {
    data[[d]] <- factor(data[[d]])
  }
  .mean <- do.call(rbind, lapply(domain, function(d) {
    do.call(rbind, lapply(levels(data[[d]]), function(z) {
      data.sub <- data[data[[d]] == z, c(stratum, variable)]
      # sample size of domain by strata
      m <- aggregate(f, data.sub, length)
      names(m) <- c("stratum", "m")
      # mean by domain and strata
      y <- aggregate(f, data.sub, mean)
      names(y) <- c("stratum", "y")
      # sum of squares by domain and strata
      ss <- aggregate(f, data.sub, function(u) sum((u - mean(u))^2))
      names(ss) <- c("stratum", "ss")
      # merge prevalences and sizes with strata
      strata.sub <- merge(strata, m, by = "stratum")
      strata.sub <- merge(strata.sub, y, by = "stratum")
      strata.sub <- merge(strata.sub, ss, by = "stratum")
      #
      M <- sum(with(strata.sub, N / n * m))
      Y <- sum(with(strata.sub, N / n * m * y)) / M
      V <- sum(with(strata.sub, N * (N - n) / (n * (n - 1)) * 
                                  (ss + m * (1 - m / n) * (y - Y)^2))) / M^2
      data.frame(domain = d, value = z, N = M, n = sum(m[2]), y = Y, v = V)
    }))
  }))
  .mean$lwr <- .mean$y - qnorm(0.975) * sqrt(.mean$v)
  .mean$upr <- .mean$y + qnorm(0.975) * sqrt(.mean$v)
  return(.mean)
}

# Sensitivity and specificity
igg <- c(npos = 343, tp = 338, nneg = 256, fp = 2)
iga <- c(npos = 343, tp = 299, nneg = 256, fp = 4)
ig_ <- c(npos = 343, tp = 339, nneg = 256, fp = 4)

# Stratum sizes
aggregate(hid ~ stratum, dta, length)

# Prevalence
set.seed(666)
nsim <- 10^6
dta$all <- 0
prev <- do.call(rbind, lapply(1:3, function(k) {
  if (k == 1) {
    ab <- "IgG"
    vr <- "serol_igg"
    ig <- igg
  } else if (k == 2) {
    ab <- "IgA"
    vr <- "serol_iga"
    ig <- iga
  } else {
    ab <- "IgG or IgA"
    vr <- "serol_any"
    ig <- ig_
  }
  se <- ig[["tp"]] / ig[["npos"]]
  sp <- 1 - ig[["fp"]] / ig[["nneg"]]
  v_se <- se * (1 - se) / ig[["npos"]]
  v_sp <- sp * (1 - sp) / ig[["nneg"]]
  dom <- c("all", "stratum")
  prev <- Mean(data = dta, variable = vr, domain = dom, pop = pop3)
  prev$prev <- (prev$y + sp - 1) / (se + sp - 1)
  prev_sim <- do.call(rbind, lapply(1:nrow(prev), function(i) {
    if (!is.na(prev[i, "v"]) & !is.na(prev[i, "v"])) {
      y_sim <- rnorm(nsim, prev[i, "y"], sqrt(prev[i, "v"]))
    } else {
      y_sim <- rep(NA, nsim)
    }
    if (!is.na(se) & !is.na(v_se)) {
      se_sim <- rnorm(nsim, se, sqrt(v_se))
    } else {
      se_sim <- rep(NA, nsim)
    }
    if (!is.na(sp) & !is.na(v_sp)) {
      sp_sim <- rnorm(nsim, sp, sqrt(v_sp))
    } else {
      sp_sim <- rep(NA, nsim)
    }
    prev <- (y_sim + sp_sim - 1) / (se_sim + sp_sim - 1)
    quantile(prev, c(.5, .025, .975), na.rm = TRUE)
  }))
  colnames(prev_sim) <- c("prev_sim", "prev_lwr", "prev_upr")
  prev <- cbind(prev, prev_sim)
  prev <- cbind(Antibody = ab, prev)
  prev <- prev[colnames(prev) != "v"]
  prev <- prev[colnames(prev) != "prev_sim"]
  names(prev)[names(prev) == "y"] <- "ppos"
  names(prev)[names(prev) == "lwr"] <- "2.5%"
  names(prev)[names(prev) == "upr"] <- "97.5%"
  names(prev)[names(prev) == "prev_lwr"] <- "2.5%"
  names(prev)[names(prev) == "prev_upr"] <- "97.5%"
  prev
}))
prev$domain <- NULL
prev$value <- factor(prev$value, 0:5,
                     c("all", "15-19", "20-39", "40-64", "65-74", ">=75"))
names(prev)[names(prev) == "value"] <- "stratum"

# Export results
sink("results/prev_3rd_visit_sessionInfo_20210215.txt")
print(sessionInfo(), locale = FALSE)
sink()
write_xlsx(prev, "results/prev_3rd_visit_20210215.xlsx")
