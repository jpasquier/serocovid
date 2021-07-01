library(REDCapR)
library(readxl)
library(writexl)

# Working directory
setwd("~/Projects/SerocoViD")

# ------------------------------- REDCap data ------------------------------- #

# API
uri   <- "https://redcap.unisante.ch/api/"

# Tokens
tokens <- c("corona_immunitas", "personal_data", "research_data")
tokens <- lapply(setNames(tokens, tokens), function(z) {
  z <- paste0("misc/redcap_", z, ".token")
  readChar(z, file.info(z)$size - 1)
})

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
rm(uri, tokens, tmp_file)

# ----------------- Date of birth / Stratum - Survey 1 & 2 ------------------ #

# Date of birth
dob1 <- read.csv("data-raw/fso_sid_hid_dob_v2.csv", sep = ";")
names(dob1)[names(dob1) == "HID"] <- "hid"
i <- is.na(dob1$hid) | is.na(dob1$Date.of.birth)
if (any(i)) dob1[i, ]
i <- !grepl("^[0-3][0-9]\\.[0-1][0-9]\\.(19|20)[0-9]{2}$", dob1$Date.of.birth)
if (any(i)) dob1[i, ]
dob1$Date.of.birth <- as.Date(dob1$Date.of.birth, format = "%d.%m.%Y")
if (any(duplicated(dob1$hid))) {
  print("dob1: duplicated hid")
  # dob1 <- dob1[nrow(dob1):1, ]
  # dob1 <- dob1[!duplicated(dob1$hid), ]
}
rm(i)

# Date of birth (from last file of first survey)
load("data/dataFSOsero_24092020.RData")
tmp <- dataFSOsero[c("uc_info_participant_hid", "Date.of.birth")]
tmp$Date.of.birth <- as.Date(tmp$Date.of.birth, format = "%d.%m.%Y")
names(tmp) <- c("hid", "Date.of.birth.1")
if (any(duplicated(tmp$hid))) stop("Duplicated hid")
dob1 <- merge(dob1, tmp, by = "hid", all = TRUE, sort = FALSE)
rm(dataFSOsero, tmp)

# Date of birth (from last file of second survey)
f <- "data-raw/Intermed_results_Q2017_20201810_date_prélèvmt.xlsx"
tmp <- as.data.frame(read_xlsx(f, sheet = "màj_18.12.2020", guess_max = 10^4))
tmp <- unique(tmp[c("hid", "DDN")])
if (any(duplicated(tmp$hid))) stop("Duplicated hid")
i <- grepl("^10[0-9]{2}-[0-9]{2}-[09]{2}$", tmp$DDN)
delta <- as.numeric(as.Date("1970-01-01") - as.Date("1899-12-30"))
tmp$DDN[i] <- as.character(as.numeric(as.Date(sub("^10", "19", tmp$DDN[i])))
                             + delta)
tmp$DDN <- as.Date(as.numeric(tmp$DDN), origin = "1899-12-30")
names(tmp)[2] <- "Date.of.birth.2"
dob1 <- merge(dob1, tmp, by = "hid", all = TRUE, sort = FALSE)
rm(f, tmp, i, delta)

# Date of birth - birth date posterior to survey
i <- is.na(dob1$Date.of.birth) | dob1$Date.of.birth < as.Date("2020-01-01")
dob1[!i, ]
dob1 <- dob1[i, ]
rm(i)

# Date of birth - duplicated hid
dup <- dob1[dob1$hid %in% dob1$hid[duplicated(dob1$hid)], ]
dup <- dup[order(dup$hid), ]
dup
dob1 <- dob1[!(dob1$hid %in% dob1$hid[duplicated(dob1$hid)]), ]
rm(dup)

# Date of birth - differences
dob1[is.na(dob1$Date.of.birth), ]
subset(dob1, !is.na(Date.of.birth) & !is.na(Date.of.birth.1) &
               Date.of.birth != Date.of.birth.1)
subset(dob1, !is.na(Date.of.birth) & !is.na(Date.of.birth.2) &
               Date.of.birth != Date.of.birth.2)

# Date of birth - Strata
smpl <- read.csv("~/Projects/SerocoViD/data-fso/COVID19_VD_V1_Total.csv",
                 sep = ";")
smpl$dateOfBirth <- as.Date(smpl$dateOfBirth)
tmp <- unique(smpl[!is.na(smpl$strate), c("dateOfBirth", "strate")])
tmp <- do.call(rbind, lapply(sort(unique(tmp$strate)), function(i) {
  z <- tmp[tmp$strate == i, "dateOfBirth"]
  data.frame(dateOfBirth = seq(min(z), max(z), by = 1), stratum.0 = i)
}))
dob1 <- merge(dob1, tmp, by.x = "Date.of.birth", by.y = "dateOfBirth",
              all.x = TRUE, sort = FALSE)
names(tmp)[names(tmp) == "stratum.0"] <- "stratum.1"
dob1 <- merge(dob1, tmp, by.x = "Date.of.birth.1", by.y = "dateOfBirth",
              all.x = TRUE, sort = FALSE)
names(tmp)[names(tmp) == "stratum.1"] <- "stratum.2"
dob1 <- merge(dob1, tmp, by.x = "Date.of.birth.2", by.y = "dateOfBirth",
              all.x = TRUE, sort = FALSE)
rm(smpl, tmp)

# Date of birth - differences in stratum
dob1[is.na(dob1$stratum.0), ]
subset(dob1, !is.na(stratum.0) & !is.na(stratum.1) &
               stratum.0 != stratum.1)
subset(dob1, !is.na(stratum.0) & !is.na(stratum.2) &
               stratum.0 != stratum.2)

# Date of birth - final stratum
# hid 647737 : On retient la strate correspondant à la date de naissance
#              1983-03-31. La date de naissance 2019-03-31 ne pouvait pas
#              être incluse dans l'échantillon (condition : âge >= 6 mois)
# hid 121212 : On retient la strate correspondant à la date de naissance
#              1993-02-21 car c'est la seule disponnible
dob1$stratum <- with(dob1, ifelse(hid == 121212, stratum.2, stratum.0))
if (any(is.na(dob1$stratum))) stop("missing stratum")
if (any(duplicated(dob1$hid))) stop("duplicated hid")

# ------------------------ Research data - Survey 1 ------------------------- #

# Select data
vars <- c("uc_info_participant_hid", "uc_info_type_participant",
          "uc_labo_covigl", "uc_labo_covigl_2", "uc_labo_covigal",
          "uc_labo_covigal_2", "uc_age_in_months")
dta1 <- serocovid_data$research_data[vars]
dta1 <- dta1[dta1$uc_info_type_participant %in% 4, ]
rm(vars)

# Rename id variable
names(dta1)[names(dta1) == "uc_info_participant_hid"] <- "hid"

# Look for duplicates
if (any(duplicated(dta1$uc_info_participant_hid))) {
  stop("dta1: duplicated hid")
}

# Serology
dta1$serol_igg <-
  apply(dta1[c("uc_labo_covigl", "uc_labo_covigl_2")], 1, function(x) {
    if (any(x %in% "#p")) {
      r <- 1
    } else if (any(x %in% c("#n", "#l"))) {
      r <- 0
    } else {
      r <- NA
    }
  })
dta1$serol_iga <-
  apply(dta1[c("uc_labo_covigal", "uc_labo_covigal_2")], 1, function(x) {
    if (any(x %in% "#p")) {
      r <- 1
    } else if (any(x %in% c("#n", "#l"))) {
      r <- 0
    } else {
      r <- NA
    }
  })
dta1$serol_any <- pmin(dta1$serol_igg + dta1$serol_iga, 1)

# Select the observation for which serology is available
with(dta1, table(!is.na(serol_igg), !is.na(serol_iga)))
dta1 <- dta1[!is.na(dta1$serol_igg) | !is.na(dta1$serol_iga), ]

# Add strata
dta1 <- merge(dta1, dob1[c("hid", "stratum")], by = "hid", all.x = TRUE)
if (any(is.na(dta1$stratum))) stop("missing stratum")

# ----------------------------- Population size ----------------------------- #

# Population size - survey 1 & 2
pop1 <- read.table(header = TRUE, text = "
  stratum       N
        1   37490
        2   42538
        3   42700
        4   42171
        5  212336
        6  268898
        7   67412
        8   63097
")

# ------------------------------- Estimation -------------------------------- #

# Estimation - Function
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

# Estimate for people aged 20 to 69
dta1$All <- "All"
dta1$AgeClass <- cut(dta1$stratum, c(1, 4, 6, 8), include.lowest = TRUE,
                     labels = c("[6m-19]", "[20-64]", "[65+]"))
prop <- Mean(data = dta1, variable = "serol_any", stratum = "stratum",
             domain = c("All", "AgeClass"), pop = pop1)
prop$v <- NULL
names(prop)[names(prop) == "y"] <- "ppos"
write_xlsx(prop, "results/prev_20_64_surv1_20210628.xlsx")

