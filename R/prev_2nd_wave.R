library(readxl)
library(writexl)

# Working directory
setwd("~/Projects/SerocoViD")

# Data
data_file <- "data-raw/Intermed_results_Q2017_20201810_date_prélèvmt.xlsx"
data <- as.data.frame(read_xlsx(data_file, sheet = "màj_18.12.2020",
                                guess_max = 10^4))
if (FALSE) {
  with(data, table(!is.na(uc_ser_date_hour_v2), !is.na(uc_labo_coviggl_v2)))
  with(data, table(!is.na(uc_ser_date_hour_v2), !is.na(uc_labo_covigal_v2)))
  subset(data, is.na(uc_ser_date_hour_v2) & !is.na(uc_labo_coviggl_v2))
}
sel <- subset(data, !is.na(uc_ser_date_hour_v2) | !is.na(uc_labo_coviggl_v2) |
                      !is.na(uc_labo_covigal_v2))
if (FALSE) data[!(data$hid %in% tmp$hid), ]
data <- sel
rm(sel)

# Duplicated IDs
if (any(duplicated(data$hid))) stop("duplicated hid")

# Serology
if (FALSE) with(data, table(uc_labo_coviggl_v2, uc_labo_covigal_v2))
for (x in c("uc_labo_coviggl_v2", "uc_labo_covigal_v2")) {
  i <- !is.na(data[[x]]) & grepl("(Négatif|Ind(e|é)terminé)", data[[x]])
  data[i, x] <- 0
  i <- !is.na(data[[x]]) & grepl("Positif", data[[x]])
  data[i, x] <- 1
  if (all(is.na(data[[x]]) | data[[x]] %in% c("0", "1"))) {
    data[[x]] <- as.numeric(data[[x]])
  } else {
    stop("error in serology recoding")
  }
}
data$serol_both <- pmin(data$uc_labo_coviggl_v2 + data$uc_labo_covigal_v2, 1)
rm(i, x)

# Week number
data$weeknbr <- as.integer(strftime(data$uc_ser_date_hour_v2, format = "%V"))
if (FALSE) {
  tmp <- as.data.frame(read_xlsx(data_file, sheet = "màj_10.12.2020"))
  tmp <- merge(data[c("hid", "weeknbr")], tmp[c("hid", "weeknbr")], by = "hid")
  all(tmp$weeknbr.y == tmp$weeknbr.x)
  rm(tmp)
}

# Select observations
if (FALSE) {
  with(data, table(
    weeknbr = !is.na(weeknbr),
    serol   = !is.na(uc_labo_coviggl_v2) | !is.na(uc_labo_covigal_v2)
  ))
  with(data, table(uc_labo_coviggl_v2, uc_labo_covigal_v2, useNA = "always"))
}
if (all(is.na(data$uc_labo_coviggl_v2) == is.na(data$uc_labo_covigal_v2))) {
  data <- data[!is.na(data$uc_labo_coviggl_v2), ]
} else {
  stop("missing serology")
}

# Date of birth
i <- grepl("^10[0-9]{2}-[0-9]{2}-[09]{2}$", data$DDN)
delta <- as.numeric(as.Date("1970-01-01") - as.Date("1899-12-30"))
data$DDN[i] <- as.character(as.numeric(as.Date(sub("^10", "19", data$DDN[i])))
                              + delta)
data$DDN <- as.Date(as.numeric(data$DDN), origin = "1899-12-30")
rm(i, delta)
if (any(is.na(data$DDN))) stop("Missing date of birth")

# Contrôle de l'âge
# Problème avec hid = 745932
#
#    hid        DDN age_in_months
# 745932 2020-05-12           719
#
if (FALSE) {
  f <- "SeroCOViDResearchDat-SrologiesIgGIgAPremi_DATA_2020-11-16_1054.csv"
  tmp <- read.csv(file.path("data-raw", f), sep = ";")
  tmp <- tmp[c("uc_info_participant_hid", "uc_age_in_months")]
  names(tmp) <- c("hid", "age_in_months")
  tmp <- merge(data, tmp, by = "hid", all = TRUE)[
    c("hid", "DDN", "age_in_months")]
  tmp <- tmp[!is.na(tmp$DDN) & !is.na(tmp$age_in_months), ]
  tmp[order(tmp$age_in_months), ]
  plot(age_in_months ~ DDN, tmp)
  tmp[tmp$hid == 745932, ]
  rm(tmp)
}
#data <- data[data$hid != 745932, ]

# Strata
smpl <- read.csv("data-fso/COVID19_VD_V1_Total.csv", sep = ";")
smpl$dateOfBirth <- as.Date(smpl$dateOfBirth)
tmp <- unique(smpl[!is.na(smpl$strate), c("dateOfBirth", "strate")])
tmp <- do.call(rbind, lapply(sort(unique(tmp$strate)), function(i) {
  z <- tmp[tmp$strate == i, "dateOfBirth"]
  data.frame(DDN = seq(min(z), max(z), by = 1), stratum = i)
}))
data <- merge(data, tmp, by = "DDN", all.x = TRUE, sort = FALSE)
data[data$hid == 745932, "stratum"] <- 6
if (any(duplicated(data$hid))) stop("duplicated hid")
if (any(is.na(data$stratum))) stop("missing stratum")
rm(tmp)

# Population size
pop <- read.table(header = TRUE, text = "
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

# Estimation
# Variance: Cochran, W. G., Sampling techniques John Wiley & Sons, 1977,
#           page 144, formula 5A.75
# domain: any variable

#variable <- "uc_labo_coviggl_v2"; stratum <- "stratum"; .pop <- pop
Mean <- function(data = data, variable = "serol", stratum = "stratum",
                 domain = NULL, .pop = pop) {
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
both <- c(npos = 343, tp = 339, nneg = 256, fp = 4)

# Stratum sizes
aggregate(hid ~ stratum, data, length)

# Grouped weeks
#if (!all(data$weeknbr %in% 43:50)) stop("Check weeknkr")
if (FALSE) table(data$weeknbr)
for (k in 0:6) {
  gpa <- 43:(43 + k)
  gpb <- (43:50)[!(43:50 %in% gpa)]
  Gpa <- paste(gpa, collapse = "-")
  Gpb <- paste(gpb, collapse = "-")
  data[[paste0("weekgp", k)]] <- ifelse(data$weeknbr %in% gpa, Gpa, Gpb)
}
rm(k, gpa, gpb, Gpa, Gpb)

# Prevalence
set.seed(666)
nsim <- 10^6
data$all <- ""
prev <- lapply(1:2, function(j) {
  # Stratum 1 has only one particpant, stratum 2 is empty, stratum 3 has
  # 4 participants
  if (j == 1) {
    data <- data[!(data$stratum %in% 1:2), ]
    pop <- pop[!(pop$stratum %in% 1:2), ]
  } else {
    data <- data[!(data$stratum %in% 1:3), ]
    pop <- pop[!(pop$stratum %in% 1:3), ]
  }
  do.call(rbind, lapply(1:3, function(k) {
    if (k == 1) {
      ab <- "IgG"
      vr <- "uc_labo_coviggl_v2"
      ig <- igg
    } else if (k == 2) {
      ab <- "IgA"
      vr <- "uc_labo_covigal_v2"
      ig <- iga
    } else {
      ab <- "IgG or IgA"
      vr <- "serol_both"
      ig <- both
    }
    se <- ig[["tp"]] / ig[["npos"]]
    sp <- 1 - ig[["fp"]] / ig[["nneg"]]
    v_se <- se * (1 - se) / ig[["npos"]]
    v_sp <- sp * (1 - sp) / ig[["nneg"]]
    dom <- c("all", "stratum", "weeknbr")
    prev <- Mean(data = data, variable = vr, domain = dom)
    W <- sort(unique(data$weeknbr))
    prev_per_week <- do.call(rbind, lapply(W, function(w) {
      wdata <- data[data$weeknbr %in% w, ]
      prev <- Mean(data = wdata, variable = vr)
      prev$domain <- "weeknbr (method 2)"
      prev$value <- w
      prev
    }))
    prev <- rbind(prev, prev_per_week)
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
})
names(prev) <- c("strates_3_8", "strates_4_8")

# Number of observations per week and per stratum
Nobs <- aggregate(uc_labo_coviggl_v2 ~ weeknbr + stratum, data, length)
Nobs <- Nobs[order(Nobs$weeknbr, Nobs$stratum), ]

# week 49
tmp <- as.data.frame(read_xlsx(data_file, sheet = "màj_10.12.2020"))
merge(tmp[tmp$weeknbr == 49, c("hid", "uc_labo_coviggl_v2")],
      data[data$weeknbr %in% 49, c("hid", "uc_labo_coviggl_v2", "stratum")],
      by = "hid", all = TRUE)

# Export results
sink("results/prev_2nd_wave_sessionInfo_20201218.txt")
print(sessionInfo(), locale = FALSE)
sink()
write_xlsx(c(prev, list(nobs = Nobs, data = data)),
           "results/prev_2nd_wave_20201218.xlsx")

# --------------------------------------------------------------------------- #
# --------------------------------------------------------------------------- #
# --------------------------------------------------------------------------- #

if (FALSE) {
  # With the survey package
  library(survey)
  wgt <- aggregate(hid ~ stratum, data, length)
  names(wgt)[2] <- "n"
  wgt <- merge(wgt, pop, by = "stratum")
  wgt$wgt <- wgt$N / wgt$n
  tmp <- merge(data, wgt, by = "stratum")
  d <- svydesign(
    id = ~hid,
    strata = ~stratum,
    weights = ~wgt,
    data = tmp,
    fpc = ~N
  )
  (prev <- svymean(~uc_labo_coviggl_v2, d))
  confint(prev)
  confint(prev, df = sum(wgt$n) - nrow(wgt)) # loi de Student au lieu de la loi normale
}
