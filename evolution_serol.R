library(readxl)
library(writexl)
library(parallel)
library(survey)

options(mc.cores = parallel::detectCores())

# Working directory
setwd("~/Projects/SerocoViD")

# Data 1st wave
f_dta <- "SeroCOViDResearchDat-SrologiesIgGIgAPremi_DATA_2020-11-16_1054.csv"
data1 <- read.csv(file.path("data-raw", f_dta), sep = ";")
names(data1)[names(data1) == "uc_info_participant_hid"] <- "hid"
if (any(duplicated(data1$hid))) stop("ID not unique")
rm(f_dta)

# Serology 1st wave
data1$serol_igg_1 <-
  apply(data1[c("uc_labo_covigl", "uc_labo_covigl_2")], 1, function(x) {
    if (any(x %in% "#p")) {
      r <- 1
    } else if (any(x %in% c("#n", "#l"))) {
      r <- 0
    } else {
      r <- NA
    }
  })
data1$serol_iga_1 <-
  apply(data1[c("uc_labo_covigal", "uc_labo_covigal_2")], 1, function(x) {
    if (any(x %in% "#p")) {
      r <- 1
    } else if (any(x %in% c("#n", "#l"))) {
      r <- 0
    } else {
      r <- NA
    }
  })
if (FALSE) with(data1, table(serol_igg, serol_iga, useNA = "ifany"))

# Select variables 1st wave
X <- c("hid", "serol_igg_1", "serol_iga_1", "uc_age_in_months")
data1 <- data1[X]
rm(X)

# Data 2nd wave
f_dta <- "data-raw/Intermed_results_Q2017_20201123_date_prélèvmt.xlsx"
data2 <- read_xlsx(f_dta, sheet = "màj_23.11", range = "A1:I1217")
data2 <- as.data.frame(data2)
rm(f_dta)

# Serology 2nd wave
for (x in c("uc_labo_coviggl_v2", "uc_labo_covigal_v2")) {
  i <- !is.na(data2[[x]]) & grepl("(Négatif|Ind(e|é)terminé)", data2[[x]])
  data2[i, x] <- 0
  i <- !is.na(data2[[x]]) & grepl("Positif", data2[[x]])
  data2[i, x] <- 1
  if (all(is.na(data2[[x]]) | data2[[x]] %in% c("0", "1"))) {
    data2[[x]] <- as.numeric(data2[[x]])
  } else {
    stop("error in serology recoding")
  }
}
names(data2)[names(data2) == "uc_labo_coviggl_v2"] <- "serol_igg_2"
names(data2)[names(data2) == "uc_labo_covigal_v2"] <- "serol_iga_2"
rm(i, x)

# Select variables 2nd wave
tmp <- unique(data2[c("hid", "DDN")])
if (any(duplicated(tmp$hid))) stop("Duplicated date of birth")
i <- grepl("^10[0-9]{2}-[0-9]{2}-[09]{2}$", tmp$DDN)
delta <- as.numeric(as.Date("1970-01-01") - as.Date("1899-12-30"))
tmp$DDN[i] <- as.character(as.numeric(as.Date(sub("^10", "19", tmp$DDN[i])))
                             + delta)
tmp$DDN <- as.Date(as.numeric(tmp$DDN), origin = "1899-12-30")
data2 <- subset(data2, !is.na(serol_igg_2) | !is.na(serol_iga_2),
                c("hid", "serol_igg_2", "serol_iga_2"))
if (any(duplicated(data2$hid))) stop("Duplicated hid")
data2 <- merge(data2, tmp, by = "hid", all = TRUE, sort = FALSE)

# Merge waves 1 and 2
data <- merge(data1, data2, by = "hid")

# Date of birth
load("data/dataFSOsero_24092020.RData")
tmp <- dataFSOsero[c("uc_info_participant_hid", "Date.of.birth")]
tmp$Date.of.birth <- as.Date(tmp$Date.of.birth, format = "%d.%m.%Y")
names(tmp)[names(tmp) == "uc_info_participant_hid"] <- "hid"
data <- merge(data, tmp, by = "hid", all.x = TRUE, sort = FALSE)
if (FALSE) {
  subset(data, !is.na(DDN) & !is.na(Date.of.birth) & DDN != Date.of.birth)
}
rm(dataFSOsero, tmp)

# Duplicated IDs
if (any(duplicated(data$hid))) {
  dup_id <- function(.data, .id, .cols = NULL) {
    if (is.null(.cols)) .cols <- names(.data)
    data[data[[.id]] %in% .data[duplicated(.data[[.id]]), .id], .cols]
  }
  dup_id(data, "hid")
  stop("duplicated hid")
}

# Select observations with at least one serology result
i <- apply(!is.na(data[grep("^serol_ig(a|g)_(1|2)$", names(data))]), 1, any)
data <- data[i, ]
rm(i)

# Evolution variables
for (ig in c("igg", "iga")) {
  x1 <- paste("serol", ig, 1, sep = "_")
  x2 <- paste("serol", ig, 2, sep = "_")
  y <- paste("serol", ig, "evl", sep = "_")
  data[[y]] <- ifelse(
    is.na(data[[x1]]),
    ifelse(
      is.na(data[[x2]]),
      "--",
      ifelse(data[[x2]] == 0, "-0", "-1")
    ),
    ifelse(
      data[[x1]] == 0,
      ifelse(
        is.na(data[[x2]]),
        "0-",
        ifelse(data[[x2]] == 0, "00", "01")
      ),
      ifelse(
        is.na(data[[x2]]),
        "1-",
        ifelse(data[[x2]] == 0, "10", "11")
      )
    )
  )
}
rm(ig, x1, x2, y)
if (FALSE) {
  unique(data[c("serol_igg_1", "serol_igg_2", "serol_igg_evl")])
  unique(data[c("serol_iga_1", "serol_iga_2", "serol_iga_evl")])
  aggregate(hid ~ serol_igg_evl, data, length)
  aggregate(hid ~ serol_iga_evl, data, length)
}

# Strata
smpl <- read.csv("data-fso/COVID19_VD_V1_Total.csv", sep = ";")
smpl$dateOfBirth <- as.Date(smpl$dateOfBirth)
tmp <- unique(smpl[!is.na(smpl$strate), c("dateOfBirth", "strate")])
tmp <- do.call(rbind, lapply(sort(unique(tmp$strate)), function(i) {
  z <- tmp[tmp$strate == i, "dateOfBirth"]
  data.frame(dateOfBirth = seq(min(z), max(z), by = 1), stratum = i)
}))
data <- merge(data, tmp, by.x = "Date.of.birth", by.y = "dateOfBirth",
              all.x = TRUE, sort = FALSE)
data <- merge(data, tmp, by.x = "DDN", by.y = "dateOfBirth",
              all.x = TRUE, sort = FALSE, suffixes = c("", ".DDN"))
if (FALSE) {
  subset(data, !is.na(stratum) & !is.na(stratum.DDN) & stratum != stratum.DDN)
}
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

# Totals
tot <- lapply(c("serol_igg_evl", "serol_iga_evl"), function(x) {
  lapply(1:2, function(k) {
    data <- data[!is.na(data[[x]]), ]
    if (k == 1) {
      data <- data[!grepl("-", data[[x]]), ]
    }
    data$z <- 1
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
    sum(tmp$wgt)
    tot <- svyby(~z, as.formula(paste0("~", x)), d, svytotal)
    ci <- confint(tot)
    tot <- cbind(tot[, 1:2], confint(tot))
    names(tot)[names(tot) == "z"] <- "N"
    n <- aggregate(as.formula(paste("hid ~", x)), data, length)
    names(n)[names(n) == "hid"] <- "n"
    merge(n, tot, by = x, all = TRUE)
  })
})
tot <- unlist(tot, recursive = FALSE)
names(tot) <- c("IgG (1)", "IgG (2)", "IgA (1)", "IgA (2)")

# Export results
sink("results/evolution_serol_sessionInfo_20201124.txt")
print(sessionInfo(), locale = FALSE)
sink()
write_xlsx(tot, "results/evolution_serol_20201124.xlsx")
