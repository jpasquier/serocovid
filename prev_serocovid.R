library(REDCapR)
library(readxl)
library(writexl)
library(ggplot2)
library(ggforce)
library(wesanderson)
library(RColorBrewer)

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

# ------------------- Date of birth / Stratum - Survey 3 -------------------- #

# Date of birth
dob3 <- read_xlsx("data-fso/COVID19_VD_V2_REDCap_et_Publipostage.xlsx",
                 sheet = "Feuil1", guess_max = 10^4)
dob3 <- as.data.frame(dob3[c("HID", "dateOfBirth")])
names(dob3)[names(dob3) == "HID"] <- "hid"
dob3$dateOfBirth <- as.Date(dob3$dateOfBirth)
if (any(duplicated(dob3$hid))) stop("duplicated hid")
if (any(is.na(dob3$dateOfBirth))) stop("missing date of birth")

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
dob3 <- merge(dob3, tmp, by = "dateOfBirth", all.x = TRUE)
if (any(is.na(dob3$stratum))) stop("missing stratum")
if (any(duplicated(dob3$hid))) stop("duplicated hid")
rm(smpl, tmp)

# ------------------------ Research data - Survey 1 ------------------------- #

# Select data
vars <- c("uc_info_participant_hid", "uc_info_type_participant",
          "uc_labo_covigl", "uc_labo_covigl_2", "uc_labo_covigal",
          "uc_labo_covigal_2")
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

# ---------------------- Research data - Survey 2 & 3 ----------------------- #

# Personal data - type of participant and date of birth
vars <- c("uc_s_participant_hid", "uc_s_type_participant", "full_dob")
pdta <- unique(serocovid_data$personal_data[vars])
if (any(duplicated(pdta$uc_s_participant_hid))) {
  stop("Duplicated hid")
}
vars <- c("hid", "uc_labo_coviggl_v2", "uc_labo_covigal_v2", "bl_vac_yn")
dta <- serocovid_data$corona_immunitas[vars]
dta <- merge(dta, pdta, by.x = "hid", by.y = "uc_s_participant_hid",
               all.x = TRUE)
if (any(is.na(dta$uc_s_type_participant))) {
  stop("missing type of participant")
}
rm(vars, pdta)

# Select observation of survey 2 & 3
dta <- dta[dta$uc_s_type_participant %in% c(4, 13), ]

# Remove fake participants from survey 3
if (FALSE) table(grepl("^123123", dta$hid), dta$uc_s_type_participant == 4)
dta <- dta[!(grepl("^123123", dta$hid) & dta$uc_s_type_participant == 13), ]

# Serology
if (FALSE) with(dta, table(uc_labo_coviggl_v2, uc_labo_covigal_v2))
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
with(dta, addmargins(table(serol_igg, serol_iga,
                           uc_s_type_participant, useNA = "ifany")))

# Select survey 2 and add strata
dta2 <- dta[dta$uc_s_type_participant == 4, ]
dta2 <- merge(dta2, dob1[c("hid", "stratum")], by = "hid", all.x = TRUE)
if (any(is.na(dta2$stratum))) stop("missing stratum")

# Select survey 3, add strata and renumber this ones
dta3 <- dta[dta$uc_s_type_participant == 13, ]
dta3 <- merge(dta3, dob3[c("hid", "stratum")], by = "hid", all.x = TRUE)
dta3$stratum <- dta3$stratum + 3
if (any(is.na(dta3$stratum))) stop("missing stratum")
rm(dta)

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

# Population totals - survey 3
pop3 <- read.table(header = TRUE, text = "
  stratum       N
        4   42140
        5  214498
        6  270755
        7   67462
        8   63605
")

# Vaccinated persons - Survey 3
vac3 <-  read.table(header = TRUE, text = "
  stratum       N
        4     162
        5    4328
        6   10132
        7    5609
        8   17851
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

# Estimation by survey
# Survey 0 is survey 1 with all strata
# Survey 4 is survey 3 without vaccinated people
dta1$all <- 0
dta2$all <- 0
dta3$all <- 0
long <- do.call(rbind, lapply(1:3, function(j) {
  v <- c("serol_any", "serol_igg", "serol_iga")[j]
  w <- c("IgG or IgA", "IgG", "IgA")[j]
  do.call(rbind, lapply(0:4, function(k) {
    if (k <= 1) {
      dta <- dta1
    } else if (k == 2) {
      dta <- dta2
    } else if (k == 3) {
      dta <- dta3
    } else {
      dta <- dta3[dta3$bl_vac_yn %in% 2, ]
    }
    if (k >= 1) {
      dta <- subset(dta, stratum >= 4)
    }
    if (k <= 2) {
      pop <- pop1
    } else if (k == 3) {
      pop <- pop3
    } else {
      pop <- merge(pop3, vac3, by = "stratum", suffixes = c(".all", ".vac"))
      pop$N <- pop$N.all - pop$N.vac
    }
    p <- Mean(data = dta, variable = v, stratum = "stratum",
              domain = c("all", "stratum"), pop = pop)
     cbind(antibody = w, visit = k, p)
  }))
}))
long <- long[!(long$visit == 0 & long$value %in% c(0, 4:8)), ]
long[long$visit == 0, "visit"] <- 1
names(long)[names(long) == "y"] <- "ppos"
long$antibody <- factor(long$antibody, c("IgG or IgA", "IgG", "IgA"))
long$domain <- factor(long$domain, c("all", "stratum"),
                      c("All (15+)", "Age groups (years)"))
long$value <- factor(long$value, 0:8,
                     c("All", "6m-4", "5-9", "10-14", "15-19", "20-39",
                       "40-64", "65-74", ">=75"))
long$v <- NULL
long <- long[order(long$antibody, long$visit, long$value), ]
long_ex <- long
long_ex$visit[long_ex$visit == 4] <- "3 vacc-"
if (FALSE) write_xlsx(long_ex, "results/prev_serocovid_20210325.xlsx")
rm(long_ex)

# Figure
fig <- lapply(list(1:3, c(1:2, 4), 1:4, 3:4), function(z) {
  visits <- c(
    "May 3rd, 2020\nJuly, 7, 2020",
    "October 20, 2020\nDecember 12, 2020",
    "February 1st, 2021\nFebruary 6, 2021",
    "February 1st, 2021\nFebruary 6, 2021\nwithout vaccinated\npeople"
  )
  if (length(z) == 4) {
    cols <- c(wes_palette(n = 4, name = "Darjeeling1")[c(1:2, 4)], "grey55")
  } else {
    cols <- wes_palette(n = length(z), name = "Darjeeling1")
  }
  tmp <- subset(long, antibody == "IgG or IgA" & visit %in% z)
  tmp$visit <- factor(tmp$visit, z, visits[z])
  p <- ggplot(tmp, aes(x = value, y = ppos, color = visit)) +
    geom_point(position = position_dodge(width = 0.3)) +
    geom_errorbar(aes(ymin = pmax(0, lwr), ymax = upr), width = 0,
                  position = position_dodge(width = 0.3)) +
    scale_y_continuous(labels = scales::percent) +
    scale_color_manual(values = cols) +
    labs(x = "", y = "Prevalence", color = "Surveys") +
    facet_row(vars(domain), scales = "free_x", space = "free") +
    theme_bw() +
    theme(legend.text = element_text(size=rel(0.6)),
          legend.key.size = unit(0.8, "cm"))
  return(p)
})
for (i in 3) { # 1:length(fig)) {
  jpeg(paste0("results/prev_serocovid_fig", i, "_20210325.jpg"),
       height = 3600, width = 9000, res = 1024)
  print(fig[[i]])
  dev.off()
}
rm(i)

#
rm(serocovid_data)
if (FALSE) save.image("results/prev_serocovid_20210325.dta", compress = "xz")
