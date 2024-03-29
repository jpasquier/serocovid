library(RColorBrewer)
library(RCurl)
library(ggforce)
library(ggplot2)
library(parallel)
library(readxl)
library(wesanderson)
library(writexl)

# Working directory
setwd("~/Projects/SerocoViD")

# Today's date
today <- format(Sys.Date(), '%Y%m%d')

# ------------------------- REDCap data and strata -------------------------- #

# API
uri   <- "https://redcap.unisante.ch/api/"

# Tokens
tokens <- c("corona_immunitas", "personal_data", "research_data",
            "corona_immunitas_vague_6")
tokens <- lapply(setNames(tokens, tokens), function(z) {
  z <- paste0("misc/redcap_", z, ".token")
  readChar(z, file.info(z)$size - 1)
})

# Import
tmp_file <- "/mnt/ramdisk/serocovid_data.rda"
if (file.exists(tmp_file)) {
  load(tmp_file)
} else {
  serocovid_data <- mclapply(tokens, function(token) {
    read.csv(text = postForm(
      uri = uri,
      token = token,
      content = 'record',
      format = 'csv'
    ))
  })
  save(serocovid_data, file = tmp_file)
}
rm(uri, tokens, tmp_file)

# Strata
strata_v1 <- readRDS("data/strata_v1.rds")
strata_v3 <- readRDS("data/strata_v3.rds")
strata_v4 <- readRDS("data/strata_v4.rds")

# --------------------------- Serology - Survey 1 --------------------------- #

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
dta1 <- merge(dta1, strata_v1, by = "hid", all.x = TRUE)
if (any(is.na(dta1$stratum))) stop("missing stratum")

# ------------------------- Serology - Survey 2 & 3 ------------------------- #

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
dta2 <- merge(dta2, strata_v1, by = "hid", all.x = TRUE)
if (any(is.na(dta2$stratum))) stop("missing stratum")

# Select survey 3, add strata and renumber this ones
dta3 <- dta[dta$uc_s_type_participant == 13, ]
dta3 <- merge(dta3, strata_v3, by = "hid", all.x = TRUE)
if (any(is.na(dta3$stratum))) stop("missing stratum")
rm(dta)

# --------------------------- Serology - Survey 4 --------------------------- #

# Personal data - type of participant and date of birth
vars <- c("uc_s_participant_hid", "uc_s_type_participant")
pdta4 <- unique(serocovid_data$personal_data[vars])
if (any(duplicated(pdta4$uc_s_participant_hid))) {
  stop("Duplicated hid")
}
vars <- c("hid", "uc_labo_qc_v4", "uc_labo_coviggl_v4", "bl_vac_yn_ph4")
dta4 <- serocovid_data$corona_immunitas[vars]
dta4 <- merge(dta4, pdta4, by.x = "hid", by.y = "uc_s_participant_hid",
              all.x = TRUE)
if (any(is.na(dta4$uc_s_type_participant))) {
  stop("missing type of participant")
}
rm(vars, pdta4)

# Select observation of survey 4
dta4 <- dta4[dta4$uc_s_type_participant %in% 14, ]

# Remove fake participants from survey 4
dta4 <- dta4[!(dta4$hid %in% c("13131313", "987654321")), ]

# Serology
u0 <- "uc_labo_qc_v4"
u1 <- "uc_labo_coviggl_v4"
if (any(is.na(dta4[[u0]]) | is.na(dta4[[u1]]))) stop("missing values")
dta4$serol_igg <- 9
dta4[dta4[[u0]] != "OK" | dta4[[u1]] == "", "serol_igg"] <- NA
dta4[dta4[[u0]] == "OK" & dta4[[u1]] == "négatif", "serol_igg"] <- 0
dta4[dta4[[u0]] == "OK" & dta4[[u1]] == "positif", "serol_igg"] <- 1
if (any(!is.na(dta4$serol_igg) & dta4$serol_igg == 9)) {
  stop("error in serology recoding")
}
rm(u0, u1)

# Select the observation for which serology is available
dta4 <- dta4[!is.na(dta4$serol_igg), ]

# Add strata
dta4 <- merge(dta4, strata_v4, by = "hid", all.x = TRUE)
if (any(is.na(dta4$stratum))) stop("missing stratum")

# --------------------------- Serology - Survey 6 --------------------------- #

# Personal data - type of participant and strata
vars <- c("uc_s_participant_hid", "uc_s_type_participant", "uc_s_strate_6",
          "uc_s_age")
pdta6 <- unique(serocovid_data$personal_data[vars])
if (any(duplicated(pdta6$uc_s_participant_hid))) {
  stop("Duplicated hid")
}
vars <- c("hid", "spike_igg_qual", "nuc_igg_qual", "neutralisation_omicron",
          "bl_vac_yn")
dta6 <- serocovid_data$corona_immunitas_vague_6[vars]
rm_hid <- c("6-841282", "6-76BKAW", "6-BHSMPF", paste0("test", 1:6))
dta6 <- dta6[!(dta6$hid %in% rm_hid), ]
dta6 <- merge(dta6, pdta6, by.x = "hid", by.y = "uc_s_participant_hid",
              all.x = TRUE)
if (any(is.na(dta6$uc_s_type_participant))) {
  stop("missing type of participant")
}
if (any(is.na(dta6$uc_s_strate_6))) {
  stop("missing strata")
}
rm(vars, pdta6)

# Select observation of survey 4
dta6 <- dta6[dta6$uc_s_type_participant %in% 15, ]

# Compare spike_igg and nuc_igg
with(dta6, table(spike_igg_qual, nuc_igg_qual, useNA = "ifany"))
with(dta6, table(dta6$spike_igg_qual == 1 | dta6$nuc_igg_qual == 1,
                 neutralisation_omicron >= 50, useNA = "ifany"))

# Serology
dta6$serol_igg <- as.numeric(dta6$spike_igg_qual == 1 | dta6$nuc_igg_qual == 1)
dta6$serol_igg_nuc <- dta6$nuc_igg_qual %% 2
dta6$serol_omicron <- dta6$neutralisation_omicron >= 50

# Rename stratum variable
names(dta6)[names(dta6) == "uc_s_strate_6"] <- "stratum"

# Age groups corresponding to the preceding serocovid strata
dta6$age_group <- cut(dta6$uc_s_age, c(0, 4, 9, 14, 19, 39, 64, 74, Inf),
                      c("6m-4", "5-9", "10-14", "15-19", "20-39", "40-64",
                        "65-74", ">=75"))
dta6$age_group <- droplevels(dta6$age_group)

# Select the observation for which serology is available
dta6 <- dta6[!is.na(dta6$serol_igg), ]

# ----------------------------- Population size ----------------------------- #

# Population size - survey 1, 3 and 4
pop <- readRDS("data/population_size_by_stratum.rds")

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
long <- do.call(rbind, mclapply(c(1:4, 6), function(j) {
  dta <- get(paste0("dta", j))
  if (j == 2) dta <- dta[dta$stratum != 3, ]
  dta$all <- factor("All")
  dta$p15 <- factor(as.numeric(dta$stratum >= (if (j <= 4) 4 else 1)),
                    0:1, c("<=14", ">=15"))
  if (j %in% c(1, 4)) {
    stratum_levels <- 1:8
    stratum_labels <- c("6m-4", "5-9", "10-14", "15-19", "20-39", "40-64",
                        "65-74", ">=75")
  } else if (j %in% 2:3) {
    stratum_levels <- 4:8
    stratum_labels <- c("15-19", "20-39", "40-64", "65-74", ">=75")
  } else {
    stratum_levels <- 1:4
    stratum_labels <- c("15-29", "30-44", "45-64", ">=65")
  }
  dta$stratum <- factor(dta$stratum, stratum_levels, stratum_labels)
  pop <- subset(pop, survey == c(1, 1, 3:4, NA, 6)[j])
  pop$stratum <- factor(pop$stratum, stratum_levels, stratum_labels)
  domains <- c("all", "p15", "stratum")
  if (j >= 3) {
    names(dta)[grep("^bl_vac_yn", names(dta))] <- "vac"
    if (j == 6) {
      dta$vac <- factor(dta$vac, 1:0, c("vac", "unvac"))
    } else {
      dta$vac <- factor(dta$vac, 1:2, c("vac", "unvac"))
    }
    if (j == 6) domains <- c(domains, "age_group")
    for (d in domains) {
      x <- paste0(d, "_vac")
      dta[[x]] <- interaction(dta[[d]], dta$vac)
      domains <- c(domains, x)
    }
  }
  K <- if (j <= 3) 1:3 else if (j == 4) 2 else c(2, 4:5)
  do.call(rbind, lapply(K, function(k) {
    v <- c("serol_any", "serol_igg", "serol_iga", "serol_igg_nuc",
           "serol_omicron")[k]
    w <- c("IgG or IgA", "IgG", "IgA", "IgG (NUC only)", "IgG (Omicron)")[k]
    p <- Mean(data = dta, variable = v, stratum = "stratum", domain = domains,
              pop = pop)
     cbind(antibody = w, visit = j, p)
  }))
}))
names(long)[names(long) == "y"] <- "ppos"
long$antibody <- factor(long$antibody, c("IgG or IgA", "IgG", "IgA",
                                         "IgG (NUC only)", "IgG (Omicron)"))
long$v <- NULL
if (TRUE) {
  write_xlsx(long, paste0("results/prev_serocovid_", today, ".xlsx"))
  saveRDS(long, paste0("results/prev_serocovid_", today, ".rds"),
          compress = "xz")
}

# Figure
fig <- lapply(list(1:4), function(z) {
  visits <- c(
    "May 3rd, 2020\nJuly, 7, 2020",
    "October 20, 2020\nDecember 12, 2020",
    "February 1st, 2021\nFebruary 6, 2021",
    "October, 2021"
  )
  if (length(z) == 4) {
    cols <- c(wes_palette(n = 4, name = "Darjeeling1")[c(1:2, 4)], "grey55")
  } else {
    cols <- wes_palette(n = length(z), name = "Darjeeling1")
  }
  tmp <- subset(long, antibody == "IgG" & visit %in% z)
  tmp <- subset(tmp, !(value %in% c("All", "<=14")))
  tmp <- subset(tmp, !(visit == 2 & value == "10-14"))
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
if (FALSE) {
  for (i in 1:length(fig)) {
    jpeg(paste0("results/prev_serocovid_fig", i, "_DEV.jpg"),
         height = 3600, width = 9000, res = 1024)
    print(fig[[i]])
    dev.off()
  }
  rm(i)
}

#
rm(serocovid_data)
if (TRUE) {
  save.image(paste0("results/prev_serocovid_", today, ".rda"), compress = "xz")
}

###############################################################################
# library(survey)
# wgt <- aggregate(hid ~ stratum, dta6, length)
# names(wgt)[2] <- "n"
# wgt <- merge(wgt, subset(pop, survey == 6), by = "stratum")
# wgt$wgt <- wgt$N / wgt$n
# tmp <- merge(dta6, wgt, by = "stratum")
# d <- svydesign(
#   id = ~hid,
#   strata = ~stratum,
#   weights = ~wgt,
#   data = tmp,
#   fpc = ~N
# )
# (prev <- svymean(~serol_igg, d))
# confint(prev)
# confint(prev, df = sum(wgt$n) - nrow(wgt)) # loi de Student au lieu de la loi
#                                            # normale
# (prev_by_stratum <- svyby(~serol_igg, ~stratum, d, svymean))
# confint(prev_by_stratum)
