library(REDCapR)

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

# ------------------------ Research data - Survey 2 ------------------------- #

# Personal data - type of participant and date of birth
vars <- c("uc_s_participant_hid", "uc_s_type_participant", "full_dob")
pdta2 <- unique(serocovid_data$personal_data[vars])
if (any(duplicated(pdta2$uc_s_participant_hid))) {
  stop("Duplicated hid")
}
vars <- c("hid", "uc_labo_coviggl_v2", "uc_labo_covigal_v2", "bl_vac_yn")
dta2 <- serocovid_data$corona_immunitas[vars]
dta2 <- merge(dta2, pdta2, by.x = "hid", by.y = "uc_s_participant_hid",
               all.x = TRUE)
if (any(is.na(dta2$uc_s_type_participant))) {
  stop("missing type of participant")
}
rm(vars, pdta2)

# Select observation of survey 2
dta2 <- dta2[dta2$uc_s_type_participant %in% 4, ]

# Serology
if (FALSE) with(dta2, table(uc_labo_coviggl_v2, uc_labo_covigal_v2))
for (k in 1:2) {
  u <- c("uc_labo_coviggl_v2", "uc_labo_covigal_v2")[k]
  v <- c("serol_igg", "serol_iga")[k]
  dta2[[v]] <- NA
  i <- !is.na(dta2[[u]]) & grepl("(Négatif|Ind(e|é)terminé)", dta2[[u]])
  dta2[i, v] <- 0
  i <- !is.na(dta2[[u]]) & grepl("Positif", dta2[[u]])
  dta2[i, v] <- 1
  if (all(is.na(dta2[[v]]) | dta2[[v]] %in% c("0", "1"))) {
    dta2[[v]] <- as.numeric(dta2[[v]])
  } else {
    stop("error in serology recoding")
  }
}
dta2$serol_any <- pmin(dta2$serol_igg + dta2$serol_iga, 1)
rm(k, i, u, v)

# Select the observation for which serology is available
dta2 <- dta2[!is.na(dta2$serol_igg) | !is.na(dta2$serol_iga), ]
with(dta2, addmargins(table(serol_igg, serol_iga,
                           uc_s_type_participant, useNA = "ifany")))

# --------------------------------------------------------------------------- #
# Taux de participation à SerocoViD 2 en fonction du résultat de SerocoViD 1  #
# --------------------------------------------------------------------------- #

V <- c("hid", "serol_igg", "serol_iga", "serol_any")
cmp <- merge(dta1[dta1$uc_age_in_months >= 210, ], dta2[V], by = "hid",
             all = TRUE, suffixes = c(".1", ".2"))
lapply(c("igg", "iga", "any"), function(x) {
  fml <- as.formula(paste0("!is.na(serol_igg.2) ~ serol_", x, ".1"))
  L <- lapply(list(length, sum), function(fct) aggregate(fml, cmp, fct))
  R <- merge(L[[1]], L[[2]], by = paste0("serol_", x, ".1"))
  names(R)[2:3] <- c("N", "n")
  R$p <- with(R, n / N)
  R$se <- with(R, sqrt(p * (1 - p) / n))
  prop_test <- with(R, prop.test(n, N))
  fisher_test <- fisher.test(cbind(R$n, R$N - R$n))
  return(list(R, prop_test, fisher_test))
})
