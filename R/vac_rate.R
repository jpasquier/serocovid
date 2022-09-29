library(RCurl)
library(dplyr)
library(survey)
library(writexl)

# Working directory
setwd("~/Projects/SerocoViD")

# ------------------------- REDCap data and strata -------------------------- #

# API
uri   <- "https://redcap.unisante.ch/api/"

# Tokens
tokens <- c("personal_data", "corona_immunitas", "corona_immunitas_vague_6")
tokens <- lapply(setNames(tokens, tokens), function(z) {
  z <- paste0("misc/redcap_", z, ".token")
  readChar(z, file.info(z)$size - 1)
})

# Import
tmp_file <- "/mnt/ramdisk/serocovid_data.rda"
if (file.exists(tmp_file)) {
  load(tmp_file)
} else {
  serocovid_data <- lapply(tokens, function(token) {
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
strata_v3 <- readRDS("data/strata_v3.rds")
strata_v4 <- readRDS("data/strata_v4.rds")

# Population size
pop <- readRDS("data/population_size_by_stratum.rds") %>%
  mutate(survey = if_else(survey == 6L, 5L, survey))

# --------------------------- Vaccination status ---------------------------- #

age_grps <- 
  c("6m-4", "5-9", "10-14", "15-19", "20-39", "40-64", "65-74", ">=75")
dta <- bind_rows(
  select(serocovid_data$corona_immunitas, hid, bl_vac_yn, bl_vac_yn_ph4,
         uc_labo_coviggl_v2, uc_labo_qc_v4, uc_labo_coviggl_v4),
  select(serocovid_data$corona_immunitas_vague_6, hid, bl_vac_yn,
         spike_igg_qual, nuc_igg_qual)
) %>%
  left_join(
    serocovid_data$personal_data %>%
      unique() %>%
      select(hid = uc_s_participant_hid, type = uc_s_type_participant,
             stratum_v6 = uc_s_strate_6, age = uc_s_age),
    by = "hid"
  ) %>%
  left_join(bind_rows(strata_v3, strata_v4), by = "hid") %>%
  filter(type %in% 13:15 & !(
    type == 13 & grepl("^123123", hid) |
    type == 14 & hid %in% c("13131313", "987654321") |
    type == 15 & hid %in%
      c("6-841282", "6-76BKAW", "6-BHSMPF", paste0("test", 1:6))
  )) %>%
  mutate(
    vac = case_when(
      type == 13 ~ as.integer(2 - bl_vac_yn),
      type == 14 ~ as.integer(2 - bl_vac_yn_ph4),
      type == 15 ~ as.integer(bl_vac_yn),
    ),
    serol_igg = case_when(
      type == 13 ~ case_when(
        grepl("(Négatif|Ind(e|é)terminé)", uc_labo_coviggl_v2) ~ 0,
        grepl("Positif", uc_labo_coviggl_v2) ~ 1,
        TRUE ~ NA_real_
      ),
      type == 14 ~ case_when(
                     uc_labo_qc_v4 != "OK" ~ NA_real_,
                     uc_labo_coviggl_v4 == "négatif" ~ 0,
                     uc_labo_coviggl_v4 == "positif" ~ 1,
                   ),
      type == 15 ~ as.numeric(spike_igg_qual == 1 | nuc_igg_qual == 1)
    ),
    age_group = case_when(
      type == 15 ~ cut(age, c(0, 4, 9, 14, 19, 39, 64, 74, Inf),
                    c("6m-4", "5-9", "10-14", "15-19", "20-39", "40-64",
                      "65-74", ">=75")),
      TRUE ~ factor(stratum, 1:8, age_grps)
    )
  ) %>%
  filter(vac %in% 0:1, serol_igg %in% 0:1) %>%
  mutate(
    stratum = ifelse(type == 15, stratum_v6, stratum),
    survey = type - 10
  ) %>%
  left_join(select(pop, -age), by = c("survey", "stratum")) %>%
  select(hid, survey, stratum, age_group, vac, N) %>%
  group_by(survey, stratum) %>%
  mutate(n = n(), w = N / n) %>%
  ungroup()

# ------------------- Estimation with the survey package -------------------- #

N <- bind_rows(
  dta %>%
    group_by(survey, age_group) %>%
    summarise(n = n(), .groups = 'drop'),
  dta %>%
    group_by(survey) %>%
    summarise(n = n(), .groups = 'drop') %>%
    mutate(age_group = 'all'),
  dta %>%
    filter(survey == 4, stratum >= 4) %>%
    group_by(survey) %>%
    summarise(n = n(), .groups = 'drop') %>%
    mutate(age_group = '>=15')
)
vac_rate_tbl <- do.call(rbind, lapply(3:5, function(s) {
  d <- svydesign(
    id = ~hid,
    strata = ~stratum,
    weights = ~w,
    data = filter(dta, survey == s),
    fpc = ~N
  )
  cn <- c("age_group", "vac_rate", "ci_lwr", "ci_upr")
  r0 <- svymean(~vac, d)
  r0 <- setNames(cbind(data.frame("all", r0[1]), confint(r0)), cn)
  r1 <- svyby(~vac, ~age_group, d, svymean)
  r1 <- setNames(cbind(r1[, 1:2], confint(r1)), cn)
  if (s == 4) {
    r2 <- svyby(~vac, ~factor(I(stratum >= 4), c(FALSE, TRUE),
                              c("<15", ">=15")),
                d, svymean)
    r2 <- setNames(cbind(r2[, 1:2], confint(r2)), cn)[2, ]
  } else {
    r2 <- NULL
  }
  r <- cbind(survey = s, bind_rows(r0, r2, r1))
  rownames(r) <- NULL
  return(r)
})) %>%
  left_join(N, by = c("survey", "age_group"))
write_xlsx(vac_rate_tbl, paste0("results/vac_rate_tbl_",
                                format(Sys.time(), "%Y%m%d"), ".xlsx"))


