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
strata_v3 <- readRDS("data/strata_v3.rds")
strata_v4 <- readRDS("data/strata_v4.rds")

# Population size
pop <- readRDS("data/population_size_by_stratum.rds")

# --------------------------- Vaccination status ---------------------------- #

age_grps <- data.frame(stratum = 1:8, age_group =
  c("6m-4", "5-9", "10-14", "15-19", "20-39", "40-64", "65-74", ">=75"))
age_grps_v6 <- data.frame(stratum_v6 = 1:4, age_group_v6 =
  c("15-29", "30-44", "45-64", ">=65"))
dta <- bind_rows(
  select(serocovid_data$corona_immunitas, hid, bl_vac_yn, bl_vac_yn_ph4),
  select(serocovid_data$corona_immunitas_vague_6, hid, bl_vac_yn)
) %>%
  left_join(
    serocovid_data$personal_data %>%
      unique() %>%
      select(hid = uc_s_participant_hid, type = uc_s_type_participant,
             stratum_v6 = uc_s_strate_6),
    by = "hid"
  ) %>%
  filter(type %in% 13:15 & !(
    type == 13 & grepl("^123123", hid) |
    type == 14 & hid %in% c("13131313", "987654321") |
    type == 15 & hid %in%
      c("6-841282", "6-76BKAW", "6-BHSMPF", paste0("test", 1:6))
  )) %>%
  mutate(vac = case_when(
    type == 13 ~ as.integer(2 - bl_vac_yn),
    type == 14 ~ as.integer(2 - bl_vac_yn_ph4),
    type == 15 ~ as.integer(bl_vac_yn),
  )) %>%
  filter(vac %in% 0:1) %>%
  left_join(bind_rows(strata_v3, strata_v4), by = "hid") %>%
  left_join(age_grps, by = "stratum") %>%
  left_join(age_grps_v6, by = "stratum_v6") %>%
  mutate(
    stratum = ifelse(type == 15, stratum_v6, stratum),
    age_group = ifelse(type == 15, age_group_v6, age_group),
    age_group =
      factor(age_group, c(age_grps$age_group, age_grps_v6$age_group_v6)),
    survey = ifelse(type == 15, 6, type - 10)
  ) %>%
  left_join(select(pop, -age), by = c("survey", "stratum")) %>%
  select(hid, survey, stratum, age_group, vac, N) %>%
  group_by(survey, stratum) %>%
  mutate(n = n(), w = N / n) %>%
  ungroup()

# ------------------- Estimation with the survey package -------------------- #

N <- dta %>%
  select(survey, age_group, N, n) %>%
  unique()
N <- N %>%
  group_by(survey) %>%
  summarise(age_group = "all", N = sum(N), n = sum(n)) %>%
  bind_rows(N)
vac_rate_tbl <- do.call(rbind, lapply(c(3:4, 6), function(s) {
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
  r <- cbind(survey = s, bind_rows(r0, r1))
  rownames(r) <- NULL
  return(r)
})) %>%
  left_join(N, by = c("survey", "age_group") ) 
write_xlsx(vac_rate_tbl, paste0("results/vac_rate_tbl_",
                                format(Sys.time(), "%Y%m%d"), ".xlsx"))


