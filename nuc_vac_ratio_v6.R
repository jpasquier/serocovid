library(RCurl)
library(dplyr)
library(parallel)
library(survey)
library(writexl)

# Working directory
setwd("~/Projects/SerocoViD")

# --------------------- REDCap data and population size --------------------- #

# API
uri   <- "https://redcap.unisante.ch/api/"

# Tokens
tokens <- c("personal_data", "corona_immunitas_vague_6")
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

# Population size
pop <- readRDS("data/population_size_by_stratum.rds")

# ------------------------------- Estimations ------------------------------- #

# Select data
dta <- serocovid_data$corona_immunitas %>%
  select(hid, spike_igg_qual, nuc_igg_qual, bl_vac_yn) %>%
  left_join(
    serocovid_data$personal_data %>%
      unique() %>%
      select(hid = uc_s_participant_hid, uc_s_type_participant,
             stratum = uc_s_strate_6),
    by = "hid"
  ) %>%
  filter(
    uc_s_type_participant == 15,
    !(hid %in% c("6-841282", "6-76BKAW", "6-BHSMPF", paste0("test", 1:6))),
    nuc_igg_qual %in% 0:1,
    bl_vac_yn %in% 0:1
  )

# Calculate weights
W <- dta %>%
  group_by(stratum) %>%
  summarise(
    n = n(),
    n_vac = sum(bl_vac_yn),
    n_nonvac = sum(1 - bl_vac_yn),
    n_nuc = sum(nuc_igg_qual),
    n_nuc_vac = sum(nuc_igg_qual * bl_vac_yn),
    n_nuc_nonvac = sum(nuc_igg_qual * (1 - bl_vac_yn)),
    .groups = "drop"
  ) %>%
  left_join(filter(pop, survey == 6), by = "stratum") %>%
  mutate(w = N / n)

# Survey design
d <- svydesign(
  id = ~hid,
  strata = ~stratum,
  weights = ~w,
  data = left_join(dta, W, by = "stratum"),
  fpc = ~N
)

# Estimations
N_pop <- do.call(rbind, mclapply(1:7, function(k) {
  if (k %in% 1:5) {
    fml <- as.formula(case_when(
      k == 1 ~ "~bl_vac_yn",
      k == 2 ~ "~I(1 - bl_vac_yn)",
      k == 3 ~ "~nuc_igg_qual",
      k == 4 ~ "~I(nuc_igg_qual * bl_vac_yn)",
      k == 5 ~ "~I(nuc_igg_qual * (1 - bl_vac_yn))"
    ))
    Z0 <- svytotal(fml, d)
    Z1 <- svyby(fml, ~stratum, d, svytotal)
  } else {
    if (k == 6) {
      num <- ~I(nuc_igg_qual * bl_vac_yn)
      den <- ~bl_vac_yn
    } else {
      num <- ~I(nuc_igg_qual * (1 - bl_vac_yn))
      den <- ~I(1 - bl_vac_yn)
    }
    Z0 <- svyratio(num, den, d)
    Z1 <- svyby(num, by = ~stratum, denominator = den, design = d, svyratio)
  }
  cn <- c("stratum", "estimate", "2.5 %", "97.5 %")
  Z0 <- setNames(cbind(data.frame(0, Z0[1]), confint(Z0)), cn)
  Z1 <- setNames(cbind(Z1[, 1:2], confint(Z1)), cn)
  v <- c("N_vac", "N_unvac", "N_nuc", "N_nuc_vac", "N_nuc_unvac", "R_nuc_vac",
         "R_nuc_unvac")[k]
  cbind(variable = c(v, rep(NA, 4)), bind_rows(Z0, Z1))
})) %>%
  mutate(stratum = 
    factor(stratum, 0:4, c("all", "15-29", "30-44", "45-64", ">=65"))) %>%
  rename(age_group = stratum)

# Export
list(sample = select(W, age_group = age, N, starts_with("n")),
     population = N_pop) %>%
  write_xlsx(paste0("results/nuc_vac_ratio_v6_",
                    format(Sys.time(), "%Y%m%d"), ".xlsx"))


