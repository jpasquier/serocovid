library(RCurl)
library(XML)
library(dplyr)
library(ggplot2)
library(gridExtra)

# Working directory
setwd("~/Projects/SerocoViD")

# ------------------------------- REDCap data ------------------------------- #

# API
uri   <- "https://redcap.unisante.ch/api/"

# Tokens
tokens <- c("corona_immunitas", "personal_data")
tokens <- lapply(setNames(tokens, tokens), function(z) {
  z <- paste0("misc/redcap_", z, ".token")
  readChar(z, file.info(z)$size - 1)
})

# Import data with the API (RCurl)
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

# -------------------------------- Serology --------------------------------- #

# Select variables
V_serol <- grep("covigg.+v4$", names(serocovid_data$corona_immunitas),
                value = TRUE)
V_serol <- c("hid", V_serol, "bl_vac_yn_ph4", "uc_labo_qc_v4")
pdta <- unique(serocovid_data$personal_data[
  c("uc_s_participant_hid", "uc_s_type_participant", "uc_s_strate_3")])
serol <- serocovid_data$corona_immunitas[V_serol]
serol <- merge(serol, pdta, by.x = "hid", by.y = "uc_s_participant_hid",
               all.x = TRUE)
if (any(is.na(serol$uc_s_type_participant))) {
  stop("missing type of participant")
}
rm(pdta, V_serol)

# Select observation of survey 4
serol <- serol[serol$uc_s_type_participant %in% 14, ]
serol$uc_s_type_participant <- NULL


# Remove fake participants
hid_fake_participants <- c("13131313", "987654321")
serol[serol$hid %in% hid_fake_participants, ]
serol <- serol[!(serol$hid %in% hid_fake_participants), ]
rm(hid_fake_participants)

# Serology
if (FALSE) {
  with(serol, table(uc_labo_coviggl_v4, bl_vac_yn_ph4, useNA = "ifany"))
  with(serol, table(uc_labo_coviggl_v4, uc_labo_qc_v4, useNA = "ifany"))
}
u0 <- "uc_labo_qc_v4"
u1 <- "uc_labo_coviggl_v4"
if (any(is.na(serol[[u0]]) | is.na(serol[[u1]]))) stop("missing values")
serol$serol <- 9
serol[serol[u0] != "OK" | serol[[u1]] == "", "serol"] <- NA
serol[serol[u0] == "OK" & serol[[u1]] == "négatif", "serol"] <- 0
serol[serol[u0] == "OK" & serol[[u1]] == "positif", "serol"] <- 1
if (any(!is.na(serol$serol) & serol$serol == 9)) {
  stop("error in serology recoding")
}
rm(u0, u1)

# Import uc_labo_coviggr_Zscore_v4 from xlsx file
if (FALSE) {
d1 <- unique(serocovid_data$personal_data[
  c("uc_s_participant_sid", "uc_s_participant_hid")])
names(d1) <- sub("uc_s_participant_", "", names(d1))
d2 <- readxl::read_xlsx(sheet = "ALL_CHIPS",
  "data-raw/20211130_results_all_plates_envoi_unisante_v2_Zscore.xlsx")
d2 <- d2 %>%
  select(sid = `SID participant`,
         uc_labo_coviggr_zscore_v4 = uc_labo_coviggr_Zscore_v4) %>%
  left_join(d1, by = "sid") %>%
  select(-sid)
serol <- serol %>%
  select(-uc_labo_coviggr_zscore_v4) %>%
  left_join(d2, by = "hid")
rm(d1, d2)
}

# Serology - Numeric value
if (any(is.na(serol$uc_labo_coviggr_v4) & !is.na(serol$serol))) stop()
serol$serol_num_unstandardized <-
  with(serol, ifelse(!is.na(serol), uc_labo_coviggr_v4, NA))
if (any(is.na(serol$uc_labo_coviggr_zscore_v4) & !is.na(serol$serol))) stop()
serol$serol_num <-
  with(serol, ifelse(!is.na(serol), uc_labo_coviggr_zscore_v4, NA))

# Compare original and standardized values
dta <- serol %>%
  select(serol, uc_labo_coviggl_spe_v4, uc_labo_coviggr_v4,
         uc_labo_coviggr_zscore_v4) %>%
  na.omit() %>%
  mutate(uc_labo_coviggl_spe_v4 =
    factor(uc_labo_coviggl_spe_v4, c("négatif", "faible", "moyen", "élevé")))
p1 <- ggplot(dta, aes(x = uc_labo_coviggl_spe_v4,
                      y = uc_labo_coviggr_v4)) +
  geom_boxplot() +
  labs(title = "Sans normalisation", caption = paste("N =", nrow(dta)))
p2 <- ggplot(dta, aes(x = uc_labo_coviggl_spe_v4,
                      y = uc_labo_coviggr_zscore_v4)) +
  geom_boxplot() +
  labs(title = "Avec normalistation", caption = paste("N =", nrow(dta)))
tiff("results/antibody_score_4th_visit_standardization.tiff",
     height = 3600, width = 7200, res = 768, compression = "zip")
grid.arrange(p1, p2, nrow = 1)
dev.off()
rm(dta, p1, p2)

# Boxplots
dta <- serol %>%
  select(uc_s_strate_3, serol, serol_num, bl_vac_yn_ph4) %>%
  rename(strate = uc_s_strate_3) %>%
  mutate(
    strate = factor(strate, 1:8, c("6m-4", "5-9", "10-14", "15-19", "20-39",
                                   "40-64", "65-74", ">=75")),
    vac = factor(bl_vac_yn_ph4, 1:2, c("Vaccinated", "Not vaccinated"))
  ) %>%
  na.omit() %>%
  filter(serol == 1)
tiff("results/antibody_score_4th_visit_by_vac.tiff",
     height = 3600, width = 5400, res = 1024, compression = "zip")
p1 <- ggplot(dta, aes(x = vac, y = serol_num)) +
  geom_boxplot() +
  labs(x = "", y = "IgG score", subtitle = "By vaccination status",
       title = "Antibody score", caption = paste("N =", nrow(dta)))
print(p1)
dev.off()
tiff("results/antibody_score_4th_visit_by_vac_and_age.tiff",
     height = 3600, width = 7200, res = 768, compression = "zip")
p2 <- ggplot(dta, aes(x = strate, y = serol_num)) +
  geom_boxplot() +
  facet_grid(rows = vars(vac), scales = "free") +
  labs(x = "", y = "IgG score", subtitle = "By vaccination status and age",
       title = "Antibody score", caption = paste("N =", nrow(dta)))
print(p2)
dev.off()
tiff("results/antibody_log_score_4th_visit_by_vac.tiff",
     height = 3600, width = 5400, res = 1024, compression = "zip")
p3 <- ggplot(dta, aes(x = vac, y = log(serol_num))) +
  geom_boxplot() +
  labs(x = "", y = "IgG log score", subtitle = "By vaccination status",
       title = "Antibody log score", caption = paste("N =", nrow(dta)))
print(p3)
dev.off()
tiff("results/antibody_log_score_4th_visit_by_vac_and_age.tiff",
     height = 3600, width = 7200, res = 768, compression = "zip")
p4 <- ggplot(dta, aes(x = strate, y = log(serol_num))) +
  geom_boxplot() +
  facet_grid(rows = vars(vac), scales = "free") +
  labs(x = "", y = "IgG log score", subtitle = "By vaccination status and age",
       title = "Antibody log score", caption = paste("N =", nrow(dta)))
print(p4)
dev.off()
rm(dta, p1, p2, p3, p4)


