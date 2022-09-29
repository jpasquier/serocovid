#    ____                              __     __ _  ____    _  _   
#   / ___|   ___  _ __  ___    ___  ___\ \   / /(_)|  _ \  | || |  
#   \___ \  / _ \| '__|/ _ \  / __|/ _ \\ \ / / | || | | | | || |_ 
#    ___) ||  __/| |  | (_) || (__| (_) |\ V /  | || |_| | |__   _|
#   |____/  \___||_|   \___/  \___|\___/  \_/   |_||____/     |_|  

library(REDCapR)
library(writexl)

# Set working directory
setwd("~/Projects/SerocoViD")

# Initialisation du générateur de nombres aléatoires
set.seed(9958)

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

# ------------------------ Research data - Survey 3 ------------------------- #

# Personal data - type of participant and date of birth
vars <- c("uc_s_participant_hid", "uc_s_type_participant")
pdta <- unique(serocovid_data$personal_data[vars])
if (any(duplicated(pdta$uc_s_participant_hid))) {
  stop("Duplicated hid")
}
vars <- c("hid", "uc_labo_coviggl_v2", "uc_labo_covigal_v2")
dta <- serocovid_data$corona_immunitas[vars]
dta <- merge(dta, pdta, by.x = "hid", by.y = "uc_s_participant_hid",
               all.x = TRUE)
if (any(is.na(dta$uc_s_type_participant))) {
  stop("missing type of participant")
}
rm(vars, pdta)

# Select participants of survey 3
dta <- dta[dta$uc_s_type_participant %in% 13 & !grepl("^123123", dta$hid), ]
if (any(is.na(dta$uc_labo_coviggl_v2) != is.na(dta$uc_labo_covigal_v2))) {
  stop("check serology")
}
dta <- dta[!is.na(dta$uc_labo_coviggl_v2), "hid", drop = FALSE]

# Add strata
dta <- merge(dta, read.csv("data-fso/link_fso_redcap_sample2.csv", sep = ";"),
             by = "hid", all.x = TRUE)
if (any(is.na(dta$PersonId)) | any(duplicated(dta$PersonId))) {
  stop("missing or duplicated fso id")
}
fso <- read.csv("data-fso/COVID19_VD_V2_Total.csv", sep = ";")
dta <- merge(dta, fso[c("PersonId", "strate")], by = "PersonId", all.x = TRUE)
dta$strate_v4 <- dta$strate - 1
names(dta)[names(dta) == "strate"] <- "strate_v3"
rm(fso)

# ------------------------------- Sample size ------------------------------- #

# Taille de l'échantillon net par strate
smpl_size <- read.table(header = TRUE, text = "
   strate_v4    n
           1  158
           2  165
           3  165
           4  113
")

# Taille de l'échantillon brut
response_rate <- readRDS("data/response_rate_v2.rds")
response_rate$strate_v4 <- response_rate$stratum - 4
smpl_size <- merge(smpl_size, response_rate[c("strate_v4", "response_rate")],
                   by = "strate_v4")
smpl_size$N <- ceiling(smpl_size$n / smpl_size$response_rate)
nmax <- aggregate(hid ~ strate_v4, dta, length)
names(nmax)[names(nmax) == "hid"] <- "N_max"
smpl_size <- merge(smpl_size, nmax, by = "strate_v4")
smpl_size$N <- pmin(smpl_size$N, smpl_size$N_max)
#writexl::write_xlsx(smpl_size, "~/smpl_size_v4.xlsx")

# Sélection (par strate) des IDs à l'aide de la fonction sample
smpl <- do.call(rbind, lapply(smpl_size$strate_v4, function(s) {
  smpl <- sample(
    dta[dta$strate_v4 == s, "hid"],
    smpl_size[smpl_size$strate_v4 == s, "N"],
    replace = FALSE
  )
  data.frame(hid = smpl, selected_v4 = 1)
}))
dta <- merge(dta, smpl, by = "hid", all.x = TRUE)
dta[is.na(dta$selected_v4), "selected_v4"] <- 0

# Export selection
saveRDS(dta, "data/sample_v4.rds", compress = "xz")
write.table(dta, "data/sample_v4.csv", sep = ";", row.names = FALSE)
