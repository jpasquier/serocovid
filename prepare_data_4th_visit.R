library(RCurl)
library(XML)
library(parallel)

options(mc.cores = detectCores() - 1)

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

# -------------------------------- Serology --------------------------------- #

# Select variables
V_serol <- grep("covigg.+v4$", names(serocovid_data$corona_immunitas),
                value = TRUE)
V_serol <- c("hid", V_serol, "bl_vac_yn_ph4", "uc_labo_qc_v4")
pdta <- unique(serocovid_data$personal_data[
  c("uc_s_participant_hid", "uc_s_type_participant")])
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

# Serology - level
u0 <- "uc_labo_qc_v4"
u1 <- "uc_labo_coviggl_spe_v4"
if (any(is.na(serol[[u0]]) | is.na(serol[[u1]]))) stop("missing values")
for (i in 1:3) {
  v <- paste0("serol_", c("low", "medium", "high")[i])
  w <- c("faible", "moyen", "élevé")[i]
  W <- c("négatif", "faible", "moyen", "élevé")
  W <- W[W != w]
  serol[[v]] <- 9
  serol[serol[u0] != "OK" | serol[[u1]] == "", v] <- NA
  serol[serol[u0] == "OK" & serol[[u1]] %in% W, v] <- 0
  serol[serol[u0] == "OK" & serol[[u1]] == w, v] <- 1
  if (any(!is.na(serol[[v]]) & serol[[v]] == 9)) {
    stop("error in serology recoding")
  }
}
rm(i, u0, u1, v, w, W)
with(serol, table(serol, serol_low + serol_medium * 10 + serol_high * 100,
                  useNA = "ifany"))

# Serology - z score
if (any(is.na(serol$uc_labo_coviggr_zscore_v4) & !is.na(serol$serol))) {
  stop("Missing serology z score")
}
serol$serol_num <-
  with(serol, ifelse(!is.na(serol), uc_labo_coviggr_zscore_v4, NA))

# ----------------------- FSO sample (original file) ------------------------ #

# Import sample and select variables
V_fso <- c("PersonId", "HouseholdSizeSRPH", "sex", "maritalStatus",
           "nationalityState", "strate")
smpl <- read.csv("data-fso/COVID19_VD_V4_Total.csv", sep = ";",
                  fileEncoding = "latin1")
smpl <- smpl[!is.na(smpl$SAMPLE), V_fso]
rm(V_fso)

# Key between FSO and REDCap data
key <- readRDS("data/key_fso_redcap_v4.rds")

# Add hid variable and remove PersonId
smpl <- merge(smpl, key, by = "PersonId", all.x = TRUE)
if (any(is.na(smpl$hid))) stop("missing hid")
if (any(duplicated(smpl$PersonId)) | any(duplicated(smpl$hid))) {
  stop("Duplicated id")
}
smpl$PersonId <- NULL

# Add serology data
smpl <- merge(smpl, serol, by = "hid", all.x = TRUE)

# Recode FSO variables
names(smpl)[names(smpl) == "maritalStatus"] <- "maritalstatus"
names(smpl)[names(smpl) == "nationalityState"] <- "nationality"
names(smpl)[names(smpl) == "HouseholdSizeSRPH"] <- "householdsize"
smpl$strate <- factor(smpl$strate, levels = 1:8)
smpl$sex <- factor(smpl$sex, levels = 1:2, labels = c("Male", "Female"))
smpl$maritalstatus <- factor(smpl$maritalstatus, levels = 1:7)
levels(smpl$maritalstatus) <- c("Single", "Married", "Widowed", "Divorced",
                                 rep("Others", 3))
smpl$nationality <- factor(ifelse(smpl$nationality == 8100, 1, 2),
                           levels = 1:2, labels = c("Swiss", "Foreign"))
smpl$householdsize <- ifelse(smpl$householdsize >= 5, 5, smpl$householdsize)
smpl$householdsize <- factor(smpl$householdsize, levels = 1:5,
                             labels = c("1person", paste0(2:4, "people"),
                                        "5ormore"))

# Vaccination and Strate x Vaccination variables
smpl <- within(smpl, {
  vac <- c("y", "n", NA)[bl_vac_yn_ph4]
  stratevac <- ifelse(is.na(vac), NA, paste0(strate, vac))
})
smpl$stratevac <- factor(smpl$stratevac,
                         as.vector(outer(1:8, c("n", "y"), paste0)))
if(FALSE) with(smpl, table(stratevac, bl_vac_yn_ph4, strate, useNA = "ifany"))

# Respondents
smpl$resp <- as.numeric(!is.na(smpl$serol))
smpl$respvac <- as.numeric(!is.na(smpl$serol) & !is.na(smpl$vac))
if(FALSE) with(smpl, table(resp, respvac, useNA = "ifany"))

# ---------------------------- Sampling weights ----------------------------- #

# Population totals (FSO)
N <- read.table(header = TRUE, text = "
   strate       N
        1   38206
        2   43179
        3   43666
        4   42434
        5  216839
        6  273416
        7   67019
        8   64770
")

# Gross sample totals
n <- aggregate(hid ~ strate, smpl, length)
names(n)[2] <- "n"
N <- merge(N, n, by = "strate")
rm(n)

# Sampling weigths
N$wsmpl <- N$N / N$n
smpl <- merge(smpl, N, by = "strate")

# --------------------------------- Margins --------------------------------- #

# Vaccination rates (e-mail VF 8 oct 2021)
vac_rate <- read.table(header = TRUE, text = "
   strate    rate
        1       0
        2  0.0001
        3  0.1943
        4  0.6739
        5  0.7451
        6  0.8244
        7  0.8944
        8  0.9145
")

# Vaccination margins - table
z <- within(merge(N, vac_rate, by = "strate"), {
  yes <- round(N * rate)
  no <- round(N * (1 - rate))
})
vac_margins <- data.frame(
  variable = "vac",
  value = c("n", "y"),
  N = apply(z[c("no", "yes")], 2, sum)
)
stratevac_margins <- cbind(variable = "stratevac", rbind(
  data.frame(value = paste0(z$strate, "n"), N = z$no),
  data.frame(value = paste0(z$strate, "y"), N = z$yes)
))
rm(z)

# Margins - table (FSO)
margins <- read.table(header = TRUE, text = "
        variable       value        N
             sex        Male   387686
             sex      Female   401843
     nationality       Swiss   537298
     nationality     Foreign   252231
   maritalstatus      Single   375217
   maritalstatus     Married   309291
   maritalstatus     Widowed    31012
   maritalstatus    Divorced    71583
   maritalstatus      Others     2426
   householdsize     1person   134990
   householdsize     2people   214269
   householdsize     3people   152734
   householdsize     4people   186887
   householdsize     5ormore   100649
")
s <- cbind(variable = "strate", setNames(N[c("strate", "N")], c("value", "N")))
margins <- rbind(s, margins, vac_margins, stratevac_margins)
rm(s)

# Margins - vector
tot <- unique(aggregate(N ~ variable, margins, sum)$N)
margins <- do.call(base::c, lapply(unique(margins$variable), function(v) {
  tmp <- margins[margins$variable == v, ][-1, ]
  setNames(tmp$N, paste0(v, tmp$value))
}))
margins <- c(`(Intercept)` = tot, margins)
rm(tot)
