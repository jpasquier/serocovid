library(REDCapR)
library(readxl)

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
smpl <- read.csv("data-fso/COVID19_VD_V1_Total.csv", sep = ";")
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

# Save strata
saveRDS(dob1[c("hid", "stratum")], "data/strata_serocovid_1_2.rds")
