library(parallel)
library(RCurl)

# Working directory
setwd("~/Projects/SerocoViD")

# ------------------------------- REDCap data ------------------------------- #

# API
uri   <- "https://redcap.unisante.ch/api/"

# Token
token <- paste0("misc/redcap_personal_data.token")
token <- readChar(token, file.info(token)$size - 1)

# Import personal data with the API (RCurl)
pdta3 <- read.csv(text = postForm(
      uri = uri,
      token = token,
      content = 'record',
      format = 'csv'
    ))
rm(uri, token)

# Select participants of survey 4
pdta3 <- subset(pdta3, uc_s_type_participant %in% 13)
if (any(duplicated(pdta3$uc_s_participant_sid))) stop("Duplicated sid")
if (any(duplicated(pdta3$uc_s_participant_hid))) stop("Duplicated hid")

# Remove fake participants
pdta3 <- pdta3[!grepl("^123123", pdta3$uc_s_participant_hid), ]

# Identity variables
id_pdta3 <- c("uc_s_nom", "uc_s_prenom", "full_dob")
if (any(duplicated(pdta3[id_pdta3]))) stop("Duplicated identity")

pdta3$full_dob <- as.Date(pdta3$full_dob, "%d.%m.%Y")

# ----------------------- FSO sample (original file) ------------------------ #

smpl3 <- read.csv("data-fso/COVID19_VD_V2_Total.csv", sep = ";",
                  fileEncoding = "latin1")
smpl3 <- smpl3[!is.na(smpl3$SAMPLE), ]

# Identity variables
id_smpl3 <- c("officialName", "firstName", "dateOfBirth")
if (any(duplicated(smpl3[id_smpl3]))) stop("Duplicated identity")

smpl3$dateOfBirth <- as.Date(smpl3$dateOfBirth)

# --------------------- key between FSO and REDCap data --------------------- #

key <- merge(
  x = smpl3[c(id_smpl3, "PersonId")],
  y = pdta3[c(id_pdta3, "uc_s_participant_hid")],
  by.x = id_smpl3,
  by.y = id_pdta3,
  all = TRUE,
  sort = FALSE
)[c("PersonId", "uc_s_participant_hid")]
names(key)[2] <- "hid"
if (any(is.na(key$PersonId) | is.na(key$hid))) stop("missing id(s)")
if (any(duplicated(key$PersonId)) | any(duplicated(key$hid))) {
  stop("duplicated id(s)")
}
if (nrow(key) != nrow(smpl3)) stop("missing participant(s)")

# Save key
saveRDS(key, "data/key_fso_redcap_v3.rds")
