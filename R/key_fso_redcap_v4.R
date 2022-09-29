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
pdta4 <- read.csv(text = postForm(
      uri = uri,
      token = token,
      content = 'record',
      format = 'csv'
    ))
rm(uri, token)

# Select participants of survey 4
pdta4 <- subset(pdta4, uc_s_type_participant %in% 14)
if (any(duplicated(pdta4$uc_s_participant_sid))) stop("Duplicated sid")
if (any(duplicated(pdta4$uc_s_participant_hid))) stop("Duplicated hid")

# Remove fake participants
hid_fake_participants <- c("13131313", "987654321")
pdta4[pdta4$uc_s_participant_hid %in% hid_fake_participants,
      c("uc_s_nom", "uc_s_prenom")]
pdta4 <- pdta4[!(pdta4$uc_s_participant_hid %in% hid_fake_participants), ]

# Identity variables
id_pdta4 <- c("uc_s_nom", "uc_s_prenom", "full_dob")
if (any(duplicated(pdta4[id_pdta4]))) stop("Duplicated identity")

# ----------------------- FSO sample (original file) ------------------------ #

smpl4 <- read.csv("data-fso/COVID19_VD_V4_Total.csv", sep = ";",
                  fileEncoding = "latin1")
smpl4 <- smpl4[!is.na(smpl4$SAMPLE), ]

# Identity variables
id_smpl4 <- c("officialName", "firstName", "dateOfBirth")
if (any(duplicated(smpl4[id_smpl4]))) stop("Duplicated identity")

# --------------------- key between FSO and REDCap data --------------------- #

key <- merge(
  x = smpl4[c(id_smpl4, "PersonId")],
  y = pdta4[c(id_pdta4, "uc_s_participant_hid")],
  by.x = id_smpl4,
  by.y = id_pdta4,
  all = TRUE,
  sort = FALSE
)[c("PersonId", "uc_s_participant_hid")]
names(key)[2] <- "hid"
if (any(is.na(key$PersonId) | is.na(key$hid))) stop("missing id(s)")
if (any(duplicated(key$PersonId)) | any(duplicated(key$hid))) {
  stop("duplicated id(s)")
}
if (nrow(key) != nrow(smpl4)) stop("missing participant(s)")

# Save key
saveRDS(key, "data/key_fso_redcap_v4.rds")
