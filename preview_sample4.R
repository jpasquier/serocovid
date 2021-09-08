library(RCurl)
library(XML)

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

# ----------------------- FSO sample (original file) ------------------------ #

smpl <- read.csv("./data-fso/COVID19_VD_V4_Total.csv", sep = ";",
                 fileEncoding = "latin1")
smpl <- smpl[!is.na(smpl$SAMPLE), ]
id_smpl <- c("officialName", "firstName", "dateOfBirth")
any(duplicated(smpl[id_smpl]))

# ------------------------- REDCap - Personal data -------------------------- #

pdta4 <- subset(serocovid_data$personal_data, uc_s_type_participant %in% 14)
any(duplicated(pdta4$uc_s_participant_sid))
any(duplicated(pdta4$uc_s_participant_hid))
id_pdta4 <- c("uc_s_nom", "uc_s_prenom", "full_dob")
any(duplicated(pdta4[id_pdta4]))

# --------------------- key between FSO and REDCap data --------------------- #

key <- merge(
  x = smpl[c(id_smpl, "PersonId")],
  y = pdta4[c(id_pdta4, "uc_s_participant_sid", "uc_s_participant_hid")],
  by.x = id_smpl,
  by.y = id_pdta4,
  all = TRUE,
  sort = FALSE
)
nrow(key) == nrow(smpl)
