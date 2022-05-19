setwd("~/Projects/SerocoViD")
smpl3 <- read.csv("data-fso/COVID19_VD_V2_Total.csv", sep = ";",
                  fileEncoding = "latin1")
smpl3 <- smpl3[!is.na(smpl3$SAMPLE), ]
smpl4 <- read.csv("data-fso/COVID19_VD_V4_Total.csv", sep = ";",
                  fileEncoding = "latin1")
smpl4 <- smpl4[!is.na(smpl4$SAMPLE), ]
key3 <- readRDS("data/key_fso_redcap_v3.rds")
key4 <- readRDS("data/key_fso_redcap_v4.rds")
gender <- rbind(
  cbind(survey = 3, merge(key3, smpl3, by = "PersonId")[c("hid", "sex")]),
  cbind(survey = 4, merge(key4, smpl4, by = "PersonId")[c("hid", "sex")])
)
gender$sex2 <- factor(gender$sex, 1:2, c("Male", "Female"))
write.table(gender, "~/gender_serocovid.csv", sep = ";", row.names = FALSE)
