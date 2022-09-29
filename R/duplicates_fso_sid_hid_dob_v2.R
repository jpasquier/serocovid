setwd("~/Projects/SerocoViD")
dob1 <- read.csv("data-raw/fso_sid_hid_dob_v2.csv", sep = ";")
i <- is.na(dob1$HID) | is.na(dob1$Date.of.birth)
if (any(i)) dob1[i, ]
i <- !grepl("^[0-3][0-9]\\.[0-1][0-9]\\.(19|20)[0-9]{2}$", dob1$Date.of.birth)
if (any(i)) dob1[i, ]
dob1$Date.of.birth <- as.Date(dob1$Date.of.birth, format = "%d.%m.%Y")
dup <- dob1[dob1$HID %in% dob1$HID[duplicated(dob1$HID)], ]
write.table(dup[order(dup$HID), ], sep = ";", row.names = FALSE,
            file = "~/duplicates_fso_sid_hid_dob_v2.csv")
rm(i, dup)
