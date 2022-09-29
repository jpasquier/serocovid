library(readxl)

# Working directory
setwd("~/Projects/SerocoViD")

# Date of birth
dob3 <- read_xlsx("data-fso/COVID19_VD_V2_REDCap_et_Publipostage.xlsx",
                 sheet = "Feuil1", guess_max = 10^4)
dob3 <- as.data.frame(dob3[c("HID", "dateOfBirth")])
names(dob3)[names(dob3) == "HID"] <- "hid"
dob3$dateOfBirth <- as.Date(dob3$dateOfBirth)
if (any(duplicated(dob3$hid))) stop("duplicated hid")
if (any(is.na(dob3$dateOfBirth))) stop("missing date of birth")

# Stratum
smpl <- read.csv("data-fso/COVID19_VD_V2_Total.csv", sep = ";")
smpl$dateOfBirth <- as.Date(smpl$dateOfBirth)
tmp <- unique(smpl[!is.na(smpl$strate), c("dateOfBirth", "strate")])
tmp <- do.call(rbind, lapply(sort(unique(tmp$strate)), function(i) {
  z <- tmp[tmp$strate == i, "dateOfBirth"]
  data.frame(dateOfBirth = seq(min(z), max(z), by = 1), stratum = i)
}))
if (!all(aggregate(stratum ~ dateOfBirth, tmp, function(x)
           length(unique(x)))[, 2] == 1)) {
  stop("ill-defined strata")
}

# Merge: get a stratum for each hid
dob3 <- merge(dob3, tmp, by = "dateOfBirth", all.x = TRUE)
if (any(is.na(dob3$stratum))) stop("missing stratum")
if (any(duplicated(dob3$hid))) stop("duplicated hid")
rm(smpl, tmp)

# The strata are renumbered so that they correspond to those of surveys 1 and 4
dob3$stratum <- dob3$stratum + 3

# Save strata
saveRDS(dob3[c("hid", "stratum")], "data/strata_v3.rds", compress = "xz")
write.table(dob3[c("hid", "stratum")], "data/strata_v3.csv", sep = ";",
            row.names = FALSE)
