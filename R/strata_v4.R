
# Working directory
setwd("~/Projects/SerocoViD")

# Import sample and select variables
str4 <- read.csv("data-fso/COVID19_VD_V4_Total.csv", sep = ";",
                  fileEncoding = "latin1")
str4 <- str4[!is.na(str4$SAMPLE), c("PersonId", "strate")]
names(str4)[names(str4) == "strate"] <- "stratum"

# Add hid variable and remove PersonId
str4 <- merge(str4, readRDS("data/key_fso_redcap_v4.rds"),
               by = "PersonId", all.x = TRUE)
if (any(is.na(str4$hid))) stop("missing hid")
if (any(duplicated(str4$PersonId)) | any(duplicated(str4$hid))) {
  stop("Duplicated id")
}
str4$PersonId <- NULL

# Save strata
saveRDS(str4[c("hid", "stratum")], "data/strata_v4.rds", compress = "xz")
write.table(str4[c("hid", "stratum")], "data/strata_v4.csv", sep = ";",
            row.names = FALSE)
