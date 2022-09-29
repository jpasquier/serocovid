setwd("~/Projects/SerocoViD")
files <- list(
  orig = "data-fso/COVID19_VD_V1_Total.csv",
  serv = "data-fso/COVID19_VD_V1_Total_SERV.csv"
)
smpl <- lapply(files, function(f) {
  dta <- read.csv(f, sep = ";")
  if (any(is.na(dta$strate) != is.na(dta$SAMPLE))) stop(paste("error in", f))
  dta[!is.na(dta$SAMPLE), ]
})
sapply(smpl, nrow)
sapply(smpl, function(z) any(duplicated(z)))
V <- c("PersonId", "dateOfBirth", "strate", "SAMPLE")
cmp <- merge(smpl$orig[V], smpl$serv[V], by = "PersonId", all = TRUE,
             suffixes = c(".orig", ".serv"))
`%~~%` <- function(x, y) !is.na(x) & !is.na(y) & x == y | is.na(x) & is.na(y)
i <- with(cmp, !(dateOfBirth.orig %~~% dateOfBirth.serv) |
               !(strate.orig %~~% strate.serv) |
               !(SAMPLE.orig %~~% SAMPLE.serv))
cmp[i, ]
