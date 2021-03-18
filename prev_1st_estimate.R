library(readxl)
library(writexl)

# Working directory
setwd("~/Projects/SerocoViD")

# Data
data_file <- "data-raw/dataFSOsero_04.06.2020_v1.1.txt"
data <- read.table(header = TRUE, file = data_file)
data$Date.of.birth.FSO <- as.Date(data$Date.of.birth.FSO, format = "%d.%m.%Y")
smpl <- read.csv("data-fso/COVID19_VD_V1_Total.csv", sep = ";")
smpl$dateOfBirth <- as.Date(smpl$dateOfBirth)
pop <- read_xlsx("data-fso/POPULATION PAR STRATES.xlsx", col_names = FALSE,
                 range = "B7:C14")
pop <- as.data.frame(pop)
colnames(pop) <- c("stratum", "N")

# Adds stratum variable to the data
tmp <- unique(smpl[!is.na(smpl$strate), c("dateOfBirth", "strate")])
names(tmp)[2] <- "stratum"
data <- merge(data, tmp, by.x = "Date.of.birth.FSO", by.y = "dateOfBirth",
              all.x = TRUE, sort = FALSE)
if (any(duplicated(data$uc_info_participant_hid))) stop()
if (any(is.na(data$stratum))) stop()
rm(tmp)

# Strata size
strata <- aggregate(uc_info_participant_hid ~ stratum, data, length)
names(strata)[2] <- "n"
strata <- merge(strata, pop, by = "stratum")

# Prevalence by stratum
npos <- aggregate(serol ~ stratum, data, sum)
names(npos)[2] <- "npos"
prev <- aggregate(serol ~ stratum, data, mean)
names(prev)[2] <- "p"
strata <- merge(strata, npos, by = "stratum")
strata <- merge(strata, prev, by = "stratum")
strata$var <- with(strata, (N - n) / ((n - 1) * N) * p * (1 - p))
strata$half.ci <- with(strata, qt(0.975, n - 1) * sqrt(var) + 1 / (2 * n))
strata$lwr <- pmax(0, strata$p - strata$half.ci)
strata$upr <- strata$p + strata$half.ci
rm(npos, prev)

# Unweighted prevalence and variance
u_p <- mean(data$serol)
u_N <- sum(pop$N)
u_n <- nrow(data)
u_var <- (u_N - u_n) / ((u_n - 1) * u_N) * u_p * (1 - u_p)

# Weighted prevalence and variance
w_p  <- with(strata, sum(p * N / sum(N)))
w_var <- with(strata, sum(N^2 * var) / sum(N)^2)

# Total
total <- data.frame(
  type = c("unweighted", "weighted"),
  p = c(u_p, w_p),
  var = c(u_var, w_var)
)
total$lwr <- total$p - qnorm(0.975) * sqrt(total$var)
total$upr <- total$p + qnorm(0.975) * sqrt(total$var)
rm(u_p, u_N, u_n, u_var, w_p, w_var)

# Export results
L <- list(strata = strata[c("stratum", "N", "n", "npos", "p", "lwr", "upr")],
          total = total[c("type", "p", "lwr", "upr")])
write_xlsx(L, path = "results/prev_1st_estimate_20200604.xlsx")

