library(readxl)
library(writexl)

# Working directory
setwd("~/Projects/SerocoViD")

# Data
load("data/dataFSOsero_08072020.RData")
data_file <- "data-raw/extractforJP_29.6.2020.csv"
data <- read.csv(data_file)
data$X <- NULL
data <- data[!is.na(data$serol), ]
data$Date.of.birth.FSO <- as.Date(data$Date.of.birth.FSO)
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
#if (any(duplicated(data$uc_info_participant_hid))) stop()
if (any(is.na(data$stratum))) stop()
rm(tmp)

# Strata size
strata <- aggregate(uc_info_participant_hid ~ stratum, data, length)
names(strata)[2] <- "n"
strata <- merge(strata, pop, by = "stratum")

# Estimation
# Variance: Cochran, W. G., Sampling techniques John Wiley & Sons, 1977,
#           page 144, formula 5A.75
# domain: dummy variable, =1 if the observation is in the domain, =0 if not

prev <- function(.data = data, .strata = strata, domain = NULL) {
  if (is.null(domain)) {
    .data$.all <- 1
    domain <- ".all"
  }
  .prev <- do.call(rbind, lapply(domain, function(z) {
    .data.sub <- .data[.data[[z]] == 1, c("stratum", "serol")]
    # prevalence by domain and strata
    p <- aggregate(serol ~ stratum, .data.sub, mean)
    names(p)[2] <- "p"
    # sample size of domain by strata
    m <- aggregate(serol ~ stratum, .data.sub, length)
    names(m)[2] <- "m"
    # merge prevalences and sizes with strata
    .strata.sub <- merge(.strata, p, by = "stratum")
    .strata.sub <- merge(.strata.sub, m, by = "stratum")
    #
    M <- sum(with(.strata.sub, N / n * m))
    P <- sum(with(.strata.sub, N / n * m * p)) / M
    V <- sum(with(.strata.sub, N * (N - n)/ (n * (n - 1)) * m * 
                                 (p * (1 - p) + (1 - m / n) * (p - P)^2))) / M^2
    data.frame(domain = z, N = M, n = sum(m[2]), p = P, v = V)
  }))
  .prev$lwr <- .prev$p - qnorm(0.975) * sqrt(.prev$v)
  .prev$upr <- .prev$p + qnorm(0.975) * sqrt(.prev$v)
  return(.prev)
}
b <- prev()


with(dataFSOsero, table(district, stratum))

do.call(rbind, lapply(unique(dataFSOsero$district), function(d) {
  dataFSOsero[[d]] <- as.numeric(dataFSOsero$district == d)
  prev(.data = dataFSOsero, domain = d)
}))
with(dataFSOsero[dataFSOsero$district == "Lausanne District", ],
     table(stratum, serol))
.data.sub <- dataFSOsero[dataFSOsero$district == "Lausanne District", ]
