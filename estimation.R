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
                                 (p * (1 - p) + (1 - m / n) * (p - P)))) / M^2
    data.frame(domain = z, N = M, n = sum(m[2]), p = P, v = V)
  }))
  .prev$lwr <- .prev$p - qnorm(0.975) * sqrt(.prev$v)
  .prev$upr <- .prev$p + qnorm(0.975) * sqrt(.prev$v)
  return(.prev)
}

# Examples
data$stratum1 <- as.numeric(data$stratum == 1)
data$stratum2 <- as.numeric(data$stratum == 2)
data$stratum3 <- as.numeric(data$stratum == 3)
data$stratum4 <- as.numeric(data$stratum == 4)
data$stratum5 <- as.numeric(data$stratum == 5)
data$stratum6 <- as.numeric(data$stratum == 6)
data$stratum7 <- as.numeric(data$stratum == 7)
data$stratum8 <- as.numeric(data$stratum == 8)
data$genre1 <- as.numeric(data$uc_info_genre == 1)
data$genre2 <- as.numeric(data$uc_info_genre == 2)
data$strata6genre1 <- data$stratum6 * data$genre1
dlist <- c(paste0("stratum", 1:8), c("genre1", "genre2", "strata6genre1"))
prev(domain = dlist)


