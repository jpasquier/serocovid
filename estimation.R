library(readxl)
library(writexl)

# Working directory
setwd("~/Projects/SerocoViD")

# Data
load("data/dataFSOsero_08072020.RData")

# Population size
pop <- read.table(header = TRUE, text = "
  stratum       N
        1   37490
        2   42538
        3   42700
        4   42171
        5  212336
        6  268898
        7   67412
        8   63097
")

# Estimation
# Variance: Cochran, W. G., Sampling techniques John Wiley & Sons, 1977,
#           page 144, formula 5A.75
# domain: any variable

data <- dataFSOsero
variable <- "serol"
stratum <- "stratum"
domain <- "uc_info_genre"
d <- domain[1]

Mean <- function(data = data, variable = "serol", stratum = "stratum",
                 domain = NULL, .pop = pop) {
  # select observed values
  data <- data[!is.na(data[[variable]]) & !is.na(data[[stratum]]), ]
  # formula variable ~ stratum
  f <- as.formula(paste(variable, "~", stratum))
  # stratum sizes
  strata <- aggregate(f, data, length)
  names(strata) <- c("stratum", "n")
  strata <- merge(pop, strata, by = "stratum")
  # convert N and n to numeric to avoid the 'integer overflow' problem
  strata$N <- as.numeric(strata$N)
  strata$n <- as.numeric(strata$n)
  # no domain
  if (is.null(domain)) {
    data$.all <- 1
    domain <- ".all"
  }
  # convert domains to factors
  for (d in domain) {
    data[[d]] <- factor(data[[d]])
  }
  .mean <- do.call(rbind, lapply(domain, function(d) {
    do.call(rbind, lapply(levels(data[[d]]), function(z) {
      data.sub <- data[data[[d]] == z, c(stratum, variable)]
      # sample size of domain by strata
      m <- aggregate(f, data.sub, length)
      names(m) <- c("stratum", "m")
      # mean by domain and strata
      y <- aggregate(f, data.sub, mean)
      names(y) <- c("stratum", "y")
      # sum of squares by domain and strata
      ss <- aggregate(f, data.sub, function(u) sum((u - mean(u))^2))
      names(ss) <- c("stratum", "ss")
      # merge prevalences and sizes with strata
      strata.sub <- merge(strata, m, by = "stratum")
      strata.sub <- merge(strata.sub, y, by = "stratum")
      strata.sub <- merge(strata.sub, ss, by = "stratum")
      #
      M <- sum(with(strata.sub, N / n * m))
      Y <- sum(with(strata.sub, N / n * m * y)) / M
      V <- sum(with(strata.sub, N * (N - n) / (n * (n - 1)) * 
                                  (ss + m * (1 - m / n) * (y - Y)^2))) / M^2
      data.frame(domain = d, value = z, N = M, n = sum(m[2]), y = Y, v = V)
    }))
  }))
  .mean$lwr <- .mean$y - qnorm(0.975) * sqrt(.mean$v)
  .mean$upr <- .mean$y + qnorm(0.975) * sqrt(.mean$v)
  return(.mean)
}


# Examples
Mean(dataFSOsero)
Mean(dataFSOsero, domain = c("stratum", "uc_info_genre"))
dataFSOsero$strate_genre <- with(dataFSOsero, stratum * 10 + uc_info_genre)
Mean(dataFSOsero, domain = "strate_genre")
Mean(dataFSOsero, variable = "age")
Mean(dataFSOsero, variable = "age", domain = c("stratum", "uc_info_genre"))

# --------------------------------------------------------------------------- #
# --------------------------------------------------------------------------- #
# --------------------------------------------------------------------------- #

# With the survey package
library(survey)
tmp <- dataFSOsero[!is.na(dataFSOsero$serol), ]
wgt <- aggregate(serol ~ stratum, tmp, length)
names(wgt)[2] <- "n"
wgt <- merge(wgt, pop, by = "stratum")
wgt$wgt <- wgt$N / wgt$n
tmp <- merge(tmp, wgt, by = "stratum")
d <- svydesign(
  id = ~uc_info_participant_hid,
  strata = ~stratum,
  weights = ~wgt,
  data = tmp,
  fpc = ~N
)
(prev <- svymean(~serol, d))
confint(prev)
confint(prev, df = sum(wgt$n) - nrow(wgt)) # loi de Student au lieu de la loi normale
(prev_by_gender <- svyby(~serol, ~uc_info_genre, design = d, svymean))
confint(prev_by_gender)
(age_by_gender <- svyby(~age, ~uc_info_genre, design = d, svymean))
confint(age_by_gender)

# chisq test with the survey package
svychisq(~serol+uc_info_genre, d, statistic = "Chisq")

fit <- svyglm(serol ~ uc_info_genre, d, family = quasibinomial)
cbind(exp(cbind(odds_ratio = coef(fit), confint(fit))), p_value = coef(summary(fit))[, 4])[-1, , drop = FALSE]
