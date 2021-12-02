library(parallel)
library(survey)
library(writexl)

options(mc.cores = detectCores() - 1)

# Working directory
setwd("~/Projects/SerocoViD")

# Prepare data
source("R/prepare_data_4th_visit.R")

# Non-response model
nrm <- glm(resp ~ strate + sex + maritalstatus + nationality + householdsize,
           family = binomial, data = smpl)
summary(nrm)

# Prevalence with calibration (raking)
prev <- lapply(1:2, function(k) {
  smpl$all <- "All"
  smpl$grp1 <- c(rep("1: 6m-19", 4), rep("2: 20-64", 2),
                 rep("3: >=65", 2))[smpl$strate]
  smpl$grp2 <- c("6m-4", rep(">=5", 7))[smpl$strate]
  smpl$grp3 <- c(rep("1: 6m-9", 2), rep("2: >=10", 6))[smpl$strate]
  smpl$grp4 <- c(rep("1: 6m-14", 3), rep("2: >=15", 5))[smpl$strate]
  smpl$grp5 <- c(rep("1: 6m-19", 3), rep("2: >=20", 5))[smpl$strate]
  smpl$vac2 <- with(smpl, ifelse(is.na(vac), "3: missing",
                                ifelse(vac == "y", "2: yes", "1: no")))
  smpl$vac2 <- paste("vac:", smpl$vac2)
  smpl$vac0 <- c(n = 0, y = 1)[smpl$vac]
  smpl$vacpos <- smpl$serol * smpl$vac0
  smpl$serol2 <- c("seroneg", "seropos")[smpl$serol + 1]
  if (k == 1) {
    net <- aggregate(resp ~ strate, smpl, mean)
    names(net)[2] <- "resprate"
    net <- merge(subset(smpl, resp == 1), net, by = "strate")
    L <- 1:3
  } else {
    net <- aggregate(respvac ~ strate, smpl, mean)
    names(net)[2] <- "resprate"
    net <- merge(subset(smpl, respvac == 1), net, by = "strate")
    L <- 1:9
  }
  net$w0 <- net$wsmpl / net$resprate
  net <- svydesign(id = ~hid, strata = ~strate, weights = ~w0, fpc = ~N,
                   data = net)
  d <- net$variables
  strate_names <- c("6m-4", "5-9", "10-14", "15-19", "20-39", "40-64", "65-74",
                    ">=75")
  fmls <- list(~all, ~strate, ~grp1, ~grp2, ~grp3, ~grp4, ~grp5, ~vac2,
               ~stratevac, ~serol2)
  s <- do.call(rbind, mclapply(1:length(fmls), function(j) {
    fml <- update(fmls[[j]], serol ~ .)
    n <- setNames(aggregate(fml, d, length), c("domaine", "n"))
    m <- setNames(aggregate(fml, d, sum), c("domaine", "npos"))
    p <- setNames(aggregate(fml, d, mean), c("domaine", "prop"))
    s <- merge(merge(n, m, by = "domaine"), p, by = "domaine")
    if (k == 2) {
      fml <- update(fml, I(vac == "y") ~ .)
      v0 <- setNames(aggregate(fml, d, sum), c("domaine", "nvac"))
      fml <- update(fml, I(vac == "y" & serol == 1) ~ .)
      v1 <- setNames(aggregate(fml, d, sum), c("domaine", "nvacpos"))
      s <- merge(merge(s, v0, by = "domaine"), v1, by = "domaine")
    }
    if (any(grepl("^strate$", fml))) {
      s$domaine <- strate_names[as.numeric(as.character(s$domaine))]
    }
    if (any(grepl("^grp2$", fml))) {
      s <- s[s$domaine != "6m-4", ]
    }
    if (any(grepl("^grp5$", fml))) {
      s <- s[s$domaine != "1: 6m-19", ]
    }
    if (any(grepl("^stratevac$", fml))) {
      s$domaine <- paste("stratevac:", as.character(s$domaine))
    }
    cbind(domaine_index = j * 100 + 1:nrow(s), s)
  }))
  for (l in L) {
    i.mrg <- grep("Intercept", names(margins))
    if (l <= 3) {
      fml <- ~ strate
      i.mrg <- c(i.mrg, grep("strate[1-8]", names(margins)))
    } else if (l <= 6) {
      fml <- ~ strate + vac
      i.mrg <- c(i.mrg, grep("strate[1-8]", names(margins)))
      i.mrg <- c(i.mrg, grep("^vac", names(margins)))
    } else {
      fml <- ~ stratevac
      i.mrg <- c(i.mrg, grep("stratevac", names(margins)))
    }
    if (l %% 3 == 2) {
      fml <- update(fml, ~ . + nationality)
      i.mrg <- c(i.mrg, grep("nationality", names(margins)))
    } else if (l %% 3 == 0) {
      fml <-
        update(fml, ~ . + sex + nationality + maritalstatus + householdsize)
      i.mrg <- c(i.mrg, grep("sex|natio|marit|house", names(margins)))
    }
    mrg <- margins[i.mrg]
    d <- calibrate(design = net, formula = fml, population = mrg,
                   calfun = "raking")
    z <- do.call(rbind, mclapply(fmls, function(fml) {
      Merge <- function(x, y) merge(x, y, by = "domaine")
      U <- c("", "_low", "_medium", "_high")
      z <- Reduce(Merge, lapply(U, function(u) {
        z <- svyby(as.formula(paste0("~serol", u)), fml, d, svymean)
        z <- cbind(z[, 1:2], confint(z))
        v <- c("domaine", paste0("prev", u, "_cal", l, c("", "_lwr", "_upr")))
        z <- setNames(z, v)
      }))
      if (k == 2 & l <= 3) {
        z1 <- svyby(~vac0, fml, d, svymean)
        z1 <- cbind(z1[, 1:2], confint(z1))
        v1 <- c("domaine", paste0("vac_rate_cal", l, c("", "_lwr", "_upr")))
        z1 <- setNames(z1, v1)
        z2 <- svyby(~vacpos, fml, d, svymean)
        z2 <- cbind(z2[, 1:2], confint(z2))
        v2 <- c("domaine", paste0("vacpos_rate_cal", l, c("", "_lwr", "_upr")))
        z2 <- setNames(z2, v2)
        z <- Reduce(Merge, list(z, z1, z2))
      }
      if (any(grepl("^strate$", fml))) {
        z$domaine <- strate_names[as.numeric(as.character(z$domaine))]
      }
      if (any(grepl("^grp2$", fml))) {
        z <- z[z$domaine != "6m-4", ]
      }
      if (any(grepl("^grp5$", fml))) {
        z <- z[z$domaine != "1: 6m-19", ]
      }
      if (any(grepl("^stratevac$", fml))) {
        z$domaine <- paste("stratevac:", as.character(z$domaine))
      }
      z
    }))
    s <- merge(s, z, by = "domaine")
  }
  s <- s[order(s$domaine_index), ]
  s$domaine_index <- NULL
  s$domaine <- sub("^[0-9]: ", "", s$domaine)
  s
})
prev[[3]] <- data.frame(calage = 1:9, variables = c(
 "strate (pas de calage)",
 "strate + nationalité (suisse/étranger)",
 "strate + sexe + nationalité + état civil + taille du ménage",
 "strate + vaccination",
 "strate + vaccination + nationalité (suisse/étranger)",
 "strate + vaccination + sexe + nationalité + état civil + taille du ménage",
 "strate * vaccination",
 "strate * vaccination + nationalité (suisse/étranger)",
 "strate * vaccination + sexe + nationalité + état civil + taille du ménage"
))
write_xlsx(prev, paste0("results/prev_4th_visit_",
                        format(Sys.Date(), "%Y%m%d"), ".xlsx"))


###############################################################################
pop2 <- N[c("strate", "N")]
names(pop2)[1] <- "stratum"
Mean <- function(data = smpl, variable = "serol", stratum = "strate",
                 domain = NULL, pop = pop2) {
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
Mean(variable = "serol")

