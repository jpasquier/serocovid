library(parallel)
library(readxl)
library(writexl)

options(mc.cores = detectCores())

# Working directory
setwd("~/Projects/SerocoViD")

# Intermediate data of the second wave
data_file <- "data-raw/Intermed_results_Q2017_20201123_date_prélèvmt.xlsx"
data <- read_xlsx(data_file, sheet = "màj_23.11", range = "A1:I1217")
data <- as.data.frame(data)
data <- data[!is.na(data$uc_labo_coviggl_v2), ]

# Date of birth
i <- grepl("^10[0-9]{2}-[0-9]{2}-[09]{2}$", data$DDN)
delta <- as.numeric(as.Date("1970-01-01") - as.Date("1899-12-30"))
data$DDN[i] <- as.character(as.numeric(as.Date(sub("^10", "19", data$DDN[i])))
                              + delta)
data$DDN <- as.Date(as.numeric(data$DDN), origin = "1899-12-30")
rm(i, delta)

# Strata
smpl <- read.csv("data-fso/COVID19_VD_V1_Total.csv", sep = ";")
smpl$dateOfBirth <- as.Date(smpl$dateOfBirth)
tmp <- unique(smpl[!is.na(smpl$strate), c("dateOfBirth", "strate")])
tmp <- do.call(rbind, lapply(sort(unique(tmp$strate)), function(i) {
  z <- tmp[tmp$strate == i, "dateOfBirth"]
  data.frame(DDN = seq(min(z), max(z), by = 1), stratum = i)
}))
data <- merge(data, tmp, by = "DDN", all.x = TRUE, sort = FALSE)
data[data$hid == 745932, "stratum"] <- 6
if (any(duplicated(data$hid))) stop("duplicated hid")
if (any(is.na(data$stratum))) stop("missing stratum")
rm(tmp)

# Sample structure
smp_str <- aggregate(hid ~ stratum, data, length)
names(smp_str)[2] <- "n"

# remove stratum 1 which includes only 1 participant
smp_str <- smp_str[smp_str$stratum != 1, ]

# proportion of the sample in each stratum
smp_str$prop_3_8 <- smp_str$n / sum(smp_str$n)
smp_str$prop_4_8 <- c(NA, smp_str$n[-1] /  sum(smp_str$n[-1]))

# Population structure in the FSO register
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
pop$prop_3_8 <- c(NA, NA, pop$N[3:8] /  sum(pop$N[3:8]))
pop$prop_4_8 <- c(rep(NA, 3), pop$N[4:8] /  sum(pop$N[4:8]))


# Function for the estimation of the prevalence
# Variance: Cochran, W. G., Sampling techniques John Wiley & Sons, 1977,
#           page 144, formula 5A.75
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
  .mean$lwr <- .mean$y - qnorm(.975) * sqrt(.mean$v)
  .mean$upr <- .mean$y + qnorm(.975) * sqrt(.mean$v)
  return(.mean)
}

# Naive estimate
N_naive <- do.call(rbind, lapply(c(.02, .04, .06), function(prec) {
  do.call(rbind, lapply(seq(.15, .50, .05), function(prev) {
    n <- ceiling(qnorm(.975)^2 * prev * (1 - prev) / prec^2)
    data.frame(prev = prev, prec = prec, n_naive = n)
  }))
}))

# Simulations that use the same sample structure as the intermediate data of
# the second wave or the structure of the population
N_sim <- lapply(c("smp", "pop"), function(struct) {
  lapply(3:4, function(j) {
    do.call(rbind, lapply(seq(.15, .50, .05), function(prev) {
      p <- do.call(rbind, mclapply(100:3500, function(n) {
        if (struct == "smp") {
          smp <- smp_str
        } else {
          smp <- pop
        }
        smp$n <- round(n * smp[[paste0("prop_", j, "_8")]])
        sim <- do.call(rbind, lapply(j:8, function(k) {
          n <- smp[smp$stratum == k, "n"] 
          npos <- round(n * prev)
          nneg <- n - npos
          serol <- c(rep(1, npos), rep(0, nneg))
          data.frame(stratum = k, serol = serol)
        }))
        p <- Mean(sim)[c("y", "lwr", "n")]
        p$lwr <- p$y - p$lwr
        names(p)[1:2] <- c("prev", "prec")
        p$trg_prev <- prev
        p
      }))
      gc()
      p <- do.call(rbind, lapply(c(.02, .04, .06), function(trg_prec) {
        diff <- abs(p$prec - trg_prec)
        p <- p[!is.na(diff) & diff == min(diff, na.rm = TRUE), ]
        if (nrow(p) > 1) {
          diff <- abs(p$prev - p$trg_prev)
          p <- p[which(diff == min(diff)), ]
        }
        if (nrow(p) > 1) {
          p <- p[1, ]
        }
        p$trg_prec <- trg_prec
        p
      }))
      names(p)[1:3] <- paste(names(p)[1:3], struct, j, "8", sep = "_")
      p
    }))
  })
})
N_sim <- Reduce(
  function(x, y) merge(x, y, by = c("trg_prev", "trg_prec"), all = TRUE),
  unlist(N_sim, recursive = FALSE)
)

# Export results
N <- merge(N_naive, N_sim, by.x = c("prev", "prec"),
           by.y = c("trg_prev", "trg_prec"), all = TRUE)
N <- N[order(N$prec, N$prev), ]
write.table(N, "results/sample_size_20201203.csv", sep = ";",
            row.names = FALSE)
write_xlsx(N, "results/sample_size_20201203.xlsx")
