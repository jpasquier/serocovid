library(parallel)
library(survey)
library(RCurl)
library(writexl)
library(XML)

# Working directory
setwd("~/Projects/SerocoViD")

# ------------------------------- REDCap data ------------------------------- #

# API
uri   <- "https://redcap.unisante.ch/api/"

# Tokens
tokens <- c("corona_immunitas", "personal_data")
tokens <- lapply(setNames(tokens, tokens), function(z) {
  z <- paste0("misc/redcap_", z, ".token")
  readChar(z, file.info(z)$size - 1)
})

# Import data with the API (RCurl)
tmp_file <- "/mnt/ramdisk/serocovid_data.rda"
if (file.exists(tmp_file)) {
  load(tmp_file)
} else {
  serocovid_data <- mclapply(tokens, function(token) {
    read.csv(text = postForm(
      uri = uri,
      token = token,
      content = 'record',
      format = 'csv'
    ))
  })
  save(serocovid_data, file = tmp_file)
}
rm(uri, tokens, tmp_file)

# -------------------------------- Serology --------------------------------- #

# Select variables
V_serol <- grep("covigg.+v4$", names(serocovid_data$corona_immunitas),
                value = TRUE)
V_serol <- c("hid", V_serol, "bl_vac_yn_ph4", "uc_labo_qc_v4")
pdta <- unique(serocovid_data$personal_data[
  c("uc_s_participant_hid", "uc_s_type_participant")])
serol <- serocovid_data$corona_immunitas[V_serol]
serol <- merge(serol, pdta, by.x = "hid", by.y = "uc_s_participant_hid",
               all.x = TRUE)
if (any(is.na(serol$uc_s_type_participant))) {
  stop("missing type of participant")
}
rm(pdta, V_serol)

# Select observation of survey 4
serol <- serol[serol$uc_s_type_participant %in% 14, ]
serol$uc_s_type_participant <- NULL

# Remove fake participants
hid_fake_participants <- c("13131313", "987654321")
serol[serol$hid %in% hid_fake_participants, ]
serol <- serol[!(serol$hid %in% hid_fake_participants), ]
rm(hid_fake_participants)

# Serology
if (FALSE) {
  with(serol, table(uc_labo_coviggl_v4, bl_vac_yn_ph4, useNA = "ifany"))
  with(serol, table(uc_labo_coviggl_v4, uc_labo_qc_v4, useNA = "ifany"))
}
u0 <- "uc_labo_qc_v4"
u1 <- "uc_labo_coviggl_v4"
if (any(is.na(serol[[u0]]) | is.na(serol[[u1]]))) stop("missing values")
serol$serol <- 9
serol[serol[u0] != "OK" | serol[[u1]] == "", "serol"] <- NA
serol[serol[u0] == "OK" & serol[[u1]] == "négatif", "serol"] <- 0
serol[serol[u0] == "OK" & serol[[u1]] == "positif", "serol"] <- 1
if (any(!is.na(serol$serol) & serol$serol == 9)) {
  stop("error in serology recoding")
}
rm(u0, u1)

###############################################################################
# TEMPORARY !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
#serol$serol <- rbinom(nrow(serol), 1, 0.8)
###############################################################################

# ----------------------- FSO sample (original file) ------------------------ #

# Import sample and select variables
V_fso <- c("PersonId", "HouseholdSizeSRPH", "sex", "maritalStatus",
           "nationalityState", "strate")
smpl <- read.csv("data-fso/COVID19_VD_V4_Total.csv", sep = ";",
                  fileEncoding = "latin1")
smpl <- smpl[!is.na(smpl$SAMPLE), V_fso]
rm(V_fso)

# Key between FSO and REDCap data
key <- readRDS("data/key_fso_redcap_v4.rds")

# Add hid variable and remove PersonId
smpl <- merge(smpl, key, by = "PersonId", all.x = TRUE)
if (any(is.na(smpl$hid))) stop("missing hid")
if (any(duplicated(smpl$PersonId)) | any(duplicated(smpl$hid))) {
  stop("Duplicated id")
}
smpl$PersonId <- NULL

# Add serology data
smpl <- merge(smpl, serol, by = "hid", all.x = TRUE)

# Recode FSO variables
names(smpl)[names(smpl) == "maritalStatus"] <- "maritalstatus"
names(smpl)[names(smpl) == "nationalityState"] <- "nationality"
names(smpl)[names(smpl) == "HouseholdSizeSRPH"] <- "householdsize"
smpl$strate <- factor(smpl$strate, levels = 1:8)
smpl$sex <- factor(smpl$sex, levels = 1:2, labels = c("Male", "Female"))
smpl$maritalstatus <- factor(smpl$maritalstatus, levels = 1:7)
levels(smpl$maritalstatus) <- c("Single", "Married", "Widowed", "Divorced",
                                 rep("Others", 3))
smpl$nationality <- factor(ifelse(smpl$nationality == 8100, 1, 2),
                           levels = 1:2, labels = c("Swiss", "Foreign"))
smpl$householdsize <- ifelse(smpl$householdsize >= 5, 5, smpl$householdsize)
smpl$householdsize <- factor(smpl$householdsize, levels = 1:5,
                             labels = c("1person", paste0(2:4, "people"),
                                        "5ormore"))

# Vaccination and Strate x Vaccination variables
smpl <- within(smpl, {
  vac <- c("y", "n", NA)[bl_vac_yn_ph4]
  stratevac <- ifelse(is.na(vac), NA, paste0(strate, vac))
})
smpl$stratevac <- factor(smpl$stratevac,
                         as.vector(outer(1:8, c("n", "y"), paste0)))
if(FALSE) with(smpl, table(stratevac, bl_vac_yn_ph4, strate, useNA = "ifany"))

# Respondents
smpl$resp <- as.numeric(!is.na(smpl$serol))
smpl$respvac <- as.numeric(!is.na(smpl$serol) & !is.na(smpl$vac))
if(FALSE) with(smpl, table(resp, respvac, useNA = "ifany"))

# ---------------------------- Sampling weights ----------------------------- #

# Population totals (FSO)
N <- read.table(header = TRUE, text = "
   strate       N
        1   38206
        2   43179
        3   43666
        4   42434
        5  216839
        6  273416
        7   67019
        8   64770
")

# Gross sample totals
n <- aggregate(hid ~ strate, smpl, length)
names(n)[2] <- "n"
N <- merge(N, n, by = "strate")
rm(n)

# Sampling weigths
N$wsmpl <- N$N / N$n
smpl <- merge(smpl, N, by = "strate")

# --------------------------- Non-response model ---------------------------- #

# Model
nrm <- glm(resp ~ strate + sex + maritalstatus + nationality + householdsize,
           family = binomial, data = smpl)
summary(nrm)

# ------------------------------- Calibration ------------------------------- #

# Vaccination rates (e-mail VF 8 oct 2021)
vac_rate <- read.table(header = TRUE, text = "
   strate    rate
        1       0
        2  0.0001
        3  0.1943
        4  0.6739
        5  0.7451
        6  0.8244
        7  0.8944
        8  0.9145
")

# Vaccination margins - table
z <- within(merge(N, vac_rate, by = "strate"), {
  yes <- round(N * rate)
  no <- round(N * (1 - rate))
})
vac_margins <- data.frame(
  variable = "vac",
  value = c("n", "y"),
  N = apply(z[c("no", "yes")], 2, sum)
)
stratevac_margins <- cbind(variable = "stratevac", rbind(
  data.frame(value = paste0(z$strate, "n"), N = z$no),
  data.frame(value = paste0(z$strate, "y"), N = z$yes)
))
rm(z)

# Margins - table (FSO)
margins <- read.table(header = TRUE, text = "
        variable       value        N
             sex        Male   387686
             sex      Female   401843
     nationality       Swiss   537298
     nationality     Foreign   252231
   maritalstatus      Single   375217
   maritalstatus     Married   309291
   maritalstatus     Widowed    31012
   maritalstatus    Divorced    71583
   maritalstatus      Others     2426
   householdsize     1person   134990
   householdsize     2people   214269
   householdsize     3people   152734
   householdsize     4people   186887
   householdsize     5ormore   100649
")
s <- cbind(variable = "strate", setNames(N[c("strate", "N")], c("value", "N")))
margins <- rbind(s, margins, vac_margins, stratevac_margins)
rm(s)

# Margins - vector
tot <- unique(aggregate(N ~ variable, margins, sum)$N)
margins <- do.call(base::c, lapply(unique(margins$variable), function(v) {
  tmp <- margins[margins$variable == v, ][-1, ]
  setNames(tmp$N, paste0(v, tmp$value))
}))
margins <- c(`(Intercept)` = tot, margins)
rm(tot)


#################################



# net sample - all participants with a serology
net <- aggregate(resp ~ strate, smpl, mean)
names(net)[2] <- "resprate"
net <- merge(subset(smpl, resp == 1), net, by = "strate")
net$w0 <- net$wsmpl / net$resprate
net <- svydesign(id = ~hid, strata = ~strate, weights = ~w0, fpc = ~N,
                 data = net)

# net sample - participants with a serology and a vaccination status
vnet <- aggregate(respvac ~ strate, smpl, mean)
names(vnet)[2] <- "resprate"
vnet <- merge(subset(smpl, respvac == 1), vnet, by = "strate")
vnet$w0 <- vnet$wsmpl / vnet$resprate
vnet <- svydesign(id = ~hid, strata = ~strate, weights = ~w0, fpc = ~N,
                  data = vnet)

# prévalence avec calage (raking)
prev <- lapply(1:2, function(k) {
  if (k == 1) {
    d <- net$variables
    L <- 1:3
  } else {
    d <- vnet$variables
    L <- 1:9
  }
  n <- rbind(data.frame(strate = 0, n = nrow(d)),
             setNames(aggregate(serol ~ strate, d, length), c("strate", "n")))
  m <- rbind(data.frame(strate = 0, npos = sum(d$serol)),
             setNames(aggregate(serol ~ strate, d, sum),
                      c("strate", "npos")))
  p <- rbind(data.frame(strate = 0, prop = mean(d$serol)),
             setNames(aggregate(serol ~ strate, d, mean),
                      c("strate", "prop")))
  s <- merge(merge(n, m, by = "strate"), p, by = "strate")
  if (k == 2) {
    v0 <- rbind(data.frame(strate = 0, nvac = sum(d$vac == "y")),
                setNames(aggregate(I(vac == "y") ~ strate, d, sum),
                         c("strate", "nvac")))
    v1 <- rbind(data.frame(strate = 0,
                           nvacpos = sum(d$vac == "y" & d$serol == 1)),
                setNames(aggregate(I(vac == "y" & serol == 1) ~ strate, d, sum),
                         c("strate", "nvacpos")))
    s <- merge(merge(s, v0, by = "strate"), v1, by = "strate")
  }
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
    d <- calibrate(design = list(net, vnet)[[k]], formula = fml,
                   population = mrg, calfun = "raking")
    z0 <- svymean(~serol, d)
    z0  <- c(0, z0[[1]], confint(z0))
    z1 <- svyby(~serol, ~strate, d, svymean)
    z1 <- cbind(z1[, 1:2], confint(z1))
    z1$strate <- as.numeric(as.character(z1$strate))
    v <- c("strate", paste0("prev_cal", l, c("", "_lwr", "_upr")))
    s <- merge(s, rbind(setNames(z0, v), setNames(z1, v)), by = "strate")
  }
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
write_xlsx(prev, paste0("results/prev_4th_visit_", format(Sys.Date(), "%Y%m%d"), ".xlsx"))


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

