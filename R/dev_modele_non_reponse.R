library(REDCapR)
library(readxl)
library(survey)

# Working directory
setwd("~/Projects/SerocoViD")


# ------------------------------- REDCap data ------------------------------- #

# API
uri   <- "https://redcap.unisante.ch/api/"

# Tokens
tokens <- c("corona_immunitas", "personal_data", "research_data")
tokens <- lapply(setNames(tokens, tokens), function(z) {
  z <- paste0("misc/redcap_", z, ".token")
  readChar(z, file.info(z)$size - 1)
})

# Import
tmp_file <- "/tmp/serocovid_data.rda"
if (file.exists(tmp_file)) {
  load(tmp_file)
} else {
  serocovid_data <- lapply(tokens, function(token) {
    redcap_read(redcap_uri = uri, token = token)$data
  })
  save(serocovid_data, file = tmp_file)
}
rm(uri, tokens, tmp_file)

# ------------------------------- Gros sample ------------------------------- #

smpl <- read.csv("data-fso/COVID19_VD_V2_Total.csv", sep = ";",
                 encoding = "latin1")
if (any(is.na(smpl$strate) != is.na(smpl$SAMPLE))) stop()
smpl <- subset(smpl, !is.na(strate))

# Add hid to gross sample
hid <- read.csv("data-fso/link_fso_redcap_sample2.csv", sep = ";")
smpl <- merge(smpl, hid, by = "PersonId", all.x = TRUE)
if (any(is.na(smpl$hid))) stop("missing hid")
rm(hid)

# ---------------------------- Sampling weights ----------------------------- #

# Population totals
N <- read.table(header = TRUE, text = "
   strate       N
        1   42140
        2  214498
        3  270755
        4   67462
        5   63605
")

# Sample totals
n <- aggregate(hid ~ strate, smpl, length)
names(n)[2] <- "n"

# Sampling weigths
w <- merge(N, n, by = "strate")
w$w0 <- w$N / w$n
smpl <- merge(smpl, w[c("strate", "w0")], by = "strate")
rm(w)

# ------------------------------- Net sample -------------------------------- #

# Select IDs in personnal data
net <- subset(serocovid_data$personal_data, uc_s_type_participant %in% 13,
              "uc_s_participant_hid")
names(net) <- "hid"
if (any(duplicated(net$hid))) stop("duplicated hid")
vars <- c("hid", "uc_labo_coviggl_v2", "uc_labo_covigal_v2")
net <- merge(net, serocovid_data$corona_immunitas[vars], by = "hid")
net <- subset(net, !is.na(uc_labo_coviggl_v2) | !is.na(uc_labo_covigal_v2))
if (any(duplicated(net$hid))) stop("duplicated hid")
nrow(net)
net$respondent <- 1

# Add response status to gross sample
smpl <- merge(smpl, net, by = "hid", all = TRUE)
if (any(is.na(smpl$strate))) stop("missing strata")
smpl[is.na(smpl$respondent), "respondent"] <- 0

# ------------------------------ Naive weights ------------------------------ #

# Sample totals
m <- aggregate(respondent ~ strate, smpl, sum)
names(m)[2] <- "m"

# Sampling weigths
w <- merge(N, m, by = "strate")
w$w_naive <- w$N / w$m
smpl <- merge(smpl, w[c("strate", "w_naive")], by = "strate")
rm(w)

# --------------------------- Non-response model ---------------------------- #

# Recode categorial variables
smpl$agegroup <- factor(smpl$strate, levels = 1:5,
                        labels = c("15-19", "20-39", "40-64", "65-74", ">=75"))
smpl$sex <- factor(smpl$sex, levels = 1:2, labels = c("Male", "Female"))
smpl$maritalstatus <- factor(smpl$maritalStatus, levels = 1:7)
levels(smpl$maritalstatus) <- c("Single", "Married", "Widowed", "Divorced",
                                 rep("Others", 3))
smpl$nationality <- factor(ifelse(smpl$nationalityState == 8100, 1, 2),
                           levels = 1:2, labels = c("Swiss", "Foreign"))
smpl$householdsize <-
  ifelse(smpl$HouseholdSizeSRPH >= 5, 5, smpl$HouseholdSizeSRPH)
smpl$householdsize <- factor(smpl$householdsize, levels = 1:5,
                             labels = c("1_person", paste0(2:4, "_people"),
                                        "5_or_more"))

# Model
nrm <- glm(respondent ~ agegroup + sex + maritalstatus + nationality +
             householdsize, family = binomial, data = smpl)

# Latex table of the non-response model coefficients
if (F) {
  require(xtable)
  t <- exp(cbind(OR = coef(nrm), confint.default(nrm)))[-1, ]
  rownames(t) <- do.call(base:::c, lapply(nrm$xlevels, function(x) x[-1]))
  v <- c("age group", "sex", "marital status", "nationality", "household size")
  c <- "Coefficients of the non-response model."
  l <- "tab:nrm_coef"
  x <- xtable(t, caption = c, label = l, digits = c(0, 2, 2, 2))
  align(x) <- c("l", "c", "c", "c")
  h <- "\\toprule \n & Odds ratio & 2.5\\% & 97.5\\% \\\\\n"
  print(x, include.rownames = T, include.colnames = F, hline.after=NULL,
        table.placement = "t",
        add.to.row = list(
          pos = as.list(c(-1, c(0, cumsum(sapply(nrm$xlevels, length) - 1)))),
          command = c(h,
                      paste0("\\midrule \n \\multicolumn{4}{l}{\\it ", v,
                             "} \\\\\n"), 
                      "\\bottomrule \n")), 
        file = file.path("~/nrm_coef.tex"))
  rm(t, v, c, l, x, h)
}

# Response probability (used as a score to determine the response probability
# class)
smpl$rps <- predict(nrm, type = "response")

# Determine the number of response probability classes
if (F) {
  pdf(file.path("~/nrm_kmeans.pdf"))
  plot(1:10, sapply(1:10, function(i) {
    set.seed(666)
    fit <- kmeans(smpl$rps, centers = i, iter.max = 50)
    fit$betweenss / fit$totss
  }), type = "b", xlab = "Number of response probability classes", 
  ylab = "Between sum of squares / total sum of squares", cex.lab = 1.3)
  dev.off()
}

set.seed(666)
smpl$rpc <- fitted(kmeans(smpl$rps, 4))[, 1]

# Weight after non-response correction
smpl$w_nr <- ifelse(smpl$respondent == 1, smpl$w0 / smpl$rpc, NA)

# ------------------------------- Calibration ------------------------------- #

# Margins

##########################3333



margins <- read.table(header = TRUE, text = "
        variable       value        N
        agegroup       15-19    42140
        agegroup       20-39   214498
        agegroup       40-64   270755
        agegroup       65-74    67462
        agegroup        >=75    63605
             sex        Male   320392
             sex      Female   338068
     nationality       Swiss   450498
     nationality     Foreign   207962
   maritalstatus      Single   243386
   maritalstatus     Married   311728
   maritalstatus     Widowed    31110
   maritalstatus    Divorced    69900
   maritalstatus      Others     2336
   householdsize    1_person   130420
   householdsize    2_people   205189
   householdsize    3_people   124702
   householdsize    4_people   131077
   householdsize   5_or_more    67072
")
#margins <- margins[!(margins$variable %in% "agegroup"), ]
tot <- unique(aggregate(N ~ variable, margins, sum)$N)
margins <- do.call(base::c, lapply(unique(margins$variable), function(v) {
  tmp <- margins[margins$variable == v, ][-1, ]
  setNames(tmp$N, paste0(v, tmp$value))
}))
margins <- c(`(Intercept)` = tot, margins)
rm(tot)


#################################



# margins <- c(`(Intercept)` = sum(N$N),
#              setNames(N$N, paste0("strate", N$strate))[-1])

# Calibration weights
# tmp <- subset(smpl, respondent == 1, c("hid", "w_nr", "age"))
# tmp$strate <- factor(tmp$strate)
calw <- svydesign(id = ~hid, weights = ~w_nr,
                  data = subset(smpl, respondent == 1))
calw <- calibrate(calw, ~ agegroup + sex + nationality + maritalstatus + 
                    householdsize, margins, calfun = "raking")
calw <- data.frame(hid = calw$variables$hid, w_nr_rak = 1 / calw$prob)
smpl <- merge(smpl, calw, by = "hid", all.x = TRUE)
rm(calw)


# calage direct - raking

calw <- svydesign(id = ~hid, weights = ~w_naive,
                  data = subset(smpl, respondent == 1))
calw <- calibrate(calw, ~ agegroup + sex + nationality + maritalstatus + 
                    householdsize, margins, calfun = "raking")
calw <- data.frame(hid = calw$variables$hid, w_rak = 1 / calw$prob)
smpl <- merge(smpl, calw, by = "hid", all.x = TRUE)
rm(calw)

plot(w_rak ~ w_nr_rak, smpl)


# calage direct - linear

calw <- svydesign(id = ~hid, weights = ~w_naive,
                  data = subset(smpl, respondent == 1))
calw <- calibrate(calw, ~ agegroup + sex + nationality + maritalstatus + 
                    householdsize, margins, calfun = "linear")
calw <- data.frame(hid = calw$variables$hid, w_lin = 1 / calw$prob)
smpl <- merge(smpl, calw, by = "hid", all.x = TRUE)
rm(calw)


# --------

with(subset(smpl, respondent == 1), summary(w0))
with(subset(smpl, respondent == 1), summary(w_nr))
with(subset(smpl, respondent == 1), summary(w_nr_rak))
with(smpl, sum(w0))
with(subset(smpl, respondent == 1), sum(w_nr))
with(subset(smpl, respondent == 1), sum(w_nr_rak))

if (F) {
  x11(); hist(smpl$w0, breaks = 50)
  x11(); hist(smpl$w_nr, breaks = 50)
  x11(); hist(smpl$w_nr_rak, breaks = 50)
}

# ---------

# Serology
if (FALSE) with(smpl, table(uc_labo_coviggl_v2, uc_labo_covigal_v2))
for (k in 1:2) {
  u <- c("uc_labo_coviggl_v2", "uc_labo_covigal_v2")[k]
  v <- c("serol_igg", "serol_iga")[k]
  smpl[[v]] <- NA
  i <- !is.na(smpl[[u]]) & grepl("(Négatif|Ind(e|é)terminé)", smpl[[u]])
  smpl[i, v] <- 0
  i <- !is.na(smpl[[u]]) & grepl("Positif", smpl[[u]])
  smpl[i, v] <- 1
  if (all(is.na(smpl[[v]]) | smpl[[v]] %in% c("0", "1"))) {
    smpl[[v]] <- as.numeric(smpl[[v]])
  } else {
    stop("error in serology recoding")
  }
}
smpl$serol_any <- pmin(smpl$serol_igg + smpl$serol_iga, 1)
rm(k, i, u, v)

# ---

X <- c("naive", "nr_rak", "rak", "lin")
lapply(X, function(x) {
  smpl$w_current <- smpl[[paste0("w_", x)]]
  d <- svydesign(
    id = ~hid,
    strata = ~strate,
    weights = ~w_current,
    data = merge(subset(smpl, respondent == 1), N, by = "strate"),
    fpc = ~N
  )
  prev0 <- svymean(~serol_any, d)
  #confint(prev, df = sum(wgt$n) - nrow(wgt)) # loi de Student au lieu de la loi normale
  prev1 <- svyby(~serol_any, ~strate, d, svymean)
  confint(prev1)
  prev <- cbind(
    rbind(c(strate = 0, prev0[1]), prev1[, 1:2]),
    rbind(confint(prev0), confint(prev1))
  )
  names(prev)[2:4] <- paste0(x, c("", "_lwr", "_upr"))
  return(prev)
})

plot(w_lin ~ w_nr_rak, smpl)


d0 <- svydesign(id = ~hid, weights = ~w_nr, strata=~strate, fpc=~N,
                data = merge(subset(smpl, respondent == 1), N, by = "strate"))
margins <- margins[!grepl("^agegroup", names(margins))]
d1 <- calibrate(d0, ~ sex + nationality + maritalstatus + 
                  householdsize, margins, calfun = "linear")
(prev0 <- svymean(~serol_any, d1))
confint(prev0)


v <- function(y, h, N, dta) {
  fml <- as.formula(paste(y, "~", h))
  s <- aggregate(fml, dta, function(u) sum((u - mean(u))^2) / (length(u) - 1))
  names(s) <- c("h", "s")
  n <- aggregate(fml, dta, length)
  names(n) <- c("h", "n")
  N <- setNames(unique(dta[c(h, N)]), c("h", "N"))
  dta2 <- merge(N, merge(n, s, by = "h"), by = "h")
  dta2$f <- dta2$n / dta2$N
  with(dta2, sum(N^2 * (1 - f) / n * s) / (sum(N))^2)
}
tmp <- do.call(rbind, lapply(unique(smpl$strate), function(h) {
  tmp <- smpl[smpl$strate == h & smpl$respondent == 1, ]
  X <- c("sex", "nationality",  "maritalstatus", "householdsize")
  for (x in X) tmp[[x]] <- droplevels(tmp[[x]])
  X <- X[sapply(tmp[X], nlevels) > 1]
  fml <- as.formula(paste("serol_any ~", paste(X, collapse = " + ")))
  m <- lm(fml, tmp)
  data.frame(hid = tmp$hid, residual = residuals(m), strate = h)
}))
tmp <- merge(tmp, N, by = "strate")
tmp <- merge(by = "strate", tmp,
             setNames(aggregate(hid ~ strate, tmp, length), c("strate", "n")))
tmp <- merge(tmp, data.frame(hid = d1$variables$hid, w_lin = 1 / d1$prob),
             by = "hid")
tmp$g <- with(tmp, n / N * w_lin)
tmp$e <- tmp$residual * tmp$g
sqrt(v("e", "strate", "N", tmp))


m <- aggregate(tmp$serol_any, list(tmp$strate), function(z) {
 
})
m <- setNames(m[[2]], m[[1]])[tmp$strate]




