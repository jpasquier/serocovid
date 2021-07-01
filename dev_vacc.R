library(REDCapR)
library(survey)
library(writexl)

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

# ------------------------ Research data - Survey 3 ------------------------- #

# Personal data - type of participant and date of birth
vars <- c("uc_s_participant_hid", "uc_s_type_participant")
pdta <- unique(serocovid_data$personal_data[vars])
if (any(duplicated(pdta$uc_s_participant_hid))) {
  stop("Duplicated hid")
}
vacv <- grep("bl_vac", names(serocovid_data$corona_immunitas), value = TRUE)
vacv <- vacv[sapply(serocovid_data$corona_immunitas[vacv],
             function(x) class(x)[1] == "numeric")]
dta <- serocovid_data$corona_immunitas[c("hid", vacv)]
dta <- merge(dta, pdta, by.x = "hid", by.y = "uc_s_participant_hid",
               all.x = TRUE)
if (any(is.na(dta$uc_s_type_participant))) {
  stop("missing type of participant")
}
rm(pdta, vars)

# Select participants of survey 3
dta <- dta[dta$uc_s_type_participant %in% 13 & !grepl("^123123", dta$hid), ]
dta <- dta[!is.na(dta$bl_vac_yn), , drop = FALSE]

# Add strata
dta <- merge(dta, read.csv("data-fso/link_fso_redcap_sample2.csv", sep = ";"),
             by = "hid", all.x = TRUE)
if (any(is.na(dta$PersonId)) | any(duplicated(dta$PersonId))) {
  stop("missing or duplicated fso id")
}
fso <- read.csv("data-fso/COVID19_VD_V2_Total.csv", sep = ";")
dta <- merge(dta, fso[c("PersonId", "strate")], by = "PersonId", all.x = TRUE)
#dta$strate_v4 <- dta$strate - 1
#names(dta)[names(dta) == "strate"] <- "strate_v3"
rm(fso)

# ------------------------------ Raw weigthing ------------------------------ #

# Population totals - survey 3
pop3 <- read.table(header = TRUE, text = "
  strate       N
       1   42140
       2  214498
       3  270755
       4   67462
       5   63605
")

# Weighting with the survey package
wgt <- aggregate(hid ~ strate, dta, length)
names(wgt)[2] <- "n"
wgt <- merge(wgt, pop3, by = "strate")
wgt$wgt <- wgt$N / wgt$n
d <- svydesign(
  id = ~hid,
  strata = ~strate,
  weights = ~wgt,
  data = merge(dta, wgt, by = "strate"),
  fpc = ~N
)
rm(wgt)

# ------------------------------- Proportions ------------------------------- #

# Proportion of missing values
miss_tbl <- t(sapply(vacv, function(v) c(
  mean(is.na(dta[[v]])),
  svymean(as.formula(paste0("~I(is.na(", v, "))")), d)[[2]]
)))
colnames(miss_tbl) <- c("smpl", "pop")

# Proportions
prop_tbl <- do.call(rbind, lapply(vacv, function(v) {
  if (length(unique(na.omit(d$variables[[v]]))) == 1) return(NULL)
  d$variables[[v]] <- factor(d$variables[[v]])
  pop_prop <- svymean(reformulate(v), d, na.rm = TRUE)
  rwnm <- sub(v, paste0(v, " "), names(pop_prop))
  tbl <- cbind(smpl_prop = prop.table(table(d$variables[[v]])),
               pop_prop, confint(pop_prop))
  rownames(tbl) <- rwnm
  return(tbl)
}))

# ----------------------------- Export results ------------------------------ #

L <- list(
  proportions = cbind(data.frame(variable = rownames(prop_tbl)), prop_tbl),
  `valeurs manquantes` = 
     cbind(data.frame(variable = rownames(miss_tbl)), miss_tbl)
)
write_xlsx(L, path = "results/bl_vac_proportions_20210520.xlsx")
rm(L)
