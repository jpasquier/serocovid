library(readxl)
library(lme4)

# Working directory
setwd("~/Projects/SerocoViD")

# Data
sero <- as.data.frame(read_xlsx("data-raw/data_20200707.xlsx"))
sero <- subset(sero, grepl("même toit", uc_info_type_participant))
sero <- sero[c("sero_result", "age_cat_2", "index_case_age_cat_2",
               "index_case_hid")]
apply(is.na(sero), 2, sum)
sero <- na.omit(sero)

# Recoding
sero$sero_result <- as.integer(sero$sero_result == "positive")
sero$age_cat_2 <- factor(sero$age_cat_2)
sero$index_case_age_cat_2 <- factor(sero$index_case_age_cat_2)

# Fixed effects models
m1 <- glm(sero_result ~ age_cat_2, family = binomial, data = sero)
summary(m1)
exp(coef(m1))  # odds ratios
m2 <- glm(sero_result ~ index_case_age_cat_2, family = binomial, data = sero)
summary(m2)
exp(coef(m2))  # odds ratios
m3 <- glm(sero_result ~ age_cat_2 + index_case_age_cat_2,
          family = binomial, data = sero)
summary(m3)
exp(coef(m3))  # odds ratios

# Mixed effect models
m4 <- glmer(sero_result ~ age_cat_2 + (1 | index_case_hid),
            family = binomial, data = sero)
summary(m4)
exp(fixef(m4))  # odds ratios
m5 <- glmer(sero_result ~ index_case_age_cat_2 + (1 | index_case_hid),
            family = binomial, data = sero)
summary(m5)
exp(fixef(m5))  # odds ratios
m6 <- glmer(sero_result ~ age_cat_2 + index_case_age_cat_2 +
              (1 | index_case_hid),
            family = binomial, data = sero)
summary(m6)
exp(fixef(m6))  # odds ratios

# sero by cluster
u <- aggregate(sero_result ~ index_case_hid, sero, mean)[, 2]
table(ifelse(u == 0, "[0]", ifelse(u == 1, "[1]", "(0,1)")))


