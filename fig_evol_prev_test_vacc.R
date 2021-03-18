library(readxl)
library(ISOweek)
library(ggplot2)
library(dplyr)
library(tidyr)
library(wesanderson)
library(writexl)

# Working directory
setwd("~/Projects/SerocoViD")

# Import data
dta <- list(
  pcrpos = "data-misc/dgs/Cas par semaine, cl. d'âge et sexe.xlsx",
  prev   = "results/prev_serocovid_20210318.xlsx"
)
dta <- lapply(dta, read_xlsx)
dta <- lapply(dta, as.data.frame)
f <- "Tableau statistique des doses par centre-data-15-03-2021 15 08 59.csv"
dta$vacc <- read.csv(file.path("data-misc/vacovid", f))
rm(f)

# REMARQUE DE LA DGS : La définition de « semaine » est celle que nous
# utilisons dans EPICOVID (première semaine avec au moins quatre jours dans la
# même année civile, du lundi au dimanche).
pcrpos <- aggregate(count ~ year + week, dta$pcrpos, sum) %>%
  arrange(year, week) %>%
  mutate(cumcount = cumsum(count),
         date = paste0(year, "-W", sprintf("%02d", week), "-4"),
         date = ISOweek2date(date))

# Median date of surveys
visits <- read.table(header = TRUE, text = "
  visit         date
      1   2020-06-07
      2   2020-11-15
      3   2021-02-03
")
visits$date <- as.Date(visits$date)

# Prevalence (proportion of positive people)
prev <- subset(dta$prev, antibody == "IgG or IgA" & domain == "All (15+)" &
                           visit %in% 1:3)
prev <- merge(prev, visits, by = "visit")

# Number of vaccinated people (source: VaCoViD)
vacc <- dta$vacc
vacc$date <- ISOweek2date(paste0(ISOweek(as.Date(vacc$Date.injection)), "-4"))
vacc <- aggregate(Nombre.de.doses ~ date + Dose.1.ou.2, vacc, sum) %>%
  rename(dose = Dose.1.ou.2, n = Nombre.de.doses) %>%
  mutate(dose = sub("Dose ", "", dose)) %>%
  filter(date <= as.Date("2021-02-04")) %>%
  pivot_wider(names_from = dose, values_from = n, names_prefix = "n",
              values_fill = 0) %>%
  arrange(date) %>%
  mutate(N1 = cumsum(n1), N2 = cumsum(n2))

# Figure
k <- 2 * 10^5
lgd <- c("Positive PCR/TDR tests",
         "Cumulative number of vaccinated people (one or two doses)")
clr <- setNames(wes_palette(n = 2, name = "Darjeeling1"), lgd)
fig <- ggplot(prev, aes(x = date, y = ppos)) +
  geom_point(size = 5) +
  geom_errorbar(aes(ymin = pmax(0, lwr), ymax = upr), size = 1.2, width = 0) +
  geom_bar(aes(y = count / k, fill = names(clr)[1]),
           data = pcrpos, stat = "identity",
           alpha = 0.7) +
  geom_bar(aes(y = N1 / k, fill = names(clr)[2]),
           data = vacc, stat = "identity", alpha = 0.7) +
  scale_fill_manual(values = clr) +
  scale_y_continuous(labels = scales::percent,
                     sec.axis = sec_axis(trans =~ . * k)) +
  labs(x = "", y = "Prevalence", fill = "") +
  theme(legend.position = "bottom", legend.text = element_text(size=rel(1.2)),
        axis.title = element_text(size = rel(1.2)),
        axis.text = element_text(size = rel(1.2)))
jpeg("results/fig_evol_prev_test_vacc_20210318.jpg", height = 3600,
     width = 7200, res = 512)
print(fig)
dev.off()

# Export data
write_xlsx(list(prev = prev, pcrpos = pcrpos, vacc = vacc),
           "results/fig_evol_prev_test_vacc_data_20210318.xlsx")
