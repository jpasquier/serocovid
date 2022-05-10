library(ISOweek)
library(dplyr)
library(ggplot2)
library(readxl)

# Working directory
setwd("~/Projects/SerocoViD")

# Number of positive cases
pcrpos0 <- "data-misc/dgs/Cas par semaine, cl. d'âge et sexe.xlsx" %>%
  read_xlsx() %>%
  group_by(year, week) %>%
  summarise(count = sum(count), .groups = "drop") %>%
  arrange(year, week) %>%
  mutate(src = "previous data",
         date = paste0(year, "-W", sprintf("%02d", week), "-4"),
         date = ISOweek2date(date))
pcrpos1 <- local({
  z <- "data_cases_weekly_gender_agegr_serocovid";
  f <- paste0("data-misc/dgs/", z, ".RData");
  load(f);
  get(z)
}) %>%
  group_by(year, week) %>%
  summarise(count = sum(count, na.rm = TRUE), .groups = "drop")  %>%
  mutate(date = paste0(year, "-W", sprintf("%02d", week), "-4"),
         date = ISOweek2date(date),
         src = "newdata")
x11()
bind_rows(pcrpos0, pcrpos1) %>%
  filter(date <= as.Date("2021-02-04")) %>%
  ggplot(aes(x = date, y = count, colour = src)) +
  geom_point() +
  geom_line() +
  labs(title = "Number of positive PCR/TDR cases") +
  theme(legend.title = element_blank(), legend.position = "bottom")

# Number of vaccinated people (source: VaCoViD)
vacc0 <- paste("data-misc/vacovid/Tableau statistique des doses par",
               "centre-data-15-03-2021 15 08 59.csv") %>%
  read.csv() %>%
  filter(Dose.1.ou.2 == "Dose 1") %>%
  mutate(date = ISOweek(as.Date(Date.injection)),
         date = ISOweek2date(paste0(date, "-4"))) %>%
  group_by(date) %>%
  summarise(n1 = sum(Nombre.de.doses)) %>%
  mutate(src = "previous data")
vacc1 <- read_xlsx("data-misc/vacovid/Statistiques_premieres_doses.xlsx",
                      range = cell_cols("A:D")) %>%
  mutate(date = paste0(annee_injection, "-W",
                       sprintf("%02d", semaine_injection), "-4"),
         date = ISOweek2date(date)) %>%
  group_by(date) %>%
  summarize(n1 = sum(Count_all)) %>%
  arrange(date) %>%
  mutate(N1 = cumsum(n1), src = "new data")
x11()
bind_rows(vacc0, vacc1) %>%
  filter(date <= as.Date("2021-03-18")) %>%
  ggplot(aes(x = date, y = n1, colour = src)) +
  geom_point() +
  geom_line() +
  labs(title = "Number of vaccinated people (one dose or more)", x = "count") +
  theme(legend.title = element_blank(), legend.position = "bottom")
