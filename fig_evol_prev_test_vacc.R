library(ISOweek)
library(dplyr)
library(ggplot2)
library(ggpubr)
library(grid)
library(readxl)
library(scales)
library(tidyr)
library(wesanderson)
library(writexl)

# Working directory
setwd("~/Projects/SerocoViD")

# Age groups
age_grps <- c("6m-4", "5-9", "10-14", "15-19",
              "20-39", "40-64", "65-74", "75+")

# Pop VD (stat VD)
N_pop <- read_xlsx("data-misc/stat_vd/pop_2021_stat_vd.xlsx") %>%
  mutate(Age = as.numeric(sub("\\+? ans?$", "", Age)),
         age_group = cut(Age, c(0, 4, 9, 14, 19, 39, 64, 74, Inf),
                         labels = age_grps, include.lowest = TRUE)) %>%
  group_by(age_group) %>%
  summarise(N_pop = sum(Total))

# Positive tests
# REMARQUE DE LA DGS : La définition de « semaine » est celle que nous
# utilisons dans EPICOVID (première semaine avec au moins quatre jours dans la
# même année civile, du lundi au dimanche).
pcrpos <- local({
  z <- "data_cases_weekly_gender_agegr_serocovid";
  f <- paste0("data-misc/dgs/", z, ".RData");
  load(f);
  get(z)
}) %>%
  filter(!(age_group %in% c("0-5m", "Inconnu"))) %>%
  group_by(year, week, age_group) %>%
  summarise(count = sum(count, na.rm = TRUE), .groups = "drop") %>%
  mutate(
    age_group = factor(sub("\\s?ans", "", age_group), age_grps),
    date = ISOweek2date(paste0(year, "-W", sprintf("%02d", week), "-4"))
  ) %>%
  left_join(N_pop, by = "age_group") %>%
  mutate(incidence = count * 10^5 / N_pop)

# Median date of surveys
visits <- read.table(header = TRUE, text = "
  visit         date
      1   2020-06-07
      2   2020-11-15
      3   2021-02-03
      4   2021-10-01
")
visits$date <- as.Date(visits$date)

# Prevalence (proportion of positive people)
prev <- readRDS("results/prev_serocovid_20220509.rds") %>%
  filter(antibody == "IgG", domain == "stratum",
         !(visit == 2 & value == "10-14")) %>%
  rename(prev = ppos) %>%
  mutate(age_group = factor(sub(">=75", "75+", value), age_grps)) %>%
  left_join(visits, by = "visit")

# Number of vaccinated people (source: VaCoViD)
sem <- expand.grid(annee_injection = 2020:2022, semaine_injection = 1:53,
                   Group_Age = c("4 ans et moins", "Entre 5 et 9 ans",
                                 "Entre 10 et 14 ans", "Entre 15 et 19 ans",
                                 "Entre 20 et 39 ans", "Entre 40 et 64 ans",
                                 "Entre 65 et 74 ans", "75 ans et +")) %>%
  filter(annee_injection == 2020 | semaine_injection %in% 1:52) %>%
  filter(annee_injection %in% 2020:2021 | semaine_injection <= 14)
vacc <- read_xlsx("data-misc/vacovid/Statistiques_premieres_doses.xlsx",
                  range = cell_cols("A:D")) %>%
  full_join(sem, by = c("annee_injection", "semaine_injection",
                        "Group_Age")) %>%
  mutate(
    Count_all = ifelse(is.na(Count_all), 0, Count_all),
    age_group = sub("et \\+", "+", gsub("Entre|ans", "", Group_Age)),
    age_group = sub("et", "-", sub("4 +et moins", "6m-4", age_group)),
    age_group = factor(gsub("\\s+", "", age_group), age_grps),
    date = ISOweek2date(paste0(annee_injection, "-W",
                               sprintf("%02d", semaine_injection), "-4"))
  ) %>%
  rename(n1 = Count_all) %>%
  arrange(age_group, date) %>%
  group_by(age_group) %>%
  mutate(N1 = cumsum(n1)) %>%
  left_join(N_pop, by = "age_group") %>%
  mutate(P1 = N1 / N_pop)

# Maximum number of cases per age class
max_cases <- full_join(pcrpos, vacc, by = c("age_group", "date")) %>%
  mutate(max_cases = pmax(count, n1, na.rm = TRUE)) %>%
  group_by(age_group) %>%
  summarize(max_cases = max(max_cases))

# Figure
k <- 5600
lgd <- c("Seroprevalence",
         "Incidence of positive PCR/TDR cases (per 10'000 persons)",
         "Proportion of vaccinated people (one dose or more)")
clr <- setNames(c("black", wes_palette(n = 2, name = "Darjeeling1")), lgd)
fig1 <- ggplot(prev, aes(x = date)) +
  geom_point(aes(y = prev), size = 1) +
  geom_errorbar(aes(ymin = pmax(0, lwr), ymax = upr), size = 1.2, width = 0) +
  geom_bar(aes(y = P1, fill = names(clr)[3]),
           data = vacc, stat = "identity", alpha = 0.5) +
  geom_bar(aes(y = incidence / k, fill = names(clr)[2]),
           data = pcrpos, stat = "identity",
           alpha = 0.7) +
  scale_fill_manual(values = clr) +
  scale_y_continuous(labels = percent,
                     sec.axis = sec_axis(trans =~ . * k,
                                         name = "Incidence")) +
  facet_grid(rows = vars(age_group)) +
  labs(x = "", y = "", fill = "") +
  theme(legend.position = "bottom", legend.text = element_text(size=rel(1.2)),
        axis.title = element_text(size = rel(1.2)),
        axis.text = element_text(size = rel(1.2))) +
  guides(fill = guide_legend(override.aes = list(alpha = 0.5)))
fig1

# Figure - With a different scale on each facet for the second axis
fig2 <- lapply(age_grps, function(a) {
  k <- filter(max_cases, age_group == a)$max_cases * 1.1
  fig <- ggplot(filter(prev, age_group == a), aes(x = date, y = prev)) +
    geom_point(size = 1) +
    geom_errorbar(aes(ymin = pmax(0, lwr), ymax = upr), size = 1.2,
                  width = 0) +
    geom_bar(aes(y = count / k, fill = names(clr)[2]),
             data = filter(pcrpos, age_group == a), stat = "identity",
             alpha = 0.7) +
    geom_bar(aes(y = n1 / k, fill = names(clr)[3]),
             data = filter(vacc, age_group == a),
             stat = "identity", alpha = 0.5) +
    scale_fill_manual(values = clr[-1]) +
    scale_y_continuous(labels = percent,
                       sec.axis = sec_axis(trans =~ . * k)) +
    facet_grid(rows = vars(age_group)) +
    labs(x = "", y = "Prevalence", fill = "") +
    theme(axis.title = element_blank(),
          legend.text = element_text(size=rel(1.0))) +
    guides(fill = guide_legend(override.aes = list(alpha = 0.5)))
  if (a != age_grps[length(age_grps)]) {
    fig <- fig +
      theme(axis.text.x = element_blank(),
            axis.ticks.x = element_blank())
  }
  return(fig)
}) %>%
  ggarrange(plotlist = ., ncol = 1, align = "v", legend = "bottom",
            common.legend = TRUE) %>%
  annotate_figure(
    left = textGrob("Prevalence", rot = 90, vjust = 1, gp = gpar(cex = 1.3)),
    right = textGrob("Number of cases", rot = 270, vjust = 1,
                     gp = gpar(cex = 1.3))
  )
fig2
jpeg("results/fig_evol_prev_test_vacc_DEV.jpg", height = 7200,
     width = 3600, res = 256)
print(fig2)
dev.off()

# Export data
if (FALSE) {
write_xlsx(list(prev = prev, pcrpos = pcrpos, vacc = vacc),
           "results/fig_evol_prev_test_vacc_data_DEV.xlsx")
}
