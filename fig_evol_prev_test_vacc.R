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
  )

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
vacc <- read_xlsx("data-misc/vacovid/Statistiques_premieres_doses.xlsx",
                  range = cell_cols("A:D")) %>%
  mutate(
    age_group = sub("et \\+", "+", gsub("Entre|ans", "", Group_Age)),
    age_group = sub("et", "-", sub("4 +et moins", "6m-4", age_group)),
    age_group = factor(gsub("\\s+", "", age_group), age_grps),
    date = ISOweek2date(paste0(annee_injection, "-W",
                               sprintf("%02d", semaine_injection), "-4"))
  ) %>%
  rename(n1 = Count_all) %>%
  arrange(age_group, date) %>%
  group_by(age_group) %>%
  mutate(N1 = cumsum(n1))

# Maximum number of cases per age class
max_cases <- full_join(pcrpos, vacc, by = c("age_group", "date")) %>%
  mutate(max_cases = pmax(count, n1, na.rm = TRUE)) %>%
  group_by(age_group) %>%
  summarize(max_cases = max(max_cases))

# Figure
k <- 22000
lgd <- c("Positive PCR/TDR cases",
         "Number of vaccinated people (one dose or more)")
clr <- setNames(wes_palette(n = 2, name = "Darjeeling1"), lgd)
fig1 <- ggplot(prev, aes(x = date, y = prev)) +
  geom_point(size = 1) +
  geom_errorbar(aes(ymin = pmax(0, lwr), ymax = upr), size = 1.2, width = 0) +
  geom_bar(aes(y = count / k, fill = names(clr)[1]),
           data = pcrpos, stat = "identity",
           alpha = 0.7) +
  geom_bar(aes(y = n1 / k, fill = names(clr)[2]),
           data = vacc, stat = "identity", alpha = 0.5) +
  scale_fill_manual(values = clr) +
  scale_y_continuous(labels = percent,
                     sec.axis = sec_axis(trans =~ . * k,
                                         name = "Number of cases")) +
  facet_grid(rows = vars(age_group)) +
  labs(x = "", y = "Prevalence", fill = "") +
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
    geom_bar(aes(y = count / k, fill = names(clr)[1]),
             data = filter(pcrpos, age_group == a), stat = "identity",
             alpha = 0.7) +
    geom_bar(aes(y = n1 / k, fill = names(clr)[2]),
             data = filter(vacc, age_group == a),
             stat = "identity", alpha = 0.5) +
    scale_fill_manual(values = clr) +
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
