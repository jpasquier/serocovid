library(ISOweek)
library(dplyr)
library(ggplot2)
library(ggrepel)
library(readxl)
library(scales)
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
  summarise(N_pop = sum(Total)) %>%
  bind_rows(data.frame(
    age_group = c("all", ">=15"), 
    N_pop = c(sum(.$N_pop), sum(.[.$age_group %in% age_grps[4:8], "N_pop"]))
  )) %>%
  mutate(age_group = factor(age_group, c("all", ">=15", age_grps)))

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
  filter(!(year == 2022 & week > 30)) %>%
  filter(!(age_group %in% c("0-5m", "Inconnu"))) %>%
  group_by(year, week, age_group) %>%
  summarise(count = sum(count, na.rm = TRUE), .groups = "drop") %>%
  mutate(
    age_group = factor(sub("\\s?ans", "", age_group), age_grps),
    date = ISOweek2date(paste0(year, "-W", sprintf("%02d", week), "-4"))
  ) %>%
  bind_rows(
    group_by(., year, week, date) %>%
      summarise(count = sum(count), .groups = "drop") %>%
      cbind(age_group = "all"),
    filter(., age_group %in% age_grps[4:8]) %>%
      group_by(year, week, date) %>%
      summarise(count = sum(count), .groups = "drop") %>%
      cbind(age_group = ">=15")
  ) %>%
  left_join(N_pop, by = "age_group") %>%
  mutate(
    age_group = factor(age_group, c("all", ">=15", age_grps)),
    incidence = count * 10^5 / N_pop
  )

# Median date of surveys
visits <- read.table(header = TRUE, text = "
  visit         date
      1   2020-06-07
      2   2020-11-15
      3   2021-02-03
      4   2021-10-01
      6   2022-06-15
")
visits$date <- as.Date(visits$date)

# Prevalence (proportion of positive people)
prev <- readRDS("results/prev_serocovid_20220929.rds") %>%
  filter(
    antibody == "IgG",
    domain == "stratum" |
      domain == "p15" & value == ">=15" |
      visit %in% c(1, 4) & domain == "all" |
      visit == 6 & domain == "age_group",
    !(visit == 2 & value == "10-14"),
    !(visit == 6 & domain == "stratum")
  ) %>%
  rename(prev = ppos) %>%
  mutate(
    age_group = tolower(sub(">=75", "75+", value)), 
    age_group = factor(age_group, c("all", ">=15", age_grps))
  ) %>%
  left_join(visits, by = "visit")

# Number of vaccinated people (source: VaCoViD)
vacc <- read_xlsx(file.path(
  'data-misc/vacovid',
  '2022-09-13_vacovid_firstdose_weekly_agegr_serocovid.xlsx'
)) %>%
  mutate(
    year = as.integer(year),
    week = as.integer(week),
    agegr_serocovid = sub("0-4", "6m-4", sub(" ans", "", agegr_serocovid)),
    N = as.integer(N)
  ) %>%
  rename(age_group = agegr_serocovid, n1 = N) %>%
  filter(n1 != 0) %>%
  right_join(
    expand.grid(year = 2020L:2022L, week = 1L:53L, age_group = age_grps) %>%
      filter(year == 2020L | week %in% 1L:52L) %>%
      filter(year %in% 2020L:2021L | week <= 30L),
    by = c("year", "week", "age_group")
  ) %>%
  mutate(n1 = if_else(is.na(n1), 0L, n1)) %>%
  bind_rows(
    group_by(., year, week) %>%
      summarise(n1 = sum(n1), .groups = "drop") %>%
      cbind(age_group = "all"),
    filter(., age_group %in% age_grps[4:8]) %>%
      group_by(year, week) %>%
      summarise(n1 = sum(n1), .groups = "drop") %>%
      cbind(age_group = ">=15")
  ) %>%
  mutate(
    age_group = factor(age_group, c("all", ">=15", age_grps)),
    date = ISOweek2date(paste0(year, "-W", sprintf("%02d", week), "-4"))
  ) %>%
  arrange(age_group, date) %>%
  group_by(age_group) %>%
  mutate(N1 = cumsum(n1)) %>%
  ungroup() %>%
  left_join(N_pop, by = "age_group") %>%
  mutate(P1 = N1 / N_pop)

# Figure
k <- ceiling(max(pcrpos$incidence) * 1.01)
lgd <- c("Seroprevalence",
         "Incidence of positive PCR/TDR cases (per 10'000 persons)",
         "Proportion of vaccinated people (one dose or more)")
clr <- setNames(c("dodgerblue4", wes_palette(n = 2, name = "Darjeeling1")),
                lgd)
fig <- ggplot(prev, aes(x = date)) +
  geom_bar(aes(y = incidence / k, fill = names(clr)[2]),
           data = pcrpos, stat = "identity") +
  geom_line(aes(y = P1), data = vacc, color = clr[3],
            linetype = "dashed") +
  geom_errorbar(aes(ymin = pmax(0, lwr), ymax = upr), linewidth = 1.2,
                width = 0, color = clr[1]) +
  geom_line(aes(y = prev), linetype = "dashed", color = clr[1]) +
  geom_point(aes(y = prev), size = 1, color = clr[1]) +
  geom_text(aes(y = prev, label = paste(round(prev * 100, 1), "%"),
                vjust = ifelse(prev > .5, 2.5, -2.5)),
            hjust = .5, color = clr[1]) +
  scale_fill_manual(values = clr) +
  scale_y_continuous(labels = percent,
                     sec.axis = sec_axis(trans =~ . * k,
                                         name = "Incidence")) +
  facet_grid(rows = vars(age_group)) +
  labs(x = "", y = "", fill = "") +
  theme_bw() +
  theme(legend.position = "bottom",
        legend.text = element_text(size = rel(1.2)),
        axis.title = element_text(size = rel(1.2)),
        axis.text = element_text(size = rel(1.2)))
if (TRUE) {
jpeg(paste0("results/fig_evol_prev_test_vacc_",
            format(Sys.Date(), "%Y%m%d"), ".jpg"),
     height = 7200, width = 3600, res = 256)
print(fig)
dev.off()
}

# Export data
if (FALSE) {
write_xlsx(list(prev = prev, pcrpos = pcrpos, vacc = vacc),
           paste0("results/fig_evol_prev_test_vacc_data_",
                  format(Sys.Date(), "%Y%m%d"), ".xlsx"))
}
