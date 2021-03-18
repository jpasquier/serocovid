library(readxl)
library(ISOweek)
library(ggplot2)
library(dplyr)

setwd("~/Projects/SerocoViD")

dta <- list(
  pcrpos = "data-misc/dgs/Cas par semaine, cl. d'âge et sexe.xlsx",
  prev   = "results/prev_serocovid_20210218.xlsx"
)
dta <- lapply(dta, read_xlsx)
dta <- lapply(dta, as.data.frame)
f <- "Tableau statistique des doses par centre-data-15-03-2021 15 08 59.csv"
dta$vacc <- read.csv(file.path("data-misc/vacovid", f))
rm(f)

# La définition de « semaine » est celle que nous utilisons dans EPICOVID
# (première semaine avec au moins quatre jours dans la même année civile, du
# lundi au dimanche).

pcrpos <- aggregate(count ~ year + week, dta$pcrpos, sum)
pcrpos$date <- with(pcrpos, paste0(year, "-W", sprintf("%02d", week), "-4"))
pcrpos$date <- ISOweek2date(pcrpos$date)
pcrpos <- pcrpos[order(pcrpos$date), ]

visits <- read.table(header = TRUE, text = "
  visit         date
      1   2020-06-07
      2   2020-11-15
      3   2021-02-03
")
visits$date <- as.Date(visits$date)


prev <- subset(dta$prev, antibody == "IgG ou IgA" & domain == "Tous (15+)" &
                           visit %in% 1:3)
prev <- merge(prev, visits, by = "visit")

vacc <- dta$vacc
vacc$date <- ISOweek2date(paste0(ISOweek(as.Date(vacc$Date.injection)), "-4"))
vacc <- aggregate(Nombre.de.doses ~ date + Dose.1.ou.2, vacc, sum)
names(vacc)[names(vacc) == "Dose.1.ou.2"] <- "dose"
names(vacc)[names(vacc) == "Nombre.de.doses"] <- "n"
vacc$dose <- sub("Dose ", "", vacc$dose)
vacc <- vacc[vacc$date <= as.Date("2021-02-04"), ]
vacc <- vacc[order(vacc$date), ]
Merge <- function(x, y) merge(x, y, by = "date", all = TRUE)
cvacc <- Reduce(Merge, lapply(1:2, function(k) {
  v <- subset(vacc, dose == k)
  v[[paste0("N", k)]] <- cumsum(v$n)
  v[grep("^(date|N[1-2])", names(v))]
}))
cvacc[is.na(cvacc$N2), "N2"] <- 0

tmp <- rbind(
  cbind(type = 1, pcrpos[c("date", "count")]),
  cbind(type = 2, setNames(cvacc[c("date", "N1")], c("date", "count")))
)
tmp$type <- factor(tmp$type, 1:2, 
  c("Positive PCR/TDR tests", "Vaccinated people (one or two doses)"))

k <- 2 * 10^5

lgd <- c("Positive PCR/TDR tests",
         "Cumulative number of vaccinated people (one or two doses)")
clr <- setNames(c("red", "steelblue"), lgd)


fig <- ggplot(prev, aes(x = date, y = ppos)) +
  geom_point(size = 5) +
  geom_errorbar(aes(ymin = pmax(0, lwr), ymax = upr), size = 1.2, width = 0) +
  geom_bar(aes(y = count / k, fill = names(clr)[1]),
           data = pcrpos, stat = "identity",
           alpha = 0.7) +
  geom_bar(aes(y = N1 / k, fill = names(clr)[2]),
           data = cvacc, stat = "identity", alpha = 0.7) +
  scale_fill_manual(values = clr) +
  #geom_bar(aes(y = N2 / k), data = cvacc, stat = "identity",
  #         fill = "blue", alpha = 0.7) +
  scale_y_continuous(labels = scales::percent,
                     sec.axis = sec_axis(trans =~ . * k)) +
  labs(x = "", y = "Prevalence", fill = "") +
  theme(legend.position = "bottom")


jpeg("~/fig_evol_prev_test_vacc.jpg", height = 3600, width = 7200, res = 512)
print(fig)
dev.off()



  
