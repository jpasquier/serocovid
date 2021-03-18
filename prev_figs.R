library(readxl)
library(ggplot2)
setwd("~/Projects/SerocoViD")
prev <- lapply(1:2, function(k) {
  file_name <- paste0("results/prev_", c("1st", "2nd")[k],
                      "_wave_20201124.xlsx")
  tbl <- read_xlsx(file_name, sheet = "strates_4_8")
  tbl <- subset(tbl, Antibody == "IgG")
  tbl$value <- as.numeric(tbl$value)
  j <- which(names(tbl) == "ppos")
  tbl <- tbl[, -(j + (0:2))]
  names(tbl)[grep("^(2\\.5%)", names(tbl))] <- "lwr"
  names(tbl)[grep("^(97\\.5%)", names(tbl))] <- "upr"
  attr(tbl, "wave") <- c("First", "Second")[k]
  tbl
})
strata <- c("0-4", "5-9", "10-14", "15-19", "20-39", "40-64", "65-74", ">=75")
figs <- lapply(prev, function(tbl) {
  tbl <- subset(tbl, domain == "stratum")
  tbl$age_class <- droplevels(factor(strata[tbl$value], strata))
  ggplot(tbl, aes(x = age_class, y = prev)) +
    geom_point() +
    geom_errorbar(aes(ymin = pmax(0, lwr), ymax = upr), width = .2) +
    scale_y_continuous(labels = scales::percent) +
    labs(x = "", y = "", title = "Seroprevalence by age group",
         subtitle = paste(attr(tbl, "wave"), "wave"))
})
figs[[3]] <- ggplot(subset(prev[[2]], domain == "weeknbr"),
                    aes(x = value, y = prev)) +
  geom_point() +
  geom_errorbar(aes(ymin = pmax(0, lwr), ymax = upr), width = .2) +
  scale_y_continuous(labels = scales::percent) +
  labs(x = "", y = "", title = "Seroprevalence by week",
       subtitle = paste(attr(prev[[2]], "wave"), "wave"))
pdf("results/prev_figs_20201125.pdf")
for(fig in figs) print(fig)
dev.off()
