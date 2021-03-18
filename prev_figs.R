library(readxl)
library(ggplot2)
setwd("~/Projects/SerocoViD")
file_name <- paste0("results/prev_2nd_wave_20201218.xlsx")
prev <- read_xlsx(file_name, sheet = "strates_4_8")
prev <- subset(prev, Antibody == "IgG" & grepl("weeknbr", domain))
prev$value <- as.numeric(prev$value)
j <- which(names(prev) == "prev")
prev <- prev[, -(j + (0:2))]
names(prev)[grep("^(2\\.5%)", names(prev))] <- "lwr"
names(prev)[grep("^(97\\.5%)", names(prev))] <- "upr"
figs <- lapply(1:2, function(k) {
  if (k == 1) {
    d <- "weeknbr"
    sttl <- "Strata 4 to 8 - Method 'domain'"
  } else {
    d <- "weeknbr (method 2)"
    sttl <- "Strata 4 to 8 - Method 2"
  }
  fig <- ggplot(subset(prev, domain == d),
                aes(x = as.factor(value), y = ppos)) +
    geom_point() +
    geom_errorbar(aes(ymin = pmax(0, lwr), ymax = upr), width = .2) +
    scale_y_continuous(labels = scales::percent, limits = c(0, 1)) +
    labs(x = "", y = "", title = "Proportion of positive tests by week",
         subtitle = sttl)
  fig
})
pdf("results/prev_figs_20201218.pdf")
for (fig in figs) print(fig)
rm(fig)
dev.off()
