library(dplyr)
library(ggplot2)
library(gridExtra)
library(survey)

# Working directory
setwd("~/Projects/SerocoViD")

# Prepare data
source("R/prepare_data_4th_visit.R")

# Compare original and standardized values
dta <- serol %>%
  select(serol, uc_labo_coviggl_spe_v4, uc_labo_coviggr_v4,
         uc_labo_coviggr_zscore_v4) %>%
  na.omit() %>%
  mutate(uc_labo_coviggl_spe_v4 =
    factor(uc_labo_coviggl_spe_v4, c("négatif", "faible", "moyen", "élevé")))
p1 <- ggplot(dta, aes(x = uc_labo_coviggl_spe_v4,
                      y = uc_labo_coviggr_v4)) +
  geom_boxplot() +
  labs(title = "Sans normalisation", caption = paste("N =", nrow(dta)))
p2 <- ggplot(dta, aes(x = uc_labo_coviggl_spe_v4,
                      y = uc_labo_coviggr_zscore_v4)) +
  geom_boxplot() +
  labs(title = "Avec normalistation", caption = paste("N =", nrow(dta)))
tiff("results/antibody_score_4th_visit_standardization.tiff",
     height = 3600, width = 7200, res = 768, compression = "zip")
grid.arrange(p1, p2, nrow = 1)
dev.off()
rm(dta, p1, p2)

# Boxplots
dta <- smpl %>%
  group_by(strate) %>%
  summarise(resprate = mean(respvac)) %>%
  inner_join(filter(smpl, respvac == 1), by = "strate") %>%
  mutate(
    w = wsmpl / resprate,
    strate = factor(strate, 1:8, c("6m-4", "5-9", "10-14", "15-19", "20-39",
                                   "40-64", "65-74", ">=75")),
    vac = factor(vac, c("y", "n"), c("Vaccinated", "Not vaccinated"))
  ) %>%
  filter(serol == 1)
figs <- lapply(1:4, function(k) {
  x <- if (k %% 2 == 1) "vac" else "strate"
  n_fun <- function(x) data.frame(y = c(-20, 0.3)[(k + 1) %/% 2], 
                                  label = paste0("(", length(x), ")"))
  p <- ggplot(dta, aes_string(x = x, y = "serol_num")) +
    stat_summary(fun.data = n_fun, geom = "text", hjust = 0.5, size = 1.8)
  y_lab <- "IgG score"
  ttl <- "Antibody score"
  sttl <- "By vaccination status"
  cap <- paste("N =", nrow(dta))
  m <- list(filename = "results/antibody_score_4th_visit_by_vac.tiff",
            height = 3600, width = 5400, res = 1024, compression = "zip")
  if (k %% 2 == 1) {
    p <- p + geom_boxplot(aes(weight = w))
    cap <- paste("Weighted data,", cap)
  } else {
    p <- p + 
      geom_boxplot() +
      facet_grid(rows = vars(vac), scales = "free")
    sttl <- paste(sttl, "and age")
    m$filename <- sub("by_vac", "by_vac_and_age", m$filename)
    m$width <- 7200
    m$res <- 768
  }
  if (k >= 3) {
    b <- c(3, 10, 30, 100, 300)
    p <- p + scale_y_continuous(breaks = b, trans = "log10", labels = b)
    y_lab <- paste(y_lab, "(log scale)")
    m$filename <- sub("score", "log_score", m$filename)
  }
  p <- p + labs(x = "", y = y_lab, subtitle = sttl, title = ttl, caption = cap)
  attr(p, "metadata") <- m
  return(p)
})
for (fig in figs) {
  do.call(tiff, attr(fig, "metadata"))
  print(fig)
  dev.off()
}

# Weighted quantiles
# https://stats.stackexchange.com/questions/13169
wgt.quantile <- function(x, w, p)  {
  w <- w[order(x)]
  x <- sort(x)
  n <- length(w)
  s <- 0:(n - 1) * w + (n - 1) * c(0, cumsum(w)[-n])
  k <- max(which(s / s[n] <= p))
  x[k] + (x[k+1] - x[k]) * (p * s[n] - s[k]) / (s[k + 1] - s[k])
}
with(dta, wgt.quantile(serol_num, w, .75))

