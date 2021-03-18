library(REDCapR)
library(survival)
library(ggplot2)
library(survminer)
#library(ggforce)
#library(wesanderson)

# Working directory
setwd("~/Projects/SerocoViD")

# ------------------------------- REDCap data ------------------------------- #

# API
uri   <- "https://redcap.unisante.ch/api/"

# Tokens
tokens <- list(
  corona_immunitas = "F1BFC49C8D21C355AE9946277DD52EDC",
  personal_data    = "2553F65D05A7FB4FB8AAACE5BDE0784C",
  research_data    = "C296C589D5FB1A89DBFC11D14EFA4365"
)

# Import
tmp_file <- "/tmp/serocovid_data.rda"
if (file.exists(tmp_file)) {
  load(tmp_file)
} else {
  serocovid_data <- lapply(tokens, function(token) {
    redcap_read(redcap_uri = uri, token = token)$data
  })
  save(serocovid_data, file = tmp_file)
}
rm(uri, tokens, tmp_file)

# ------------------------ Research data - Survey 3 ------------------------- #

# Personal data - type of participant
vars <- c("uc_s_participant_hid", "uc_s_type_participant")
pdta <- unique(serocovid_data$personal_data[vars])
if (any(duplicated(pdta$uc_s_participant_hid))) {
  stop("Duplicated hid")
}

# Research data
vars <- c("hid", "pcr_date", "pcr_result", "uc_ser_date_hr_kittrans",
          "uc_labo_coviggl_v2", "uc_labo_covigal_v2", "uc_labo_coviggr_v2",
          "uc_labo_covigar_v2", "bl_hosp_start", "bl_vac_yn", "bl_vac_first")
dta3 <- serocovid_data$corona_immunitas[vars]
dta3 <- merge(dta3, pdta, by.x = "hid", by.y = "uc_s_participant_hid",
              all.x = TRUE)
if (any(is.na(dta3$uc_s_type_participant))) {
  stop("missing type of participant")
}
rm(vars, pdta)

# Select observation of survey 3
dta3 <- dta3[dta3$uc_s_type_participant %in% 13, ]

# Remove fake participants from survey 3
dta3 <- dta3[!grepl("^123123", dta3$hid), ]

# Serology
if (FALSE) with(dta3, table(uc_labo_coviggl_v2, uc_labo_covigal_v2))
for (k in 1:2) {
  u <- c("uc_labo_coviggl_v2", "uc_labo_covigal_v2")[k]
  v <- c("serol_igg", "serol_iga")[k]
  dta3[[v]] <- NA
  i <- !is.na(dta3[[u]]) & grepl("(Négatif|Ind(e|é)terminé)", dta3[[u]])
  dta3[i, v] <- 0
  i <- !is.na(dta3[[u]]) & grepl("Positif", dta3[[u]])
  dta3[i, v] <- 1
  if (all(is.na(dta3[[v]]) | dta3[[v]] %in% c("0", "1"))) {
    dta3[[v]] <- as.numeric(dta3[[v]])
  } else {
    stop("error in serology recoding")
  }
}
dta3$serol_any <- pmin(dta3$serol_igg + dta3$serol_iga, 1)
rm(k, i, u, v)

# Select observations for which serology is available
with(dta3,
 all(!is.na(serol_igg) == !is.na(serol_iga)) &
 all(!is.na(serol_igg) == !is.na(uc_labo_coviggr_v2)) &
 all(!is.na(serol_igg) == !is.na(uc_labo_covigar_v2))
)
with(dta3, table(pcr_result %in% 1, !is.na(serol_igg)))
dta3 <- dta3[!is.na(dta3$serol_igg), ]

# Remove people with unknown vaccinated status
# hid 3-37922994: bl_vac_yn is missing, serology is negative
# hid 3-53497922, 3-58277322, 3-82348728 bl_vac:
#   bl_vac_yn==3, serology is negative
dta3[is.na(dta3$bl_vac_yn) | dta3$bl_vac_yn == 3, ]
dta3 <- dta3[dta3$bl_vac_yn %in% 1:2, ]

# Time from pcr to serology
dta3$time_from_pcr <-
  as.integer(as.Date(dta3$uc_ser_date_hr_kittrans) - dta3$pcr_date)

# Inconsistencies in dates
dta3[!is.na(dta3$time_from_pcr) & dta3$time_from_pcr < 0,
      c("hid", "pcr_date", "uc_ser_date_hr_kittrans")]

# Hospitalisation
dta3$hosp <-  with(dta3, as.integer(
  !is.na(bl_hosp_start) & !is.na(pcr_date) & 
    abs(bl_hosp_start - pcr_date) <= 14))

# Time from vaccination to serology
dta3$time_from_vac <-
  as.integer(as.Date(dta3$uc_ser_date_hr_kittrans) - dta3$bl_vac_first)
dta3[dta3$bl_vac_yn == 1 & is.na(dta3$time_from_vac), ]

# Inconsistencies in dates
dta3[!is.na(dta3$time_from_vac) &
        (dta3$time_from_vac < 0 | dta3$time_from_vac > 30),
      c("hid", "bl_vac_first", "uc_ser_date_hr_kittrans", "time_from_vac",
        "serol_igg")]

# --------------------------------------------------------------------------- #

# Figures
options(scipen=10000)
K <- list(with_hosp = 1, without_hosp = 2, vacc = 3, vacc_pos = 4)
fig_log_igg <- lapply(K, function(k) {
  dta_pcr <- subset(dta3, serol_igg == 1 & bl_vac_yn == 2 & pcr_result %in% 1)
  dta_serol_nvac <- subset(dta3, serol_igg == 1 & bl_vac_yn == 2)
  dta_serol_vac <- subset(dta3, serol_igg == 1 & bl_vac_yn == 1)
  dta_vac <- subset(dta3, bl_vac_yn == 1)
  if (k == 1) {
    dta_pcr$hosp <- factor(dta_pcr$hosp, 0:1, c("ok", "hospitalized"))
  } else if (k == 2) {
    dta_pcr <- subset(dta_pcr, hosp == 0)
    dta_serol_nvac <- subset(dta_serol_nvac, hosp == 0)
    dta_serol_vac <- subset(dta_serol_vac, hosp == 0)
  } else {
    dta_pcr <- subset(dta_pcr, hosp == 0)
    dta_serol_nvac <- subset(dta_serol_nvac, hosp == 0)
    dta_serol_vac <- subset(dta_serol_vac, hosp == 0)
    dta_vac <- subset(dta_vac, !is.na(time_from_vac))
    dta_vac$status <- with(dta_vac, droplevels(
      factor(serol_igg + (pcr_result %in% 1) * 2, 0:3,
             c("serol-", "serol+", "neg+pcr", "TDR/PCR+"))))
    if (k == 4) {
      dta_vac <- subset(dta_vac, serol_igg == 1)
      dta_vac$status <- droplevels(dta_vac$status)
    }
  }
  if (k %in% 1:2) {
    fig <- ggplot(dta_pcr, aes(x = time_from_pcr, y = uc_labo_coviggr_v2))
  } else {
    fig <- ggplot(dta_vac, aes(x = time_from_vac, y = uc_labo_coviggr_v2))
  } 
  fig <- fig +
    geom_smooth(method = lm, formula = y ~ x, se = TRUE, color = "grey40")
  if (k == 1) {
    fig <- fig +
      geom_point(aes(color = hosp)) +
      scale_color_manual(breaks = c("hospitalized"),
                         values = c("black", "red"))
  } else if (k == 2) {
    fig <- fig + geom_point()
  } else if (k == 3) {
    fig <- fig +
      geom_point(aes(color = status)) +
      scale_color_manual(values = c("blue", "black", "red"))
  } else {
    fig <- fig +
      geom_point(aes(color = status)) +
      scale_color_manual(values = c("black", "red"))
  }
  if (k %in% 1:2) {
    x.breaks <- c(0, 100, 200, 300)
    x.bp <- c(350, 400, 430)
    w.bp <- 15
    fig <- fig +
      geom_boxplot(aes(x = x.bp[1], y = uc_labo_coviggr_v2), width = w.bp) +
      geom_vline(xintercept = 375, linetype = "longdash") +
      annotate("text", x = 250, y = 8, label = paste("N =", nrow(dta_pcr)),
               size = 3) +
      labs(x = "Time between positive SARS-CoV-2 TDR/PCR and serology (days)")
  } else {
    x.breaks <- c(0, 10, 20)
    x.bp <- c(33, 30, 26)
    w.bp <- 2
    fig <- fig +
      geom_boxplot(aes(x = x.bp[1], y = uc_labo_coviggr_v2), width = w.bp,
                   data = dta_pcr) +
      geom_vline(xintercept = 28, linetype = "longdash") +
      labs(x = "Time between vaccination and serology (days)")
  }
  if (k == 3) {
    fig <- fig +
      annotate("text", x = 20, y = .1, label = paste("N =", nrow(dta_vac)),
               size = 3)
  } else if (k == 4) {
    fig <- fig +
      annotate("text", x = 20, y = 8, label = paste("N =", nrow(dta_vac)),
               size = 3)
  }
  if (k == 3) {
    fig <- fig +
      scale_y_continuous(breaks = 10^(-1:2), trans = "log10",
                         labels = c("0.1", "1", "10", "100"))
  } else {
    fig <- fig + scale_y_continuous(trans = "log10")
  }
  fig <- fig +
    geom_boxplot(aes(x = x.bp[2], y = uc_labo_coviggr_v2), width = w.bp,
                 data = dta_serol_nvac) +
    geom_boxplot(aes(x = x.bp[3], y = uc_labo_coviggr_v2), width = w.bp,
                 data = dta_serol_vac) +
    scale_x_continuous(
      breaks = c(x.breaks, x.bp),
      labels = c(x.breaks,
                 paste0("TDR/PCR+\nvacc-\n(", nrow(dta_pcr), ")"),
                 paste0("serol+\nvacc-\n(", nrow(dta_serol_nvac), ")"),
                 paste0("serol+\nvacc+\n(", nrow(dta_serol_vac), ")"))
    ) +
    theme_bw()
  if (k %in% c(1, 3:4)) {
    fig <- fig +
      theme(legend.title = element_blank(),
            legend.justification = 0:1,
            legend.position = c(0.02, 0.98),
            legend.key.size = unit(0, "points"),
            legend.text = element_text(size = rel(.6)),
            legend.box.background = element_rect(),
            legend.direction = "horizontal") +
      guides(color = guide_legend(override.aes = list(size = rel(.6))))
  }
  fig <- fig +
    labs(y = "IgG antibody level (log scale)")
  return(fig)
})
x11();fig_log_igg$with_hosp
x11();fig_log_igg$without_hosp
x11();fig_log_igg$vacc
x11();fig_log_igg$vacc_pos
jpeg("~/fig_log_igg_with_hosp.jpg", height = 3600, width = 7200, res = 1024)
print(fig_log_igg$with_hosp)
dev.off()
jpeg("~/fig_log_igg_without_hosp.jpg", height = 3600, width = 7200, res = 1024)
print(fig_log_igg$without_hosp)
dev.off()
jpeg("~/fig_log_igg_vacc.jpg", height = 4800, width = 7200, res = 1024)
print(fig_log_igg$vacc)
dev.off()
jpeg("~/fig_log_igg_vacc_pos.jpg", height = 4800, width = 7200, res = 1024)
print(fig_log_igg$vacc_pos)
dev.off()

# --------------------------------------------------------------------------- #

# Survival analysis
fit_igg <- survfit(Surv(time_from_pcr, 1 - serol_igg) ~ 1, data = dta3)
surv_igg <- ggsurvplot(fit_igg, palette = "#2E9FDF",
                       legend = "none", title = "IgG")

fit_iga <- survfit(Surv(time_from_pcr, 1 - serol_iga) ~ 1, data = dta3)
surv_iga <- ggsurvplot(fit_iga, palette = "#2E9FDF",
                       legend = "none", title = "IgA")

jpeg("~/surv_igg.jpg", height = 3600, width = 7200, res = 1024)
print(surv_igg)
dev.off()
jpeg("~/surv_iga.jpg", height = 3600, width = 7200, res = 1024)
print(surv_iga)
dev.off()
