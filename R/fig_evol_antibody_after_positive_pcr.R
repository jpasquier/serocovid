library(RCurl)
library(dplyr)
library(ggplot2)
library(ggh4x)
library(parallel)

# Working directory
setwd("~/Projects/SerocoViD")

# ------------------------------- REDCap data ------------------------------- #

# API
uri   <- "https://redcap.unisante.ch/api/"

# Tokens
tokens <- c("corona_immunitas", "personal_data", "research_data",
            "corona_immunitas_vague_6")
tokens <- lapply(setNames(tokens, tokens), function(z) {
  z <- paste0("misc/redcap_", z, ".token")
  readChar(z, file.info(z)$size - 1)
})

# Import
tmp_file <- "/mnt/ramdisk/serocovid_data.rda"
if (file.exists(tmp_file)) {
  load(tmp_file)
} else {
  serocovid_data <- mclapply(tokens, function(token) {
    read.csv(text = postForm(
      uri = uri,
      token = token,
      content = 'record',
      format = 'csv'
    ))
  })
  save(serocovid_data, file = tmp_file)
}
rm(uri, tokens, tmp_file)

# -------------------------- Prepare research data -------------------------- #

dta <- bind_rows(
  select(
    serocovid_data$corona_immunitas, hid, uc_labo_qc_v4,
    matches('^uc_labo_covigg(l|r)(_zscore)?_v(2|4)'), uc_ser_date_hr_kittrans,
    matches('^pcr_(prior|date|result)(_multiple)?$'),
    matches('^bl_vac_(yn|first|second|third)(_ph4)?$'),
  ),
  select(
    serocovid_data$corona_immunitas_vague_6, hid, 
    matches('^(spike|nuc)_igg_(qual|ratio)$'), nuc_igg_qual,
    uc_ser_date_hr_ps_v6, bl_test_pos_yn, starts_with('bl_infection_date'),
    matches('^bl_vac_(yn|first|second|third|fourth)(_date)?$')
  )
) %>%
  left_join(
    serocovid_data$personal_data %>%
      unique() %>%
      select(hid = uc_s_participant_hid, type = uc_s_type_participant,
        matches('^uc_status_v4_(check_inf___(1|2)|time_prel_(centre|home))$')),
    by = "hid"
  ) %>%
  filter(
    type == 13 & !grepl("^123123", hid) |
    type == 14 & !(hid %in% c("13131313", "987654321")) |
    type == 15 & !(hid %in%
      c("6-841282", "6-76BKAW", "6-BHSMPF", paste0("test", 1:6)))
  ) %>%
  mutate(
    serol = case_when(
      type == 13 ~ case_when(
        grepl("(Négatif|Ind(e|é)terminé)", uc_labo_coviggl_v2) ~ 0L,
        grepl("Positif", uc_labo_coviggl_v2) ~ 1L,
        TRUE ~ NA_integer_
      ),
      type == 14 ~ case_when(
                     uc_labo_qc_v4 != "OK" ~ NA_integer_,
                     uc_labo_coviggl_v4 == "négatif" ~ 0L,
                     uc_labo_coviggl_v4 == "positif" ~ 1L,
                   ),
      type == 15 ~ as.integer(spike_igg_qual == 1 | nuc_igg_qual == 1)
    ),
    serol_score = case_when(
      type == 13 ~ uc_labo_coviggr_v2,
      type == 14 ~ uc_labo_coviggr_zscore_v4,
      type == 15 ~ pmax(spike_igg_ratio, nuc_igg_ratio)
    ),
    serol_date = as.Date(case_when(
      type == 13 ~ uc_ser_date_hr_kittrans,
      type == 14 ~ if_else(
        uc_status_v4_check_inf___1 == 1 | uc_status_v4_check_inf___2 == 1,
        uc_status_v4_time_prel_centre, uc_status_v4_time_prel_home),
      type == 15 ~ uc_ser_date_hr_ps_v6
    )),
    positive_prior_test = case_when(
      type %in% 13:14 ~ case_when(
        pcr_prior == 0 ~ 0L,
        pcr_prior == 1 ~ case_when(
          as.Date(pcr_date) > serol_date ~ 0L,
          pcr_result == 1 ~ 1L,
          pcr_result == 2 ~ 0L,
          TRUE ~ NA_integer_
        ),
        pcr_prior == 3 ~ case_when(
          as.Date(pcr_date_multiple) > serol_date ~ 0L,
          pcr_result_multiple == 0 ~ 0L,
          pcr_result_multiple == 1 ~ 1L,
          TRUE ~ NA_integer_
        ),
        TRUE ~ NA_integer_
      ),
      type == 15 ~ as.integer(bl_test_pos_yn)
    ),
    positive_test_date = case_when(
      type %in% 13:14 ~ as.Date(case_when(
        positive_prior_test == 1 ~ case_when(
          pcr_prior == 1 ~ pcr_date,
          pcr_prior == 3 ~ pcr_date_multiple
        ),
        TRUE ~ NA_character_
      )),
      type == 15 ~ do.call(pmax, c(across(starts_with('bl_infection_date'), 
                                          as.Date), na.rm = TRUE))
    ),
    time_from_positive_test = as.integer(serol_date - positive_test_date),
    vac = case_when(
      type == 13 ~ recode(bl_vac_yn, `2` = 0L, `3` = NA_integer_),
      type == 14 ~ recode(bl_vac_yn_ph4, `2` = 0L, `3` = NA_integer_),
      type == 15 ~ recode(bl_vac_yn, `2` = NA_integer_),
    ),
    n_vac_before_serol = {
      z <- .[grep('^bl_vac_(first|second|third|fourth)(_(ph4|date))?$',
                  names(.))]
      z <- do.call(cbind, lapply(z, function(w) as.Date(w) <= serol_date))
      ifelse(apply(is.na(z), 1, all), NA, apply(z, 1, sum, na.rm = TRUE))
    },
    vac = if_else(vac %in% 1 & n_vac_before_serol %in% 0, 0L, vac),
    vac_date = {
      z <- .[grep('^bl_vac_(first|second|third|fourth)(_(ph4|date))?$',
                  names(.))]
      z <- lapply(z, function(w) {
        w <- as.Date(w)
        if_else(w <= serol_date, w, as.Date(NA))
      })
      if_else(vac == 1L, do.call(pmax, c(z, na.rm = TRUE)), as.Date(NA)) 
    },
    time_from_vac = as.integer(serol_date - vac_date),
    visit = type - 10
  )

# --------------------------- Inconsistent dates ---------------------------- #

with(dta, table(visit, time_from_positive_test < 0 |
                         time_from_positive_test > 1095))
with(dta, table(visit, time_from_vac < 0 | time_from_vac > 1095))
dta %>%
  filter(visit == 4, time_from_vac < 0 | time_from_vac > 1095) %>%
  select(hid, visit, matches('^uc_status_v4_time_prel_(centre|home)$'),
         serol_date, matches('^bl_vac_(first|second|third)_ph4$'),
         time_from_vac)

# --------------------------------- Figures --------------------------------- #

fig_vac <- lapply(3:5 , function(v) {
  dta %>%
    filter(visit == v, serol == 1, !is.na(time_from_vac), !is.na(serol),
           !is.na(serol_score), time_from_vac >= 0 & time_from_vac <= 1095) %>%
    select(serol_score, time_from_vac, n_vac_before_serol) %>%
    mutate(doses = factor(n_vac_before_serol)) %>%
    ggplot(aes(x = time_from_vac, y = serol_score)) +
    geom_point(aes(colour = doses)) +
    geom_smooth(method = lm, formula = y ~ x, se = TRUE) +
    scale_y_continuous(trans = "log10") +
    theme_bw() +
    labs(x = "Time between vaccination (last dose) and serology (days)",
         y = "IgG antibody level")
})

fig_test <- lapply(1:2, function(k) {
  fig <- dta %>%
    filter(visit %in% list(c(3, 5), 4)[[k]], vac == 0L, !is.na(serol),
           !is.na(serol_score), !is.na(time_from_positive_test)) %>%
    select(visit, serol, serol_score, time_from_positive_test) %>%
    mutate(visit = droplevels(factor(visit, 3:5, paste('Survey', 3:5)))) %>%
    ggplot(aes(x = time_from_positive_test, y = serol_score)) +
    geom_point() +
    geom_smooth(method = lm, formula = y ~ x, se = TRUE) +
    scale_y_continuous(trans = "log10") +
    labs(x = "Time between positive SARS-CoV-2 TDR/PCR and serology (days)",
         y = "IgG antibody level (log scale)")
  if (k == 1) {
    fig <- fig +
      facet_grid(cols = vars(visit), scales = 'free_x', space = 'free_x')
  }
  return(fig)
})

d <- format(Sys.Date(), "%Y%m%d")
for (k in 1:3) {
  jpeg(paste0("results/fig_evol_igg_after_vac_", k + 2, "_", d, ".jpg"),
       height = 3600, width = 7200, res = 512)
  print(fig_vac[[k]])
  dev.off()
}
for (k in 1:2) {
  jpeg(paste0("results/fig_evol_igg_after_positive_test_",
              c("3_5", "4")[k], "_", d, ".jpg"),
       height = 3600, width = 7200, res = 768)
  print(fig_test[[k]])
  dev.off()
}
