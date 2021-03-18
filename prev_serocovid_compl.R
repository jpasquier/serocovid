library(writexl)

# Working directory
setwd("~/Projects/SerocoViD")

#
load("results/prev_serocovid_20210218.dta")

#
ppos <- do.call(rbind, lapply(1:3, function(k) {
  dta <- get(paste0("dta", k))
  dta <- subset(dta, stratum >= 4)
  dta$all <- "all"
  dta$agegroup3 <- factor(dta$stratum, 4:8,
                          c("15-19", "20-64", "20-64", ">=65", ">=65"))
  pop <- list(pop1, pop1, pop3)[[k]]
  ppos <- Mean(data = dta, variable = "serol_any", stratum = "stratum",
          domain = c("all", "agegroup3"), pop = pop)
  ppos$v <- NULL
  names(ppos)[names(ppos) == "y"] <- "ppos"
  ppos <- cbind(phase = k, antibody = "IgG or IgA", ppos)
}))

#
write_xlsx(ppos, "results/prev_serocovid_20210218_compl_20210302.xlsx")
