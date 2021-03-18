library(writexl)

# Working directory
setwd("~/Projects/SerocoViD")

#
load("results/prev_serocovid_20210218.dta")

#
dta3$all <- "all"
dta3$agegroup3 <- factor(dta3$stratum, 4:8,
                         c("15-19", "20-64", "20-64", ">=65", ">=65"))
ppos <- Mean(data = dta3, variable = "serol_any", stratum = "stratum",
        domain = c("all", "agegroup3"), pop = pop1)
ppos$v <- NULL
names(ppos)[names(ppos) == "y"] <- "ppos"
ppos <- cbind(antibody = "IgG or IgA", ppos)

#
write_xlsx(ppos, "results/prev_serocovid_20210218_compl_20210301.xlsx")
