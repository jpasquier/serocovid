library(writexl)

# Working directory
setwd("~/Projects/SerocoViD")

# Raw samples
smpl <- c(smpl1 = 1, smpl3 = 2, smpl4 = 4)
smpl <- lapply(smpl, function(s) {
  filename <- paste0("data-fso/COVID19_VD_V", s, "_Total.csv")
  subset(read.csv(filename, sep = ";"), !is.na(strate))
})

# Raw sample sizes
rss <- do.call(rbind, lapply(names(smpl), function(s) {
  d <- smpl[[s]]
  k <- substr(s, nchar(s), nchar(s))
  n <- setNames(as.data.frame(table(d$strate)), c("stratum", "n_raw_sample"))
  cbind(visit = k, n)
}))
write.table(rss, "results/raw_sample_sizes.csv", sep = ";", col.names = FALSE)
write_xlsx(rss, "results/raw_sample_sizes.xlsx")
