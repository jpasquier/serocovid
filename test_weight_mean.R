
source("R/estimation_20200616.R", echo = TRUE)
data <- merge(data, within(strata, {wgt = N / n})[c("stratum", "wgt")],
              by = "stratum", sort = FALSE)
x <- data[data$genre1 == 1, "serol"]
w <- data[data$genre1 == 1, "wgt"]
weighted.mean(x, w)
weighted.ttest.ci <- function(x, weights, conf.level = 0.95) {
    require(Hmisc)
    nx <- length(x)
    df <- nx - 1
    vx <- wtd.var(x, weights, normwt = TRUE) ## From Hmisc
    mx <- weighted.mean(x, weights)
    stderr <- sqrt(vx/nx)
    tstat <- mx/stderr ## not mx - mu
    alpha <- 1 - conf.level
    cint <- qt(1 - alpha/2, df)
    cint <- tstat + c(-cint, cint)
    cint * stderr
}
weighted.ttest.ci(x, w)
wtd.var(x, w, normwt = TRUE) - 
  1 / sum(w) * sum(w * (x - sum(x * w) / sum(w))^2)
