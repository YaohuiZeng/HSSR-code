rm(list = ls())
gc()

require(grpreg)
require(ggplot2)
require(plyr)
require(splines)

sessionInfo()

load("bcTCGA.RData")

# methods <- c("No-Active", "None", "SSR-No-Active", "SSR", "SEDPP-No-Active",
#              "SEDPP", "SSR-BEDPP-No-Active", "SSR-BEDPP")
# methods.name <- c("NS-NAC", "NS-AC", "SSR-NAC", "SSR-AC", "SEDPP-NAC", 
#                   "SEDPP-AC", "SSR-BEDPP-NAC", "SSR-BEDPP-AC")

methods <- c("No-Active", "None", "SSR-No-Active", "SEDPP-No-Active", "SSR-BEDPP-No-Active")
methods.name <- c("NS-NAC", "NS-AC", "SSR-NAC", "SEDPP-NAC", "SSR-BEDPP-NAC")

eps <- 1e-5
lambda.min <- 0.1
lam.log <- FALSE
rep <- 20
date <- Sys.Date()
# date <- "2017-04-06"

set.seed(1234)
p <- ncol(X)
n <- nrow(X)

## set up splines
# Restrict to n w/ highest variance (not really necessary)
# n <- 5000
# ind <- order(apply(X, 2, var), decreasing=TRUE)[1:n]
# names(ind) <- colnames(X)[ind]

# Create matrix of basis functions (B) and group assignments
df <- 5   # df per feature in spline basis
BB <- sapply(as.data.frame(X), ns, simplify = "array", df = df)
B <- t(apply(BB, 1, rbind))
# group <- rep(1:p, rep(df, p))
group <- rep(1:p, each = df)

# ------------------------------------------
# testing
# n <- 100
# ngrp <- c(100, 200, 500)
# rep <- 2
# -----------------------------------------

time.model <- matrix(NA, ncol = length(methods), nrow = rep)
time.all <- matrix(NA, ncol = length(methods), nrow = rep)
colnames(time.model) <- methods.name
colnames(time.all) <- methods.name

cat("Start time: ", format(Sys.time()), "\n\n")
for (i in 1:rep) {
  for (j in 1:length(methods)) {
    t <- system.time(
      fit <- grpreg(B, y, group = group, penalty="grLasso", screen = methods[j],
                    lambda.min = lambda.min, eps = eps, log.lambda = lam.log)
    )
    time.model[i, j] <- fit$time
  }
}
cat("End time: ", format(Sys.time()), "\n\n")

cat("Time.model: \n")
print(time.model)
cat("Mean of time.model: \n")
print(colMeans(time.model))
cat("SE of time.model: \n")
print(
  apply(time.model, 2, function(x) {
    sd(x) / sqrt(length(x))
  })
)

save(list = ls(), file = paste0(date, "_breast_spline_results.RData"))

