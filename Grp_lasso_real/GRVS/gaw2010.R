rm(list = ls())
gc()

require(grpreg)
require(ggplot2)
require(plyr)
sessionInfo()

load("gaw2010.RData")

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

set.seed(1234)
y.ind <- sample(ncol(Y), rep)

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
  y <- Y[, y.ind[i]]
  for (j in 1:length(methods)) {
    t <- system.time(
      fit <- grpreg(X, y, group = group, penalty="grLasso", screen = methods[j],
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

save(list = ls(), file = paste0(date, "_gaw2010_Y_results.RData"))

