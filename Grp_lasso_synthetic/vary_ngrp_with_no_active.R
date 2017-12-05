
rm(list = ls())
gc()

source("~/GitHub/grpreg_experiment/simu/vary_ngrp_compare_no_active/0413_ngrp_true_10_1e-5/utilties_with_no_active.R")
dir <- "~/GitHub/grpreg_experiment/simu/vary_ngrp_compare_no_active/0413_ngrp_true_10_1e-5/"
setwd(dir)
seed <- 1234
# set.seed(seed)
date <- Sys.Date()
# date <- "2017-02-05"

n <- 1000
ngrp <- c(100, 200, 500, 1000, 2000, 5000, 10000)
grp.size <- 10
ngrp.true <- 10
eff.nonzero <- 1
sigma <- 0.1
rep <- 20

methods <- c("No-Active", "None", "SSR-No-Active", "SEDPP-No-Active", "SSR-BEDPP-No-Active")
methods.name <- c("NS-NAC", "NS-AC", "SSR-NAC", "SEDPP-NAC", "SSR-BEDPP-NAC")

eps <- 1e-5
lam.min <- 0.1
lam.log <- FALSE

# ------------------------------------------
# testing
# n <- 100
# ngrp <- c(100, 200, 500)
# rep <- 2
# -----------------------------------------

res <- sim(n, grp.size, ngrp, ngrp.true, eff.nonzero, sigma, rep, methods, 
           methods.name, eps, lam.min, lam.log, seed)

post_analysis(res, date)
