
# ============================================================
# Case 1: vary n
# ============================================================
rm(list = ls())
gc()

## Note: to replicate, change dir accordingly.
source("~/GitHub/biglasso_experiment/linear_rules/utilities_with_nac.R")
dir <- "~/GitHub/biglasso_experiment/linear_rules/vary_n/20170117_result_eps_1e-8/"
setwd(dir)
seed <- 1234
# date <- Sys.Date()
date <- "2017-01-17"

n <- c(100, 200, 500, 1000, 2000, 5000, 10000)
p <- 10000
q <- 20
eff.nonzero <- 1
corr <- 0
rep <- 20

methods <- c("NS-NAC", "None", "SSR", "SSR-NAC", "SEDPP", "SEDPP-NAC",
             "SSR-Dome",  "SSR-Dome-NAC", "SSR-BEDPP",  "SSR-BEDPP-NAC")
eps <- 1e-8
lam.min <- 0.1
lam.log <- FALSE
sigma <- 0.1
backingfile <- 'back.bin'
descrpfile <- 'descrb.desc'
backingpath <- getwd()

# ===================================
# testing parameters
# p <- 1000
# n <- c(100, 200, 500)
# rep <- 2
# ===================================

res <- sim_with_nac(n = n, p = p, q = q, eff.nonzero = eff.nonzero, corr = corr, 
           rep = rep, methods = methods, eps = eps, lam.min = lam.min, 
           lam.log = lam.log, sigma = sigma, seed = seed,
           backingfile = backingfile, descrpfile = descrpfile, 
           backingpath = backingpath)

post_analysis_with_nac(res, date)

