rm(list = ls())

require(biglasso)

## benchmark timings
sim_real <- function(n, p, rep, methods, eps, lam.min, lambda.log, seed, sample.y = FALSE) {
  
  parms <- list(n = n, p = p, rep = rep, methods = methods, eps = eps,
                lam.min = lam.min, lam.log = lambda.log, seed = seed)
  
  cat("\nStart simulation: ", format(Sys.time()))
  cat("\n============================================================\n")
  cat("\nR sessionInfo: \n\n")
  print(sessionInfo())
  
  ## print out simulation setting
  cat("\n============================================================\n")
  cat("\t Simulation settings: \n")
  cat("\t ---------------------------------------------\n")
  cat("\t n = ", n, '\n')
  cat("\t p = ", p, '\n')
  cat("\t rep = ", rep, '\n')
  cat("\t methods = ", methods, '\n')
  cat("\t eps = ", eps, '\n')
  cat("\t lam.min = ", lam.min, '\n')
  cat("\t lam.log = ", lambda.log, '\n')
  cat("\n============================================================\n\n")
  
  time.all <- array(NA, dim = c(rep, length(methods), length(p)))
  time.model <- array(NA, dim = c(rep, length(methods), length(p)))
 
  for (i in 1:length(p)) {
    cat("\tp =", p[i], "; start time: ", format(Sys.time()), '\n')
    cat("\t---------------------------------------------\n")
    
    for (j in 1:rep) {
      new.seed <- seed + 100*i+j
      set.seed(new.seed)
      
      time.a <- NULL
      time.m <- NULL
      
      if (sample.y) {
        ## sample response verctor for different replications (NYT and MNIST data only)
        y.idx <- sample(ncol(Y), 1)
        y <- as.numeric(Y[, y.idx])
      }
      
      fit.hsr0 <- biglasso(X.bm, y, family = 'gaussian', screen = 'SSR-BEDPP', 
                           lambda.log.scale = lambda.log, lambda.min = lam.min,
                           eps = eps)
      lambda <- fit.hsr0$lambda
      
      for (k in 1:length(methods)) {
        t <- system.time(
          fit <- biglasso(X.bm, y, penalty="lasso", family = "gaussian", 
                          screen = methods[k], lambda = lambda, eps = eps)
        )
        time.m <- c(time.m, fit$time)
        time.a <- c(time.a, as.numeric(t['elapsed']))
        
        rm(fit)
        gc()
        
      }
      
      time.all[j, , i] <- time.a
      time.model[j, , i] <- time.m
      
      cat("\t\trep", j, '; time.all = ', format(time.a, nsmall=3),
          "; time.model = ", format(time.m, nsmall=3), "\n")
      
      
    }
    cat("\tp =", p[i], "; end time: ", format(Sys.time()), "\n")  
    cat("\n============================================================\n")
    
  }
  
  list(time.all = time.all,
       time.model = time.model,
       parms = parms)

}

summary_res <- function(res, dataname) {
  cat("\n============================================================\n")
  cat("\n Time.all: \n\n")
  print(res$time.all)
  
  cat("\n Time.all Mean: \n\n")
  time.mean.all <- apply(res$time.all, c(2, 3), mean, na.rm = TRUE)
  rownames(time.mean.all) <- methods
  colnames(time.mean.all) <- dataname
  print(time.mean.all)
  
  cat("\n Time.all SE: \n\n")
  time.se.all <- apply(res$time.all, c(2, 3), function(x) {
    x <- x[!is.na(x)]
    if (length(x) <= 1) {
      return(NA)
    } else {
      return(sd(x) / sqrt(length(x)))
    }
  })
  rownames(time.se.all) <- methods
  colnames(time.se.all) <- dataname
  print(time.se.all)
  
  cat("\n--------------------------------------------\n")
  cat("\n Time.model: \n\n")
  print(res$time.model)
  
  cat("\n Time.model Mean: \n\n")
  time.mean.model <- apply(res$time.model, c(2, 3), mean, na.rm = TRUE)
  rownames(time.mean.model) <- methods
  colnames(time.mean.model) <- dataname
  print(time.mean.model)
  
  cat("\n Time.model SE: \n\n")
  time.se.model <- apply(res$time.model, c(2, 3), function(x) {
    x <- x[!is.na(x)]
    if (length(x) <= 1) {
      return(NA)
    } else {
      return(sd(x) / sqrt(length(x)))
    }
  })
  rownames(time.se.model) <- methods
  colnames(time.se.model) <- dataname
  print(time.se.model)
  
}

setwd("~/GitHub/biglasso_experiment/image_data/MNIST/20170203_res_kdd/")
load('~/GitHub/biglasso_experiment/image_data/MNIST/mnist_data.RData')
n <- nrow(X)
p <- ncol(X)
X.bm <- as.big.matrix(X, type = 'double')

# --------------------------------------------------------------
## benchmark params
# --------------------------------------------------------------
lambda.log <- F
eps <- 1e-8
lam.min <- 0.1
date <- Sys.Date()
seed <- 1234
rep <- 20
dataname <- "MNIST"
methods <- c("NS-NAC", "None", "SSR", "SSR-NAC", "SEDPP", "SEDPP-NAC",
             "SSR-Dome",  "SSR-Dome-NAC", "SSR-BEDPP",  "SSR-BEDPP-NAC")

res <- sim_real(n, p, rep, methods, eps, lam.min, lambda.log, seed, sample.y = TRUE)

summary_res(res, dataname)

save(res, file = paste0(date, "_", dataname, "_", eps, "_Results.RData"))
