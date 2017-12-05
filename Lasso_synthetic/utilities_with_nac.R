
require(biglasso)
require(plyr)
require(ggplot2)
require(mvtnorm)

which_case <- function(n, p, q, eff.nonzero, corr, sigma) {
  if (length(n) > 1) {
    sim.case <- "vary_n"
    which.vary <- n
    which.vary.c <- 'n'
    
  } else if (length(p) > 1) {
    sim.case <- "vary_p"
    which.vary <- p
    which.vary.c <- 'p'
    
  } else if (length(q) > 1) {
    sim.case <- 'vary_q'
    which.vary <- q
    which.vary.c <- 'q'
    
  } else if (length(eff.nonzero) > 1) {
    sim.case <- 'vary_beta'
    which.vary <- eff.nonzero
    which.vary.c <- 'eff.nonzero'
    
  } else if (length(corr) > 1) {
    sim.case <- 'vary_corr'
    which.vary <- corr
    which.vary.c <- 'corr'
    
  } else if (length(sigma) > 1) {
    sim.case <- 'vary_sigma'
    which.vary <- sigma
    which.vary.c <- 'sigma'
    
  } else {
    sim.case <- 'vary_NA'
    which.vary <- n
    which.vary.c <- 'NA'
  }
  list(sim.case = sim.case, which.vary = which.vary, which.vary.c = which.vary.c)
}

sim_with_nac <- function(n, p, q, eff.nonzero, corr, sigma, rep, methods, eps, 
                         lam.min, lam.log, seed, backingfile, descrpfile, backingpath) {
  
  case <- which_case(n = n, p = p, q = q, eff.nonzero = eff.nonzero, corr = corr, sigma = sigma)
  sim.case <- case$sim.case
  which.vary <- case$which.vary
  which.vary.c <- case$which.vary.c
 
  parms <- list(case = case, n = n, p = p, q = q, eff.nonzero = eff.nonzero,
                corr = corr, sigma = sigma, rep = rep, methods = methods, eps = eps,
                lam.min = lam.min, lam.log = lam.log)

  cat("\nStart simulation: ", format(Sys.time()))
  cat("\n============================================================\n")
  cat("\nR sessionInfo: \n\n")
  print(sessionInfo())
  
  ## print out simulation setting
  cat("\n Simulation case: ", sim.case)
  cat("\n============================================================\n")
  cat("\t Simulation settings: \n")
  cat("\t ---------------------------------------------\n")
  cat("\t n = ", n, '\n')
  cat("\t p = ", p, '\n')
  cat("\t q = ", q, '\n')
  cat("\t eff.nonzero = ", eff.nonzero, '\n')
  cat("\t corr = ", corr, '\n')
  cat("\t sigma = ", sigma, '\n')
  cat("\t rep = ", rep, '\n')
  cat("\t methods = ", methods, '\n')
  cat("\t eps = ", eps, '\n')
  cat("\t lam.min = ", lam.min, '\n')
  cat("\t lam.log = ", lam.log, '\n')
  cat("\n============================================================\n\n")
  
  
  time.all <- array(NA, dim = c(rep, length(methods), length(which.vary)))
  time.model <- array(NA, dim = c(rep, length(methods), length(which.vary)))
  
  for (i in 1:length(which.vary)) {
    cat("\t", which.vary.c, " = ", which.vary[i], "; start time: ", format(Sys.time()), '\n')  
    cat("\t---------------------------------------------\n")
    
    for (j in 1:rep) {
      
      new.seed <- seed + 100*i+j
      set.seed(new.seed)
      
      time.a <- NULL
      time.m <- NULL
      
      if (sim.case == 'vary_n') {
        beta.nonzero <- runif(q, -eff.nonzero, eff.nonzero)
        beta <- rep(0, p)
        nonzero.id <- sample(p, q)
        beta[nonzero.id] <- beta.nonzero
        
        if (corr < 1 && corr > 0) {
          Sigma <- matrix(corr, ncol = p, nrow = p)
          diag(Sigma) <- 1
          x <- rmvnorm(n[i], sigma = Sigma, method = 'chol')
        } else {
          x <- matrix(rnorm(p * n[i]), ncol = p)
        }
        y <- x %*% beta + sigma * rnorm(n[i])
      } else if (sim.case == 'vary_p') {
       
        beta.nonzero <- runif(q, -eff.nonzero, eff.nonzero)
        beta <- rep(0, p[i])
        nonzero.id <- sample(p[i], q)
        beta[nonzero.id] <- beta.nonzero
        
        if (corr < 1 && corr > 0) {
          Sigma <- matrix(corr, ncol = p[i], nrow = p[i])
          diag(Sigma) <- 1
          x <- rmvnorm(n, sigma = Sigma, method = 'chol')
        } else {
          x <- matrix(rnorm(p[i] * n), ncol = p[i])
        }
        y <- x %*% beta + sigma * rnorm(n)
      } else if (sim.case == 'vary_q') {
        
        beta.nonzero <- runif(q[i], -eff.nonzero, eff.nonzero)
        beta <- rep(0, p)
        nonzero.id <- sample(p, q[i])
        beta[nonzero.id] <- beta.nonzero
        
        if (corr < 1 && corr > 0) {
          Sigma <- matrix(corr, ncol = p, nrow = p)
          diag(Sigma) <- 1
          x <- rmvnorm(n, sigma = Sigma, method = 'chol')
        } else {
          x <- matrix(rnorm(p * n), ncol = p)
        }
        y <- x %*% beta + sigma * rnorm(n)
      } else if (sim.case == 'vary_beta') {
        
        beta.nonzero <- runif(q, -eff.nonzero[i], eff.nonzero[i])
        beta <- rep(0, p)
        nonzero.id <- sample(p, q)
        beta[nonzero.id] <- beta.nonzero
        
        if (corr < 1 && corr > 0) {
          Sigma <- matrix(corr, ncol = p, nrow = p)
          diag(Sigma) <- 1
          x <- rmvnorm(n, sigma = Sigma, method = 'chol')
        } else {
          x <- matrix(rnorm(p * n), ncol = p)
        }
        x <- matrix(rnorm(n * p), ncol = p)
        y <- x %*% beta + sigma * rnorm(n)
      } else if (sim.case == "vary_corr") {

        beta.nonzero <- runif(q, -eff.nonzero, eff.nonzero)
        beta <- rep(0, p)
        nonzero.id <- sample(p, q)
        beta[nonzero.id] <- beta.nonzero
        ## correlation matrix
        if (corr[i] < 1 && corr[i] > 0) {
          Sigma <- matrix(corr[i], ncol = p, nrow = p)
          diag(Sigma) <- 1
          x <- rmvnorm(n, sigma = Sigma, method = 'chol')
        } else {
          x <- matrix(rnorm(p * n), ncol = p)
        }
        y <- x %*% beta + sigma * rnorm(n)
        
      } else if (sim.case == 'vary_sigma') {
        beta.nonzero <- runif(q, -eff.nonzero, eff.nonzero)
        beta <- rep(0, p)
        nonzero.id <- sample(p, q)
        beta[nonzero.id] <- beta.nonzero
        
        if (corr < 1 && corr > 0) {
          Sigma <- matrix(corr, ncol = p, nrow = p)
          diag(Sigma) <- 1
          x <- rmvnorm(n, sigma = Sigma, method = 'chol')
        } else {
          x <- matrix(rnorm(p * n), ncol = p)
        }
        y <- x %*% beta + sigma[i] * rnorm(n)
      } else if (sim.case == 'vary_NA') {
        beta.nonzero <- runif(q, -eff.nonzero, eff.nonzero)
        beta <- rep(0, p)
        nonzero.id <- sample(p, q)
        beta[nonzero.id] <- beta.nonzero
        
        if (corr < 1 && corr > 0) {
          Sigma <- matrix(corr, ncol = p, nrow = p)
          diag(Sigma) <- 1
          x <- rmvnorm(n, sigma = Sigma, method = 'chol')
        } else {
          x <- matrix(rnorm(p * n), ncol = p)
        }
        y <- x %*% beta + sigma * rnorm(n)
      }

      X.bm <- as.big.matrix(x, backingfile = backingfile, descriptorfile = descrpfile,
                            backingpath = backingpath, type = 'double')
     
      fit.hsr0 <- biglasso(X.bm, y, family = 'gaussian', screen = 'SSR-BEDPP', 
                           lambda.log.scale = lam.log, lambda.min = lam.min,
                           ncores = 4, eps = eps)
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
     
      file.remove(paste0(backingpath, '/', backingfile))
      file.remove(paste0(backingpath, '/', descrpfile))
     
      time.all[j, , i] <- time.a
      time.model[j, , i] <- time.m
      
      cat("\t\trep", j, '; time.all = ', format(time.a, nsmall=3),
          "; time.model = ", format(time.m, nsmall=3), "\n")
      
    }
    cat("\t", which.vary.c, " = ", which.vary[i], "; end time: ", format(Sys.time()), '\n') 
    cat("\n============================================================\n")
    
  }

  cat("\nEnd simulation: ", format(Sys.time()), "\n")
  cat("\n============================================================\n")
  
  list(time.all = time.all, time.model = time.model, parms = parms)
  
}


plot_res <- function(time, filename, xlab, ylab = 'Computing time (s)', 
                     width, height, resolution, size, alpha) {
 
  gp <- ggplot(time, aes(x = Which.vary, y = time, color = Method)) +
    geom_line(size = size, alpha = alpha) +
    scale_x_continuous(breaks = pretty(range(time$Which.vary)),
                       limits = range(time$Which.vary)) +
    scale_y_continuous(breaks = pretty(range(time$time))) +
    xlab(xlab) +
    ylab(ylab) +
    theme_bw() + 
    theme(legend.position = c(.2, .7))
  
  png(filename = filename, width = width, height = height,
      units = 'in', res = resolution)
  print(gp)
  dev.off() 
  
  gp
  
}


post_analysis_with_nac <- function(res, date, xlab = NULL, ylab = NULL,
                                   width = 5, height = 4, resolution = 144, 
                                   size = 1, alpha = 1) {
  methods <- res$parms$methods
  which.vary <- res$parms$case$which.vary
  sim.case <- res$parms$case$sim.case

  cat("\n============================================================\n")
  cat("\n Time.all Mean: \n\n")
  time.mean.all <- apply(res$time.all, c(2, 3), mean, na.rm = TRUE)
  rownames(time.mean.all) <- methods
  colnames(time.mean.all) <- which.vary
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
  colnames(time.se.all) <- which.vary
  print(time.se.all)
  
  cat("\n--------------------------------------------\n")
  cat("\n Time.model Mean: \n\n")
  time.mean.model <- apply(res$time.model, c(2, 3), mean, na.rm = TRUE)
  rownames(time.mean.model) <- methods
  colnames(time.mean.model) <- which.vary
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
  colnames(time.se.model) <- which.vary
  print(time.se.model)
  
  ## plots
  # -----------------------------------------------------------------------------
  if (is.null(xlab)) {
    if (sim.case == "vary_n") {
      xlab <- 'Number of observations'
    } else if (sim.case == "vary_p") {
      xlab <- 'Number of features'
    } else if (sim.case == 'vary_q') {
      xlab <- 'Number of active features'
    } else if (sim.case == 'vary_beta') {
      xlab <- 'Magnitude of beta'
    } else if (sim.case == 'vary_corr') {
      xlab <- 'Magnitude of correlation'
    } else if (sim.case == 'vary_sigma') {
      xlab <- 'Sigma'
    } else {
      xlab <- 'NA'
    }
  } 
  
  if (is.null(ylab)) {
    ylab <- 'Computing time (s)'
  }
  
  rule.name <- methods

  ##  plot time.all
  time.df.all <- data.frame(time = matrix(t(time.mean.all), ncol = 1, byrow = T),
                            Method = rep(rule.name, each = length(which.vary)),
                            Which.vary = rep(which.vary, length(methods)))
  time.df.all$Method <- factor(time.df.all$Method, methods)
  filename <- paste0(date, '_', sim.case, '_all.png')
  gp.all <- plot_res(time.df.all, filename, xlab, ylab, width=width, height, resolution, size, alpha) 

 
  ##  plot time.model
  time.df.model <- data.frame(time = matrix(t(time.mean.model), ncol = 1, byrow = T),
                              Method = rep(rule.name, each = length(which.vary)),
                              Which.vary = rep(which.vary, length(methods)))
  time.df.model$Method <- factor(time.df.model$Method, methods)
  
  # remove case n = 100
  time.df.model <- time.df.model[time.df.model$Which.vary != 100, ]
  
  filename <- paste0(date, '_', sim.case, '_model.png')
  gp.model <- plot_res(time.df.model, filename, xlab, ylab, width=width, height, resolution, size, alpha) 
  
  # ## plot time.all with active 
  # time.df.all.wac <- subset(time.df.all, Method %in% c("None", "SSR", "SEDPP", "SSR-Dome", "SSR-BEDPP"))
  # time.df.all.wac$Method <- revalue(time.df.all.wac$Method, c("None"="AC"))
  # filename <- paste0(date, '_', sim.case, '_all_wac.png')
  # gp.all.wac <- plot_res(time.df.all.wac, filename, xlab, ylab, width, height, resolution, size, alpha) 
  # 
  # ## plot time.all without active 
  # time.df.all.nac <- subset(time.df.all, Method %in% c("NS-NAC", "SSR-NAC", "SEDPP-NAC", "SSR-Dome-NAC", "SSR-BEDPP-NAC"))
  # time.df.all.nac$Method <- revalue(time.df.all.nac$Method, 
  #                                   c("NS-NAC"="NS", "SSR-NAC"='SSR', 'SEDPP-NAC'='SEDPP', 
  #                                     "SSR-Dome-NAC"="SSR-Dome", 'SSR-BEDPP-NAC'='SSR-BEDPP'))
  # filename <- paste0(date, '_', sim.case, '_all_nac.png')
  # gp.all.nac <- plot_res(time.df.all.nac, filename, xlab, ylab, width, height, resolution, size, alpha) 
  # 
  # ## plot time.model with active 
  # time.df.model.wac <- subset(time.df.model, Method %in% c("None", "SSR", "SEDPP", "SSR-Dome", "SSR-BEDPP"))
  # time.df.model.wac$Method <- revalue(time.df.model.wac$Method, c("None"="AC"))
  # filename <- paste0(date, '_', sim.case, '_model_wac.png')
  # gp.model.wac <- plot_res(time.df.model.wac, filename, xlab, ylab, width, height, resolution, size, alpha) 
  
  ## plot time.model without active 
  time.df.model.nac <- subset(time.df.model, Method %in% c("NS-NAC", "None", "SSR-NAC", "SEDPP-NAC", "SSR-Dome-NAC", "SSR-BEDPP-NAC"))
  time.df.model.nac$Method <- revalue(time.df.model.nac$Method, 
                                      c("NS-NAC"="NS", "None"="AC", "SSR-NAC"='SSR', 
                                        'SEDPP-NAC'='SEDPP', "SSR-Dome-NAC"="SSR-Dome", 
                                        'SSR-BEDPP-NAC'='SSR-BEDPP'))
  filename <- paste0(date, '_', sim.case, '_model_nac.png')
  gp.model.nac <- plot_res(time.df.model.nac, filename, xlab, ylab, width, height, resolution, size, alpha) 
  
  save(list = ls(), file = paste0(date, '_', sim.case, '_results.RData'))
  
}

