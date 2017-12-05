
require(grpreg)
require(ggplot2)
require(plyr)

which_case <- function(n, grp.size, ngrp, ngrp.true, eff.nonzero, sigma) {
  if (length(n) > 1) {
    sim.case <- "vary_n"
    which.vary <- n
    which.vary.c <- 'n'
    
  } else if (length(grp.size) > 1) {
    sim.case <- "vary_grp_size"
    which.vary <- grp.size
    which.vary.c <- 'grp.size'
    
  } else if (length(ngrp) > 1) {
    sim.case <- 'vary_ngrp'
    which.vary <- ngrp
    which.vary.c <- 'ngrp'
    
  } else if (length(ngrp.true) > 1) {
    sim.case <- 'vary_ngrp_true'
    which.vary <- ngrp.true
    which.vary.c <- 'ngrp.true'
    
  } else if (length(eff.nonzero) > 1) {
    sim.case <- 'vary_beta'
    which.vary <- eff.nonzero
    which.vary.c <- 'eff.nonzero'

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


sim <- function(n, grp.size, ngrp, ngrp.true, eff.nonzero, sigma, rep, 
                methods, methods.name, eps, lam.min, lam.log, seed) {
  
  case <- which_case(n, grp.size, ngrp, ngrp.true, eff.nonzero, sigma)
  sim.case <- case$sim.case
  which.vary <- case$which.vary
  which.vary.c <- case$which.vary.c
  
  parms <- list(case = case, n = n, grp.size = grp.size, ngrp = ngrp, 
                ngrp.true = ngrp.true, eff.nonzero = eff.nonzero, sigma = sigma, 
                rep = rep, methods = methods, methods.name = methods.name, eps = eps,
                lam.min = lam.min, lam.log = lam.log, seed = seed)

  cat("\n============================================================\n")
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
  cat("\t grp.size = ", grp.size, '\n')
  cat("\t ngrp = ", ngrp, '\n')
  cat("\t ngrp.true = ", ngrp.true, '\n')
  cat("\t eff.nonzero = ", eff.nonzero, '\n')
  cat("\t sigma = ", sigma, '\n')
  cat("\t rep = ", rep, '\n')
  cat("\t methods = ", methods, '\n')
  cat("\t methods.name = ", methods.name, '\n')
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
      time.a <- NULL
      time.m <- NULL
      
      new.seed <- seed + 100*i+j
      set.seed(new.seed)
      
      if (sim.case == 'vary_n') {
        p <- ngrp * grp.size
        beta.nonzero <- runif(ngrp.true, -eff.nonzero, eff.nonzero)
        beta <- rep(0, p)
        nonzero.id <- sample(ngrp, ngrp.true) # random true group indices
        for (k in 1:ngrp.true) {
          beta[((nonzero.id[k]-1)*grp.size+1):(nonzero.id[k]*grp.size)] <- rep(beta.nonzero[k], grp.size)
        }
        x <- matrix(rnorm(p * n[i]), ncol = p)
        y <- x %*% beta + sigma * rnorm(n[i])
        group <- rep(1:ngrp, each = grp.size)
        
      } else if (sim.case == 'vary_grp_size') {
        
        p <- ngrp * grp.size[i]
        beta.nonzero <- runif(ngrp.true, -eff.nonzero, eff.nonzero)
        beta <- rep(0, p)
        nonzero.id <- sample(ngrp, ngrp.true) # random true group indices
        for (k in 1:ngrp.true) {
          beta[((nonzero.id[k]-1)*grp.size[i]+1):(nonzero.id[k]*grp.size[i])] <- rep(beta.nonzero[k], grp.size[i])
        }
        x <- matrix(rnorm(p * n), ncol = p)
        y <- x %*% beta + sigma * rnorm(n)
        group <- rep(1:ngrp, each = grp.size[i])
        
      } else if (sim.case == 'vary_ngrp') {
        p <- ngrp[i] * grp.size
        beta.nonzero <- runif(ngrp.true, -eff.nonzero, eff.nonzero)
        beta <- rep(0, p)
        nonzero.id <- sample(ngrp[i], ngrp.true) # random true group indices
        for (k in 1:ngrp.true) {
          beta[((nonzero.id[k]-1)*grp.size+1):(nonzero.id[k]*grp.size)] <- rep(beta.nonzero[k], grp.size)
        }
        x <- matrix(rnorm(p * n), ncol = p)
        y <- x %*% beta + sigma * rnorm(n)
        group <- rep(1:ngrp[i], each = grp.size)
        
      } else if (sim.case == 'vary_ngrp_true') {
        p <- ngrp * grp.size
        beta.nonzero <- runif(ngrp.true[i], -eff.nonzero, eff.nonzero)
        beta <- rep(0, p)
        nonzero.id <- sample(ngrp, ngrp.true[i]) # random true group indices
        for (k in 1:ngrp.true[i]) {
          beta[((nonzero.id[k]-1)*grp.size+1):(nonzero.id[k]*grp.size)] <- rep(beta.nonzero[k], grp.size)
        }
        x <- matrix(rnorm(p * n), ncol = p)
        y <- x %*% beta + sigma * rnorm(n)
        group <- rep(1:ngrp, each = grp.size)
      
      } else if (sim.case == 'vary_beta') {
        p <- ngrp * grp.size
        beta.nonzero <- runif(ngrp.true, -eff.nonzero[i], eff.nonzero[i])
        beta <- rep(0, p)
        nonzero.id <- sample(ngrp, ngrp.true) # random true group indices
        for (k in 1:ngrp.true) {
          beta[((nonzero.id[k]-1)*grp.size+1):(nonzero.id[k]*grp.size)] <- rep(beta.nonzero[k], grp.size)
        }
        x <- matrix(rnorm(n * p), ncol = p)
        y <- x %*% beta + sigma * rnorm(n)
        group <- rep(1:ngrp, each = grp.size)
        
      } else if (sim.case == 'vary_sigma') {
        p <- ngrp * grp.size
        beta.nonzero <- runif(ngrp.true, -eff.nonzero, eff.nonzero)
        beta <- rep(0, p)
        nonzero.id <- sample(ngrp, ngrp.true) # random true group indices
        for (k in 1:ngrp.true) {
          beta[((nonzero.id[k]-1)*grp.size+1):(nonzero.id[k]*grp.size)] <- rep(beta.nonzero[k], grp.size)
        }
        x <- matrix(rnorm(n * p), ncol = p)
        y <- x %*% beta + sigma[i] * rnorm(n)
        group <- rep(1:ngrp, each = grp.size)

      } else if (sim.case == 'vary_NA') {
        p <- ngrp * grp.size
        beta.nonzero <- runif(ngrp.true, -eff.nonzero, eff.nonzero)
        beta <- rep(0, p)
        nonzero.id <- sample(ngrp, ngrp.true) # random true group indices
        for (k in 1:ngrp.true) {
          beta[((nonzero.id[k]-1)*grp.size+1):(nonzero.id[k]*grp.size)] <- rep(beta.nonzero[k], grp.size)
        }
        x <- matrix(rnorm(p * n), ncol = p)
        y <- x %*% beta + sigma * rnorm(n)
        group <- rep(1:ngrp, each = grp.size)
        
      }

      for (l in 1:length(methods)) {
        t <- system.time(
          fit <- grpreg(x, y, group = group, penalty="grLasso", screen = methods[l],
                        lambda.min = lam.min, eps = eps, log.lambda = FALSE)
        )
        time.m <- c(time.m, fit$time)
        time.a <- c(time.a, as.numeric(t['elapsed']))
        rm(fit)
        gc()
      }

      rm(x)
      rm(beta)
      gc()
      
      time.all[j, , i] <- time.a
      time.model[j, , i] <- time.m
      
      cat("\t\trep", j, '; time.all = ', format(time.a, nsmall=3),
          "; time.model = ", format(time.m, nsmall=3), "\n")
      # cat("\t\trep", j, '; time.model = ', time.m, "\n")

    }
    cat("\t", which.vary.c, " = ", which.vary[i], "; end time: ", format(Sys.time()), '\n') 
    cat("\n============================================================\n")
    
  }
 
  cat("\nEnd simulation: ", format(Sys.time()), "\n")
  cat("\n============================================================\n")
  
  list(time.all = time.all,
       time.model = time.model,
       parms = parms
  )

}

post_analysis <- function(res, date, width = 6, height = 4, resolution = 120) {
  methods <- res$parms$methods.name
  # methods <- res$parms$methods
  which.vary <- res$parms$case$which.vary
  sim.case <- res$parms$case$sim.case
  
  if (sim.case == "vary_n") {
    xlab <- 'Number of observations'
  } else if (sim.case == "vary_grp_size") {
    xlab <- 'Group size'
  } else if (sim.case == 'vary_ngrp') {
    xlab <- 'Number of groups'
  } else if (sim.case == 'vary_ngrp_true') {
    xlab <- 'Number of true groups'
  } else if (sim.case == 'vary_beta') {
    xlab <- 'Magnitude of beta'
  } else if (sim.case == 'vary_sigma') {
    xlab <- 'Sigma'
  } else {
    xlab <- 'NA'
  }
  
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
  
  ## plot
  # pdf(file = paste0(date, '_', sim.case, '_obj_coef_diff.pdf'))
  # method.names <- methods[-1]
  # boxplot(res$obj.diff, names = method.names,
  #         main = 'objective difference relative to glmnet true')
  # boxplot(res$coef.diff, names = method.names,
  #         main = 'distance (L2) of coef vec to glmnet coef')
  # dev.off()
  
  ## result analysis
  # -----------------------------------------------------------------------------
  rule.name <- methods
  time.df.all <- data.frame(time = matrix(t(time.mean.all), ncol = 1, byrow = T),
                        Method = rep(rule.name, each = length(which.vary)),
                        Which.vary = rep(which.vary, length(methods)))
  time.df.all$Method <- factor(time.df.all$Method, methods)
  
  time.df.model <- data.frame(time = matrix(t(time.mean.model), ncol = 1, byrow = T),
                            Method = rep(rule.name, each = length(which.vary)),
                            Which.vary = rep(which.vary, length(methods)))
  time.df.model$Method <- factor(time.df.model$Method, methods)

  ## plot time.model without active
  time.df.model.nac <- subset(time.df.model, Method %in% c("NS-NAC", "NS-AC", "SSR-NAC", "SEDPP-NAC", "SSR-BEDPP-NAC"))
  time.df.model.nac$Method <- revalue(time.df.model.nac$Method, c("NS-NAC"="Basic GD", "NS-AC"="AC",
                                                                  "SSR-NAC"='SSR', 'SEDPP-NAC'='SEDPP',
                                                                  'SSR-BEDPP-NAC'='SSR-BEDPP'))
  
  gp.model.nac <- ggplot(time.df.model.nac, aes(x = Which.vary, y = time, color = Method)) +
    geom_line(size = 1) +
    scale_x_continuous(breaks = pretty(range(time.df.model.nac$Which.vary)),
                       limits = range(time.df.model.nac$Which.vary)) +
    scale_y_continuous(breaks = pretty(range(time.df.model.nac$time))) +
    xlab(xlab) +
    ylab("Computing time (s)") +
    theme_bw() +
    theme(legend.position = c(.2, .7)
          ,legend.key.size = unit(0.5, "cm")
          )
  
  png(filename = paste0(date, '_', sim.case, '_model_nac.png'), width = width, height = height,
      units = 'in', res = resolution)
  print(gp.model.nac)
  dev.off()

  save(list = ls(), file = paste0(date, '_', sim.case, '_results.RData'))

}

