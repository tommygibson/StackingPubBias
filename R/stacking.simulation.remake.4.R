####### First stacking simulation
## mean effect theta = 0.1
## first selection function
## moderate version
## target sample sizes of 10, 20, 40, 80

library(rstan)
library(R2jags)
library(here)
library(tidyverse)
library(loo)

options(mc.cores = 4) # we have 8 cores available but running 4 chains

summary.func <- function(x){
  c(mean(x), sd(x), quantile(x, c(.025, .975)))
}
propensity.func.1 <- function(p){
  p.prop <- ifelse(p < 0.005, 1, 
                   ifelse(p < 0.2, exp(-2 * p), 
                          ifelse(p < 0.5, exp(-4 * p), .1)))
  
  return(p.prop)
}
propensity.func.2 <- function(p){
  p.prop <- ifelse(p < 0.005, 1,
                   ifelse(p < 0.2, exp(-0.5 * p),
                          ifelse(p < 0.5, exp(-1 * p), .5)))
  return(p.prop)
}

# simulation originally run on 12/17/2021

set.seed(1217)

theta0 <- 0.1
tau <- 0.2

# initial number of studies per analysis
S.full <- read_rds(here("R", "sim.1.4.initial.S.rds"))[4,]
S.target <- c(10, 20, 40, 80)
M <- 200

t <- stan_model(here("R", "Models", "step.twoside.stan"))
t.loo <- stan_model(here("R", "Models", "step.twoside.loo.stan"))
o <- stan_model(here("R", "Models", "step.oneside.stan"))
o.loo <- stan_model(here("R", "Models", "step.oneside.loo.stan"))

for(n in 1:length(S.full)){
  
  weights <- as.data.frame(matrix(nrow = M, ncol = 8)) # weight for each model (mav, bai, step functions)
  stacked.summary <- as.data.frame(matrix(nrow = M, ncol = 4)) # mean, sd, 2.5, 97.5
  model.sums <- as.data.frame(matrix(nrow = 9 * M, ncol = 6)) # model, iteration, mean, sd, 2.5%, 97.5%
  data.summary <- as.data.frame(matrix(nrow = M, ncol = 5)) # num. studies, mean.select, mean.full, sd.select, sd.full
  big_k <- as.data.frame(matrix(nrow = M, ncol = 8))
  funnel.plots <- list()
  seeds <- list()
  
  names(stacked.summary) <- c("est.mean", "sd", "ci.lower", "ci.upper")
  names(model.sums) <- c("model", "iteration", "est.mean", "sd", "ci.lower", "ci.upper")
  names(weights) <- c("Mavridis", "Bai", "twoside.1", "twoside.2", 
                      "oneside.1", "oneside.2", "oneside.3", "oneside.4")
  names(data.summary) <- c("Num.studies", "mean.select", "mean.full", "sd.select", "sd.full")
  names(big_k) <- names(weights)
  for(j in 1:M){
    
    # dat contains simulated y and z
    dat <- matrix(nrow = S.full[n], ncol = 3)
    # standard errors
    s <- runif(S.full[n], 0.1, 0.8)
    # true study-level effects
    theta <- rnorm(S.full[n], mean = theta0, sd = tau)
    
    # simulate y and calculate p-values
    for(k in 1:S.full[n]){
      dat[k, 1:2] <- c(rnorm(1, theta[k], s[k]), s[k])
      dat[k, 3] <- 1 - pnorm(dat[k, 1] / s[k])
    }
    
    # complete data
    dat.full <- as.data.frame(dat)
    # names(dat.full) <- c('y', 'z', 's')
    names(dat.full) <- c('y', 's', 'p')
    
    # observed studies
    
    dat.full$propensity <- with(dat.full,
                                propensity.func.2(p))
    dat.full$select <- NA
    for(k in 1:S.full[n]){
      dat.full$select[k] <- rbinom(1, 1, dat.full$propensity[k]) # selection indicator
    }
    #subset data by selection indicator
    dat.select <- dat.full %>%
      filter(select == 1)
    # prepare data for jags/stan
    S <- dim(dat.select)[1]
    y <- dat.select$y
    s <- dat.select$s
    p <- dat.select$p
    
    L1 <- 0
    L2 <- 0.5
    U1 <- 0.5
    U2 <- 1
    
    # summarize selected data and full data
    data.summary[j,] <- c(S, mean(dat.select$y), mean(dat.full$y), mean(dat.select$s), mean(dat.full$s))
    if(j <= 10) {
      # save a few funnel plots
      funnel.plots[[j]] <- ggplot(dat.full, aes(x = y, y = s, color = (select == 1))) +
        geom_point()
    }
    
    # initial values
    init.gen.std <- function(){
      list(
        theta0 = rnorm(1, 0, 0.25),
        tau = runif(1, 0.1, 0.4)
      )
    }
    init.gen.bai.adj <- function(){
      list(
        z = runif(S, 0, 1),
        theta0 = rnorm(1, 0, 0.25),
        rho = runif(1, -0.5, 0.5),
        tau = runif(1, 0.1, 0.5),
        gamma0 = runif(1, -1, 1),
        gamma1 = runif(1, 0, max(s))
      )
    }
    
    init.gen.mav.adj <- function(){
      list(
        z = runif(S, 0, 1),
        theta0 = rnorm(1, 0, 0.25),
        rho = runif(1, -0.5, 0.5),
        tau = runif(1, 0.1, 0.5),
        p.low = runif(1, 0.38, 0.42),
        p.high = runif(1, 0.78, 0.82)
      )
    }
    
    
    stack.params <- c("theta0", "loglik")
    std.params <- c("theta0")
    std.dat <- list(S = S,
                    y = y,
                    s = s)
    
    bai.dat <- std.dat <- list(S = S,
                               y = y,
                               s = s)
    
    mav.dat <- list(S = S,
                    y = y,
                    s = s,
                    L1 = L1, L2 = L2,
                    U1 = U1, U2 = U2)
    
    twoside.dat.1 <- list(S = S,
                          y = y,
                          s = s,
                          steps = array(c(.05), dim = 1),
                          M = 1)
    twoside.dat.2 <- list(S = S,
                          y = y,
                          s = s,
                          steps = c(.1, .01),
                          M = 2)
    
    oneside.dat.1 <- list(S = S,
                          y = y,
                          s = s,
                          steps = array(c(.025), dim = 1),
                          M = 1)
    
    oneside.dat.2 <- list(S = S,
                          y = y,
                          s = s,
                          steps = c(.5, .025),
                          M = 2)
    
    oneside.dat.3 <- list(S = S,
                          y = y,
                          s = s,
                          steps = c(.05, .005),
                          M = 2)
    oneside.dat.4 <- list(S = S,
                          y = y,
                          s = s,
                          steps = c(.1, .025),
                          M = 2)
    
    fit.std <- do.call(jags.parallel,
                       list(data = names(std.dat), parameters.to.save = std.params, inits = init.gen.std,
                            model.file = here("R", "Models", "std.meta.txt"), 
                            n.iter = 2000, n.chains = 4, DIC = FALSE))
    copas.mav <- do.call(jags.parallel,
                         list(data = names(mav.dat), inits = init.gen.mav.adj, parameters.to.save = stack.params,
                              model.file = here("R", "Models", "copas.jags.mavridis.adj.txt"),
                              n.iter = 4000, n.chains = 4, n.burnin = 3000, DIC = FALSE))
    copas.bai <- do.call(jags.parallel,
                         list(data = names(bai.dat), inits = init.gen.bai.adj, parameters.to.save = stack.params,
                              model.file = here("R", "Models", "copas.jags.bai.adj.txt"),
                              n.iter = 4000, n.chains = 4, n.burnin = 3000, DIC = FALSE))
    
    twoside.1 <- sampling(t, data = twoside.dat.1,
                          iter = 2000, chains = 4, cores = 1)
    twoside.2 <- sampling(t, data = twoside.dat.2,
                          iter = 2000, chains = 4, cores = 1)
    oneside.1 <- sampling(o, data = oneside.dat.1,
                          iter = 2000, chains = 4, cores = 1)
    oneside.2 <- sampling(o, data = oneside.dat.2,
                          iter = 2000, chains = 4, cores = 1)
    oneside.3 <- sampling(o, data = oneside.dat.3,
                          iter = 2000, chains = 4, cores = 1)
    oneside.4 <- sampling(o, data = oneside.dat.4,
                          iter = 2000, chains = 4, cores = 1)
    
    
    model.sums[(9 * j - 8):(9 * j),] <- cbind(c("std", "Mavridis", "Bai", "twoside.1", "twoside.2", 
                                                "oneside.1", "oneside.2", "oneside.3", "oneside.4"),
                                              rep(j, 9),
                                              rbind(fit.std$BUGSoutput$summary[c(1:3, 7)],
                                                    copas.mav$BUGSoutput$summary[(S + 1), c(1:3, 7)], #1-3, 7 are mean, sd, 2.5, 97.5
                                                    copas.bai$BUGSoutput$summary[(S + 1), c(1:3, 7)],
                                                    summary(twoside.1, pars = "theta")$summary[c(1, 3, 4, 8)],
                                                    summary(twoside.2, pars = "theta")$summary[c(1, 3, 4, 8)],
                                                    summary(oneside.1, pars = "theta")$summary[c(1, 3, 4, 8)],
                                                    summary(oneside.2, pars = "theta")$summary[c(1, 3, 4, 8)],
                                                    summary(oneside.3, pars = "theta")$summary[c(1, 3, 4, 8)],
                                                    summary(oneside.4, pars = "theta")$summary[c(1, 3, 4, 8)]))
    
    loglik <- list()
    loglik[[1]] <- copas.mav$BUGSoutput$sims.list$loglik
    loglik[[2]] <- copas.bai$BUGSoutput$sims.list$loglik
    loglik[[3]] <- extract_log_lik(twoside.1, parameter_name = "loglik")
    loglik[[4]] <- extract_log_lik(twoside.2, parameter_name = "loglik")
    loglik[[5]] <- extract_log_lik(oneside.1, parameter_name = "loglik")
    loglik[[6]] <- extract_log_lik(oneside.2, parameter_name = "loglik")
    loglik[[7]] <- extract_log_lik(oneside.3, parameter_name = "loglik")
    loglik[[8]] <- extract_log_lik(oneside.4, parameter_name = "loglik")
    
    r_eff <- lapply(loglik, function(x){
      relative_eff(exp(x), chain_id = rep(1:4, each = 1000)) # 4 chains, each with 1000 draws
    })
    
    loo_list <- lapply(1:length(loglik), function(j){
      loo::loo(loglik[[j]], r_eff = r_eff[[j]])
    })
    
    ### loo (and stacking) will be unreliable if there are big pareto k
    ### diagnosic values for any points i.
    ### extract indices of big pareto k's for each model
    big_k_indices <- lapply(loo_list, function(x){
      which(x$diagnostics$pareto_k > 0.7)
    })
    num_big_k <- sapply(big_k_indices, length)
    
    big_k[j,] <- num_big_k
    ### if there aren't any big pareto k's then we can just stack!
    if(sum(num_big_k) == 0){
      weights[j,] <- loo_model_weights(loo_list, method = 'stacking', r_eff_list = r_eff)
    }
    
    else{ ### matrix of pointwise loo-elpd values
      ### need these for when we have to perform exact loo
      lpds <- do.call(cbind, (lapply(1:length(loo_list), function(j){
        unlist(loo_list[[j]]$pointwise[,1])
      })))
      
      for(m in which(num_big_k != 0)){
        ### set up data for leave-one-out
        S <- dim(dat.select)[1] - 1
        H <- 1
        
        for(i in 1:num_big_k[m]){
          y_h <- dat.select$y[big_k_indices[[m]][i]]
          s_h <- dat.select$s[big_k_indices[[m]][i]]
          y <- dat.select$y[-big_k_indices[[m]][i]]
          s <- dat.select$s[-big_k_indices[[m]][i]]
          
          ## different data structure, inits, etc for each model (ugh)
          if(m == 1){
            dat.loo <- list(S = S, H = H,
                            y = y, y_h = y_h,
                            s = s, s_h = s_h,
                            L1 = L1, L2 = L2,
                            U1 = U1, U2 = U2)
            init.gen.loo <- function(){
              list(
                z = runif(S, 0, 1),
                theta0 = rnorm(1, 0, 0.25),
                rho = runif(1, -0.5, 0.5),
                tau = runif(1, 0.1, 0.5),
                p.low = runif(1, 0.38, 0.42),
                p.high = runif(1, 0.78, 0.82)
              )
            }
            
            params.loo <- c("loglik_h")
            ### Getting weird errors using jags.parallel
            fit.loo <- do.call(jags.parallel,
                               list(data = names(dat.loo), inits = init.gen.loo, parameters.to.save = params.loo,
                                    model.file = here("R", "Models", "copas.mav.loo.txt"),
                                    n.iter = 4000, n.chains = 4, n.burnin = 3000, DIC = FALSE))
            # fit.loo <- jags(data = dat.loo, inits = init.gen.loo, parameters.to.save = params.loo,
            #                 model.file = here("R", "Models", "copas.mav.loo.txt"),
            #                 n.iter = 5000, n.chains = 4, n.burnin = 4000, DIC = FALSE)
            
            # occasionally jags can't handle the log-likelihood of real bad observations
            # extract the regular likelihood instead, log it afterwards, seems to work
            elpd_holdout <- elpd(log(fit.loo$BUGSoutput$sims.array))$estimates[1]
            
            lpds[big_k_indices[[m]][i], m] <- elpd_holdout
          }
          else if(m == 2){
            dat.loo <- list(S = S, H = H,
                            y = y, y_h = y_h,
                            s = s, s_h = s_h)
            init.gen.loo <- function(){
              list(
                z = runif(S, 0, 1),
                theta0 = rnorm(1, 0, 0.25),
                rho = runif(1, -0.5, 0.5),
                tau = runif(1, 0.1, 0.5),
                gamma0 = runif(1, -1, 1),
                gamma1 = runif(1, 0, max(s))
              )
            }
            
            params.loo <- c("loglik_h")
            fit.loo <- do.call(jags.parallel,
                               list(data = names(dat.loo), inits = init.gen.loo, parameters.to.save = params.loo,
                                    model.file = here("R", "Models", "copas.bai.loo.txt"),
                                    n.iter = 4000, n.chains = 4, n.burnin = 3000, DIC = FALSE))
            # fit.loo <- jags(data = dat.loo, inits = init.gen.loo, parameters.to.save = params.loo,
            #                 model.file = here("R", "Models", "copas.bai.loo.txt"),
            #                 n.iter = 4000, n.chains = 4, n.burnin = 3000, DIC = FALSE)
            
            # occasionally jags can't handle the log-likelihood of real bad observations
            # extract the regular likelihood instead, log it afterwards, seems to work
            elpd_holdout <- elpd(log(fit.loo$BUGSoutput$sims.array))$estimates[1]
            
            lpds[big_k_indices[[m]][i], m] <- elpd_holdout  
          }
          else if(m == 3){ ## two-side selection, p = .05
            dat.loo <- list(S = S, H = H,
                            y = y, y_h = array(y_h, dim = 1),
                            s = s, s_h = array(s_h, dim = 1),
                            steps = array(c(.05), dim = 1),
                            M = 1)
            fit.loo <- sampling(t.loo, data = dat.loo,
                                iter = 2000, chains = 4, cores = 1)
            
            elpd_holdout <- elpd(extract_log_lik(fit.loo, parameter_name = "loglik_h", merge_chains = FALSE))$estimates[1]
            
            lpds[big_k_indices[[m]][i], m] <- elpd_holdout
          }
          else if(m == 4){ ## two-side, p = .01, .1
            dat.loo <- list(S = S, H = H,
                            y = y, y_h = array(y_h, dim = 1),
                            s = s, s_h = array(s_h, dim = 1),
                            steps = c(0.1, 0.01),
                            M = 2)
            fit.loo <- sampling(t.loo, data = dat.loo,
                                iter = 2000, chains = 4, cores = 1)
            
            elpd_holdout <- elpd(extract_log_lik(fit.loo, parameter_name = "loglik_h", merge_chains = FALSE))$estimates[1]
            
            lpds[big_k_indices[[m]][i], m] <- elpd_holdout
          }
          else if(m == 5){ ## one-side, p = .025
            dat.loo <- list(S = S, H = H,
                            y = y, y_h = array(y_h, dim = 1),
                            s = s, s_h = array(s_h, dim = 1),
                            steps = array(c(.025), dim = 1),
                            M = 1)
            fit.loo <- sampling(o.loo, data = dat.loo,
                                iter = 2000, chains = 4, cores = 1)
            
            elpd_holdout <- elpd(extract_log_lik(fit.loo, parameter_name = "loglik_h", merge_chains = FALSE))$estimates[1]
            
            lpds[big_k_indices[[m]][i], m] <- elpd_holdout
          }
          else if(m == 6){ ## one-side, p = .025, .5
            dat.loo <- list(S = S, H = H,
                            y = y, y_h = array(y_h, dim = 1),
                            s = s, s_h = array(s_h, dim = 1),
                            steps = c(.5, .025),
                            M = 2)
            fit.loo <- sampling(o.loo, data = dat.loo,
                                iter = 2000, chains = 4, cores = 1)
            
            elpd_holdout <- elpd(extract_log_lik(fit.loo, parameter_name = "loglik_h", merge_chains = FALSE))$estimates[1]
            
            lpds[big_k_indices[[m]][i], m] <- elpd_holdout
          }
          else if(m == 7){ ## one-side, p = .005, .05
            dat.loo <- list(S = S, H = H,
                            y = y, y_h = array(y_h, dim = 1),
                            s = s, s_h = array(s_h, dim = 1),
                            steps = c(0.05, 0.005),
                            M = 2)
            fit.loo <- sampling(o.loo, data = dat.loo,
                                iter = 2000, chains = 4, cores = 1)
            
            elpd_holdout <- elpd(extract_log_lik(fit.loo, parameter_name = "loglik_h", merge_chains = FALSE))$estimates[1]
            
            lpds[big_k_indices[[m]][i], m] <- elpd_holdout
          }
          else if(m == 8){ ## one-side, p = .025, .1
            dat.loo <- list(S = S, H = H,
                            y = y, y_h = array(y_h, dim = 1),
                            s = s, s_h = array(s_h, dim = 1),
                            steps = c(0.1, 0.025),
                            M = 2)
            fit.loo <- sampling(o.loo, data = dat.loo,
                                iter = 2000, chains = 4, cores = 1)
            
            elpd_holdout <- elpd(extract_log_lik(fit.loo, parameter_name = "loglik_h", merge_chains = FALSE))$estimates[1]
            
            lpds[big_k_indices[[m]][i], m] <- elpd_holdout
          }
          
        }
      }
      
      ### having re-calculated elpd for problematic points, we calculate stacking weights
      ### using an S x M matrix of elpd's rather than with loo objects
      weights[j,] <- stacking_weights(lpds)
    }
    
    num.sims <- round(unlist(weights[j,]) * 4000)
    sims <- list()
    
    sims[[1]] <- sample(copas.mav$BUGSoutput$sims.list$theta0, size = num.sims[1], replace = TRUE)
    sims[[2]] <- sample(copas.bai$BUGSoutput$sims.list$theta0, num.sims[2], replace = TRUE)
    sims[[3]] <- sample(rstan::extract(twoside.1)$theta, num.sims[3], replace = TRUE)
    sims[[4]] <- sample(rstan::extract(twoside.2)$theta, num.sims[4], replace = TRUE)
    sims[[5]] <- sample(rstan::extract(oneside.1)$theta, num.sims[5], replace = TRUE)
    sims[[6]] <- sample(rstan::extract(oneside.2)$theta, num.sims[6], replace = TRUE)
    sims[[7]] <- sample(rstan::extract(oneside.3)$theta, num.sims[7], replace = TRUE)
    sims[[8]] <- sample(rstan::extract(oneside.4)$theta, num.sims[8], replace = TRUE)
    sims <- unlist(sims)
    # summary
    stacked.summary[j, ] <- summary.func(sims)
    
    if(j > 10){
      rm(sims, loo_list, r_eff, loglik)
      gc()
    }
    # seeds[[j]] <- .Random.seed
    
  }
  
  func1.moderate.small <- list(stacked.summary, model.sums, weights, theta0, data.summary, funnel.plots, big_k)
  names(func1.moderate.small) <- c("stacked", "models", "weights", "theta0", "data.summary", "funnel.plots", "big.k")
  saveRDS(func1.moderate.small, file = here("R", "Results", paste("func1.moderate.small.", S.target[n], ".rds", sep = "")))
  
}




