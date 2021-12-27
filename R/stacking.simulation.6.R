####### Sixth stacking simulation
## selection function is "moderate"
## same general form as simulations 1-3

library(rstan)
library(R2jags)
library(here)
library(tidyverse)
library(loo)

options(mc.cores = 4) # we have 8 cores available but running 4 chains

summary.func <- function(x){
  c(mean(x), sd(x), quantile(x, c(.025, .975)))
}
propensity.func.2 <- function(p){
  p.prop <- ifelse(p < 0.005, 1,
                   ifelse(p < 0.2, exp(-0.5 * p),
                          ifelse(p < 0.5, exp(-1 * p), .5)))
  return(p.prop)
}

set.seed(1222)

theta0 <- 0.4
tau <- 0.2

S.full <- 70
M <- 100

weights <- matrix(nrow = M, ncol = 5) # weight for each model (mav, bai, 1-step, 2-step, 3-step)
stacked.summary <- matrix(nrow = M, ncol = 4) # mean, sd, 2.5, 97.5
model.sums <- matrix(nrow = 6 * M, ncol = 6) # model, iteration, mean, sd, 2.5%, 97.5%
data.summary <- matrix(nrow = M, ncol = 5) # num. studies, mean.select, mean.full, sd.select, sd.full
funnel.plots <- list()
seeds <- list()
big_k <- matrix(nrow = M, ncol = 5)


for(j in 1:M){
  
  # dat contains simulated y and z
  dat <- matrix(nrow = S.full, ncol = 3)
  # standard errors
  s <- runif(S.full, 0.2, 0.8)
  # true study-level effects
  theta <- rnorm(S.full, mean = theta0, sd = tau)
  
  # simulate y and calculate p-values
  for(k in 1:S.full){
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
  for(k in 1:S.full){
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
  step1.dat <- list(S = S,
                    y = y,
                    s = s,
                    steps = array(c(.05), dim = 1),
                    M = 1)
  step2.dat <- list(S = S,
                    y = y,
                    s = s,
                    steps = c(.1, .01),
                    M = 2)
  step2.dat.2 <- list(S = S,
                      y = y,
                      s = s,
                      steps = c(.2, .05),
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
  
  # copas.mav <- stan(file = here("R", "mav.stan"), data = mav.dat,
  #                   iter = 5000, warmup = 3000)
  # copas.bai <- stan(file = here("R", "bai.2.stan"), data = bai.dat,
  #                   iter = 5000, warmup = 3000)
  
  step.1 <- stan(file = here("R", "Models", "stepfunc.stan"), data = step1.dat,
                 iter = 2000, chains = 4)
  step.2 <- stan(file = here("R", "Models", "stepfunc.stan"), data = step2.dat,
                 iter = 2000, chains = 4)
  step.3 <- stan(file = here("R", "Models", "stepfunc.stan"), data = step2.dat.2,
                 iter = 2000, chains = 4)
  
  model.sums[(6 * j - 5):(6 * j),] <- cbind(c("std", "Mavridis", "Bai", "one.step", "two.step", "three.step"),
                                            rep(j, 6),
                                            rbind(fit.std$BUGSoutput$summary[c(1:3, 7)],
                                                  copas.mav$BUGSoutput$summary[(S + 1), c(1:3, 7)], #1-3, 7 are mean, sd, 2.5, 97.5
                                                  copas.bai$BUGSoutput$summary[(S + 1), c(1:3, 7)],
                                                  summary(step.1, pars = "theta")$summary[c(1, 3, 4, 8)],
                                                  summary(step.2, pars = "theta")$summary[c(1, 3, 4, 8)],
                                                  summary(step.3, pars = "theta")$summary[c(1, 3, 4, 8)]))
  
  loglik <- list()
  loglik[[1]] <- copas.mav$BUGSoutput$sims.list$loglik
  loglik[[2]] <- copas.bai$BUGSoutput$sims.list$loglik
  loglik[[3]] <- extract_log_lik(step.1, parameter_name = "loglik")
  loglik[[4]] <- extract_log_lik(step.2, parameter_name = "loglik")
  loglik[[5]] <- extract_log_lik(step.3, parameter_name = "loglik")
  
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
  
  else{
    lpds <- matrix(nrow = S, ncol = length(loo_list))
    
    for(m in 1:length(loo_list)){
      lpds[,m] <- unlist(loo_list[[m]]$pointwise[,1])
      
      if(num_big_k[m] == 0){
        next
      }
      else{
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
          else if(m == 3){
            dat.loo <- list(S = S, H = H,
                            y = y, y_h = array(y_h, dim = 1),
                            s = s, s_h = array(s_h, dim = 1),
                            steps = array(c(.05), dim = 1),
                            M = 1)
            fit.loo <- stan(file = here("R", "Models", "stepfunc_loo.stan"), data = dat.loo,
                            iter = 2000, chains = 4, cores = 1)
            
            elpd_holdout <- elpd(extract_log_lik(fit.loo, parameter_name = "loglik_h", merge_chains = FALSE))$estimates[1]
            
            lpds[big_k_indices[[m]][i], m] <- elpd_holdout
          }
          else if(m == 4){
            dat.loo <- list(S = S, H = H,
                            y = y, y_h = array(y_h, dim = 1),
                            s = s, s_h = array(s_h, dim = 1),
                            steps = c(0.01, 0.10),
                            M = 2)
            fit.loo <- stan(file = here("R", "Models", "stepfunc_loo.stan"), data = dat.loo,
                            iter = 2000, chains = 4, cores = 1)
            
            elpd_holdout <- elpd(extract_log_lik(fit.loo, parameter_name = "loglik_h", merge_chains = FALSE))$estimates[1]
            
            lpds[big_k_indices[[m]][i], m] <- elpd_holdout
          }
          else if(m == 5){
            dat.loo <- list(S = S, H = H,
                            y = y, y_h = array(y_h, dim = 1),
                            s = s, s_h = array(s_h, dim = 1),
                            steps = c(0.05, 0.20),
                            M = 2)
            fit.loo <- stan(file = here("R", "Models", "stepfunc_loo.stan"), data = dat.loo,
                            iter = 2000, chains = 4, cores = 1)
            
            elpd_holdout <- elpd(extract_log_lik(fit.loo, parameter_name = "loglik_h", merge_chains = FALSE))$estimates[1]
            
            lpds[big_k_indices[[m]][i], m] <- elpd_holdout
          }
          
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
  sims[[3]] <- sample(rstan::extract(step.1)$theta, num.sims[3], replace = TRUE)
  sims[[4]] <- sample(rstan::extract(step.2)$theta, num.sims[4], replace = TRUE)
  sims[[5]] <- sample(rstan::extract(step.3)$theta, num.sims[5], replace = TRUE)
  sims <- unlist(sims)
  # summary
  stacked.summary[j, ] <- summary.func(sims)
  
  if(j > 10){
    rm(sims, loo_list, r_eff, loglik)
    gc()
  }
  seeds[[j]] <- .Random.seed
}

# summaryize results
rm(sims, loo_list, r_eff, loglik)
gc()
stacked.summary <- as.data.frame(stacked.summary)
model.sums <- as.data.frame(model.sums)
weights <- as.data.frame(weights)
data.summary <- as.data.frame(data.summary)
big_k <- as.data.frame(big_k)
names(stacked.summary) <- c("est.mean", "sd", "ci.lower", "ci.upper")
names(model.sums) <- c("model", "iteration", "est.mean", "sd", "ci.lower", "ci.upper")
names(weights) <- c("mav", "bai", "one.step", "two.step", "three.step")
names(data.summary) <- c("Num.studies", "mean.select", "mean.full", "sd.select", "sd.full")
names(big_k) <- names(weights)

sim.results.6 <- list(stacked.summary, model.sums, weights, theta0, data.summary, funnel.plots, big_k)
names(sim.results.6) <- c("stacked", "models", "weights", "theta0", "data.summary", "funnel.plots", "big.k")
saveRDS(sim.results.6, file = here("R", "Results", "sim.results.6.rds"))



model.sums %>%
  filter(iteration %in% which(!is.na(stacked.summary[,1]))) %>%
  group_by(model) %>%
  summarize(m = mean(as.numeric(est.mean)),
            s = sd(as.numeric(est.mean)),
            ms = mean(as.numeric(sd)),
            cover = sum(ci.lower < 0.4 & ci.upper > 0.4),
            length = mean(as.numeric(ci.upper) - as.numeric(ci.lower)),
            rmse = sqrt(mean((as.numeric(est.mean) - 0.4)^2)))




