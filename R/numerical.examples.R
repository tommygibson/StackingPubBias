## Numerical examples using stacking

library(rstan)
library(R2jags)
library(here)
library(tidyverse)
library(loo)
library(RobustBayesianCopas)
library(gghighlight)
library(metafor)

### Data setup for each model

stack_analysis <- function(y, s){
  
  S <- S.full <- length(y)
  y.full <- y
  s.full <- s
  bai.dat <- std.dat <- list(S = S,
                             s = s,
                             y = y)
  mav.dat <- list(S = S,
                  s = s,
                  y = y,
                  L1 = 0, L2 = 0.5,
                  U1 = 0.5, U2 = 1)
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
  
  # Inits for copas models run in jags
  
  init.gen.std <- function(){
    list(
      theta0 = rnorm(1, 0, 0.25),
      tau = runif(1, 0.1, 0.4)
    )
  }
  init.gen.bai <- function(){
    list(
      z = runif(S, 0, 1),
      theta0 = rnorm(1, 0, 0.25),
      rho = runif(1, -0.5, 0.5),
      tau = runif(1, 0.1, 0.5),
      gamma0 = runif(1, -1, 1),
      gamma1 = runif(1, 0, max(s))
    )
  }
  
  init.gen.mav <- function(){
    list(
      z = runif(S, 0, 1),
      theta0 = rnorm(1, 0, 0.25),
      rho = runif(1, -0.5, 0.5),
      tau = runif(1, 0.1, 0.5),
      p.low = runif(1, 0.1, 0.4),
      p.high = runif(1, 0.6, 0.9)
    )
  }
  
  params <- c("loglik", "theta0")
  
  std <- do.call(jags.parallel,
                 list(data = names(bai.dat), inits = init.gen.std, 
                      parameters.to.save = c("theta0"), model.file = here("R", "Models", "std.meta.txt"),
                      n.iter = 10000, n.burnin = 8000, n.thin = 1, n.chains = 4, DIC = FALSE))
  bai <- do.call(jags.parallel,
                 list(data = bai.dat, inits = init.gen.bai, 
                      parameters.to.save = params, model.file = here("R", "Models", "copas.jags.bai.adj.txt"),
                      n.iter = 20000, n.burnin = 16000, n.thin = 2, n.chains = 4, DIC = FALSE))
  mav <- do.call(jags.parallel,
                 list(data = mav.dat, inits = init.gen.mav, 
                      parameters.to.save = params, model.file = here("R", "Models", "copas.jags.mavridis.adj.txt"),
                      n.iter = 20000, n.burnin = 16000, n.thin = 2, n.chains = 4, DIC = FALSE))
  
  twoside.1 <- stan(file = here("R", "Models", "step.twoside.stan"), data = twoside.dat.1,
                    iter = 4000, chains = 4)
  twoside.2 <- stan(file = here("R", "Models", "step.twoside.stan"), data = twoside.dat.2,
                    iter = 4000, chains = 4)
  oneside.1 <- stan(file = here("R", "Models", "step.twoside.stan"), data = oneside.dat.1,
                    iter = 4000, chains = 4)
  oneside.2 <- stan(file = here("R", "Models", "step.twoside.stan"), data = oneside.dat.2,
                    iter = 4000, chains = 4)
  oneside.3 <- stan(file = here("R", "Models", "step.twoside.stan"), data = oneside.dat.3,
                    iter = 4000, chains = 4)
  oneside.4 <- stan(file = here("R", "Models", "step.twoside.stan"), data = oneside.dat.4,
                    iter = 4000, chains = 4)
  
  
  loglik <- list()
  loglik[[1]] <- mav$BUGSoutput$sims.list$loglik
  loglik[[2]] <- bai$BUGSoutput$sims.list$loglik
  loglik[[3]] <- extract_log_lik(twoside.1, parameter_name = "loglik")
  loglik[[4]] <- extract_log_lik(twoside.2, parameter_name = "loglik")
  loglik[[5]] <- extract_log_lik(oneside.1, parameter_name = "loglik")
  loglik[[6]] <- extract_log_lik(oneside.2, parameter_name = "loglik")
  loglik[[7]] <- extract_log_lik(oneside.3, parameter_name = "loglik")
  loglik[[8]] <- extract_log_lik(oneside.4, parameter_name = "loglik")
  
  
  r_eff <- lapply(loglik, function(x){
    relative_eff(exp(x), chain_id = rep(1:4, each = 2000)) # 4 chains, each with 1000 draws
  })
  
  loo_list <- lapply(1:length(loglik), function(j){
    loo::loo(loglik[[j]], r_eff = r_eff[[j]])
  })
  
  big_k_indices <- lapply(loo_list, function(x){
    which(x$diagnostics$pareto_k > 0.5)
  })
  num_big_k <- sapply(big_k_indices, length)
  
  ### if there aren't any big pareto k's then we can just stack!
  if(sum(num_big_k) == 0){
    weights[j,] <- loo_model_weights(loo_list, method = 'stacking', r_eff_list = r_eff)
  }
  
  else{
    lpds <- do.call(cbind, (lapply(1:length(loo_list), function(j){
      unlist(loo_list[[j]]$pointwise[,1])
    })))
    for(m in which(num_big_k != 0)){
      lpds[,m] <- unlist(loo_list[[m]]$pointwise[,1])

      ### set up data for leave-one-out
      S <- S.full - 1
      H <- 1
      
      for(i in 1:num_big_k[m]){
        y_h <- y.full[big_k_indices[[m]][i]]
        s_h <- s.full[big_k_indices[[m]][i]]
        y <- y.full[-big_k_indices[[m]][i]]
        s <- s.full[-big_k_indices[[m]][i]]
        
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
              p.low = runif(1, 0.1, 0.4),
              p.high = runif(1, 0.6, 0.9)
            )
          }
          
          params.loo <- c("loglik_h")
          ### Getting weird errors using jags.parallel
          fit.loo <- do.call(jags.parallel,
                             list(data = names(dat.loo), inits = init.gen.loo, parameters.to.save = params.loo,
                                  model.file = here("R", "Models", "copas.mav.loo.txt"),
                                  n.iter = 20000, n.chains = 4, n.thin = 2, n.burnin = 16000, DIC = FALSE))
          
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
                                  n.iter = 20000, n.chains = 4, n.thin = 2, n.burnin = 16000, DIC = FALSE))
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
          fit.loo <- stan(file = here("R", "Models", "step.twoside.loo.stan"), data = dat.loo,
                          iter = 4000, chains = 4, cores = 1)
          
          elpd_holdout <- elpd(extract_log_lik(fit.loo, parameter_name = "loglik_h", merge_chains = FALSE))$estimates[1]
          
          lpds[big_k_indices[[m]][i], m] <- elpd_holdout
        }
        else if(m == 4){ ## two-side, p = .01, .1
          dat.loo <- list(S = S, H = H,
                          y = y, y_h = array(y_h, dim = 1),
                          s = s, s_h = array(s_h, dim = 1),
                          steps = c(0.1, 0.01),
                          M = 2)
          fit.loo <- stan(file = here("R", "Models", "step.twoside.loo.stan"), data = dat.loo,
                          iter = 4000, chains = 4, cores = 1)
          
          elpd_holdout <- elpd(extract_log_lik(fit.loo, parameter_name = "loglik_h", merge_chains = FALSE))$estimates[1]
          
          lpds[big_k_indices[[m]][i], m] <- elpd_holdout
        }
        else if(m == 5){ ## one-side, p = .025
          dat.loo <- list(S = S, H = H,
                          y = y, y_h = array(y_h, dim = 1),
                          s = s, s_h = array(s_h, dim = 1),
                          steps = array(c(.025), dim = 1),
                          M = 1)
          fit.loo <- stan(file = here("R", "Models", "step.oneside.loo.stan"), data = dat.loo,
                          iter = 4000, chains = 4, cores = 1)
          
          elpd_holdout <- elpd(extract_log_lik(fit.loo, parameter_name = "loglik_h", merge_chains = FALSE))$estimates[1]
          
          lpds[big_k_indices[[m]][i], m] <- elpd_holdout
        }
        else if(m == 6){ ## one-side, p = .025, .5
          dat.loo <- list(S = S, H = H,
                          y = y, y_h = array(y_h, dim = 1),
                          s = s, s_h = array(s_h, dim = 1),
                          steps = c(.5, .025),
                          M = 2)
          fit.loo <- stan(file = here("R", "Models", "step.oneside.loo.stan"), data = dat.loo,
                          iter = 4000, chains = 4, cores = 1)
          
          elpd_holdout <- elpd(extract_log_lik(fit.loo, parameter_name = "loglik_h", merge_chains = FALSE))$estimates[1]
          
          lpds[big_k_indices[[m]][i], m] <- elpd_holdout
        }
        else if(m == 7){ ## one-side, p = .005, .05
          dat.loo <- list(S = S, H = H,
                          y = y, y_h = array(y_h, dim = 1),
                          s = s, s_h = array(s_h, dim = 1),
                          steps = c(0.05, 0.005),
                          M = 2)
          fit.loo <- stan(file = here("R", "Models", "step.oneside.loo.stan"), data = dat.loo,
                          iter = 4000, chains = 4, cores = 1)
          
          elpd_holdout <- elpd(extract_log_lik(fit.loo, parameter_name = "loglik_h", merge_chains = FALSE))$estimates[1]
          
          lpds[big_k_indices[[m]][i], m] <- elpd_holdout
        }
        else if(m == 8){ ## one-side, p = .025, .1
          dat.loo <- list(S = S, H = H,
                          y = y, y_h = array(y_h, dim = 1),
                          s = s, s_h = array(s_h, dim = 1),
                          steps = c(0.1, 0.025),
                          M = 2)
          fit.loo <- stan(file = here("R", "Models", "step.oneside.loo.stan"), data = dat.loo,
                          iter = 4000, chains = 4, cores = 1)
          
          elpd_holdout <- elpd(extract_log_lik(fit.loo, parameter_name = "loglik_h", merge_chains = FALSE))$estimates[1]
          
          lpds[big_k_indices[[m]][i], m] <- elpd_holdout
        }
          
      }

      
    }
    ### having re-calculated elpd for problematic points, we calculate stacking weights
    ### using an S x M matrix of elpd's rather than with loo objects
    weights <- stacking_weights(lpds)
  }
  
  names(weights) = c("Mavridis", "Bai", "twoside.1", "twoside.2", "oneside.1", "oneside.2", "oneside.3", "oneside.4")
  # 
  num.sims <- round(unlist(weights) * 8000)
  sims <- list()
  
  # take w_k * T samples from each model for the stacked posterior
  sims[[1]] <- sample(mav$BUGSoutput$sims.list$theta0, size = num.sims[1], replace = TRUE)
  sims[[2]] <- sample(bai$BUGSoutput$sims.list$theta0, num.sims[2], replace = TRUE)
  sims[[3]] <- sample(rstan::extract(twoside.1)$theta, num.sims[3], replace = TRUE)
  sims[[4]] <- sample(rstan::extract(twoside.2)$theta, num.sims[4], replace = TRUE)
  sims[[5]] <- sample(rstan::extract(oneside.1)$theta, num.sims[5], replace = TRUE)
  sims[[6]] <- sample(rstan::extract(oneside.2)$theta, num.sims[6], replace = TRUE)
  sims[[7]] <- sample(rstan::extract(oneside.3)$theta, num.sims[7], replace = TRUE)
  sims[[8]] <- sample(rstan::extract(oneside.4)$theta, num.sims[8], replace = TRUE)
  sims <- unlist(sims)
  
  # summary
  models.sums <- cbind.data.frame(c("std", "Mavridis", "Bai", "twoside.1", "twoside.2", "oneside.1", "oneside.2", "oneside.3", "oneside.4", "stacked"),
                            rbind(std$BUGSoutput$summary[c(1:3, 7)],
                                  mav$BUGSoutput$summary[(S.full + 1), c(1:3, 7)], #1-3, 7 are mean, sd, 2.5, 97.5
                                  bai$BUGSoutput$summary[(S.full + 1), c(1:3, 7)],
                                  summary(twoside.1, pars = "theta")$summary[c(1, 3, 4, 8)],
                                  summary(twoside.2, pars = "theta")$summary[c(1, 3, 4, 8)],
                                  summary(oneside.1, pars = "theta")$summary[c(1, 3, 4, 8)],
                                  summary(oneside.2, pars = "theta")$summary[c(1, 3, 4, 8)],
                                  summary(oneside.3, pars = "theta")$summary[c(1, 3, 4, 8)],
                                  summary(oneside.4, pars = "theta")$summary[c(1, 3, 4, 8)],
                                  summary.func(sims)))
  names(models.sums) <- c("Model", "Mean", "SD", "2.5%", "97.5%")
  
  # posterior samples for each model
  all.sims <- cbind.data.frame(c(std$BUGSoutput$sims.list$theta0, mav$BUGSoutput$sims.list$theta0, bai$BUGSoutput$sims.list$theta0,
                                 rstan::extract(twoside.1)$theta, rstan::extract(twoside.2)$theta, rstan::extract(oneside.1)$theta,
                                 rstan::extract(oneside.2)$theta, rstan::extract(oneside.3)$theta, rstan::extract(oneside.4)$theta,
                                 sims),
                               c(rep(c("standard", "Mavridis", "Bai", "twoside.1", "twoside.2", 
                                     "oneside.1", "oneside.2", "oneside.3", "oneside.4"), each = length(mav$BUGSoutput$sims.list$theta0)),
                                 rep("stacked", length(sims))))
                    
  names(all.sims) <- c("theta", "model")
  
  full.anal <- list(models.sums, weights, all.sims)
  names(full.anal) <- c("Summary", "Weights", "posterior.samples")
  return(full.anal)
}

#### numerical example 1: hackshaw data on lung cancer from second-hand smoking
dat1 <- Hackshaw1997
y1 <- dat1$log_OR
s1 <- dat1$SE
S1 <- length(y1)


numex.1 <- stack_analysis(y1, s1)

############ numerical example 2: odds of grant by gender
dat.2 <- dat.bornmann2007 %>%
  mutate(a = maward,
         b = waward,
         c = mtotal - maward,
         d = wtotal - waward,
         nm = maward,
         nw = waward,
         log.or = log((a * d) / (b * c)),
         se.log.or = sqrt(1 / a + 1 / b + 1 / c + 1 / d))
y2 <- dat.2$log.or
s2 <- dat.2$se.log.or


numex.2 <- stack_analysis(y2, s2)


############ numerical example 3: landenberger 2005, odds of recidivism with/without cognitive behavioral therapy

dat.3 <- dat.landenberger2005 %>%
  mutate(log.or = log((n.cbt.non * n.ctrl.rec) / (n.cbt.rec * n.ctrl.non)),
         se.log.or = sqrt(1 / n.cbt.non + 1 / n.ctrl.rec + 1 / n.cbt.rec + 1 / n.ctrl.non))
y3 <- dat.3$log.or
s3 <- dat.3$se.log.or

numex.3 <- stack_analysis(y3, s3)

## plots

hackshaw.post <- numex.1$posterior.samples %>%
  mutate(weight = c(rep(numex.1$Weights, each = 8000), rep(1, 8000))) %>%
  ggplot(aes(x = theta, color = model, linetype = (model == "stacked"), size = model == "stacked")) + 
  geom_density() +
  gghighlight(weight > 0.01) +
  xlim(c(-0.15, 0.4)) +
  scale_linetype_manual(values = c(1, 2)) +
  scale_size_manual(values = c(.25, .75)) +
  scale_color_manual(values = c("blue", "green", "red"),
                     labels = c("Mavridis", "Stacked", "Two-side (1)")) +
  theme_bw() +
  labs(title = "Log-relative risk of cancer from second-hand smoke",
       x = "log(RR)") +
  guides(size = FALSE, linetype = FALSE)

bornmann.post <- numex.2$posterior.samples %>%
  mutate(weight = c(rep(numex.2$Weights, each = 8000), rep(1, 8000))) %>%
  ggplot(aes(x = theta, color = model, linetype = (model == "stacked"), size = (model == "stacked"))) + 
  geom_density() +
  gghighlight(weight > 0.01) +
  theme_bw() +
  scale_linetype_manual(values = c(1, 2)) +
  scale_size_manual(values = c(.25, .75)) +
  scale_color_manual(name = "Model", values = c("blue", "red", "turquoise1", "green"),
                     labels = c("Bai", "Mavridis", "One-side (2)", "Stacked")) +
  labs(title = "Posterior distributions for log-odds ratio of grant acceptance by gender",
       x = "log(OR)") +
  guides(size = FALSE, linetype = FALSE)

landenberger.post <- numex.3$posterior.samples %>%
  mutate(weight = c(rep(numex.2$Weights, each = 8000), rep(1, 8000))) %>%
  ggplot(aes(x = theta, color = model, linetype = (model == "stacked"), size = (model == "stacked"))) + 
  geom_density() +
  gghighlight(weight > 0.01) +
  theme_bw() +
  scale_linetype_manual(values = c(1, 2)) +
  scale_size_manual(values = c(.25, .75)) +
  scale_color_manual(values = c("blue", "red", "turquoise1", "green"),
                     labels = c("Bai", "Mavridis", "One-side (2)", "Stacked")) +
  labs(title = "Log-odds ratio of recidivism by CBT/no intervention",
       x = "log(OR)") +
  guides(size = FALSE, linetype = FALSE)

## save plots
  ggsave(here("R", "Results", "hackshaw_post.pdf"), plot = hackshaw.post,
         width = 5, height = 4, units = "in")
  ggsave(here("R", "Results", "bornmann_post.pdf"), plot = bornmann.post,
         width = 5, height = 4, units = "in")
  ggsave(here("R", "Results", "landenberger_post.pdf"), plot = landenberger.post,
         width = 5, height = 4, units = "in")
  

