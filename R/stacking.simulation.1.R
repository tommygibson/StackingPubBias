####### First stacking simulation
## selection function is just by p-values

library(rstan)
library(R2jags)
library(here)
library(tidyverse)
library(loo)

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
  p.prop <- ifelse(p < 0.5, 1, 0.3)
  
  return(p.prop)
}

set.seed(1212)

theta0 <- 0.4
tau <- 0.2

S <- 30
M <- 1000

weights <- matrix(nrow = M, ncol = 5) # weight for each model (mav, bai, 1-step, 2-step, 3-step)
stacked.summary <- matrix(nrow = M, ncol = 4) # mean, sd, 2.5, 97.5
model.sums <- matrix(nrow = 6 * M, ncol = 6) # model, iteration, mean, sd, 2.5%, 97.5%
data.summary <- matrix(nrow = M, ncol = 5) # num. studies, mean.select, mean.full, sd.select, sd.full
funnel.plots <- list()
seeds <- list()

for(j in 1:M){
  
  # dat contains simulated y and z
  dat <- matrix(nrow = S, ncol = 3)
  # standard errors
  s <- runif(S, 0.2, 0.8)
  # true study-level effects
  theta <- rnorm(S, mean = theta0, sd = tau)
  
  # simulate y and calculate p-values
  for(k in 1:S){
    dat[k, 1:2] <- c(rnorm(1, theta[k], s[k]), s[k])
    dat[k, 3] <- 1 - pnorm(dat[k, 1] / s[k])
  }
  
  # complete data
  dat.full <- as.data.frame(dat)
  # names(dat.full) <- c('y', 'z', 's')
  names(dat.full) <- c('y', 's', 'p')
  
  # observed studies

  dat.full$propensity <- with(dat.full,
                              propensity.func.1(p))
  dat.full$select <- NA
  for(k in 1:S){
    dat.full$select[k] <- rbinom(1, 1, dat.full$propensity[k]) # selection indicator
  }
  #subset data by selection indicator
  dat.select <- dat.full %>%
    filter(select == 1)
  # prepare data for jags/stan
  S.select <- dim(dat.select)[1]
  y <- dat.select$y
  s <- dat.select$s
  p <- dat.select$p
  
  # summarize selected data and full data
  data.summary[j,] <- c(S.select, mean(dat.select$y), mean(dat.full$y), mean(dat.select$s), mean(dat.full$s))
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
      z = runif(S.select, 0, 1),
      theta0 = rnorm(1, 0, 0.25),
      rho = runif(1, -0.5, 0.5),
      tau = runif(1, 0.1, 0.5),
      gamma0 = runif(1, -1, 1),
      gamma1 = runif(1, 0, max(s))
    )
  }
  
  init.gen.mav.adj <- function(){
    list(
      z = runif(S.select, 0, 1),
      theta0 = rnorm(1, 0, 0.25),
      rho = runif(1, -0.5, 0.5),
      tau = runif(1, 0.1, 0.5),
      p.low = runif(1, 0.38, 0.42),
      p.high = runif(1, 0.78, 0.82)
    )
  }
  
  
  stack.params <- c("theta0", "loglik")
  std.params <- c("theta0")
  std.dat <- list(S = S.select,
                  y = y,
                  s = s)
  
  bai.dat <- log.p.dat <- std.dat <- list(S = S.select,
                                          y = y,
                                          s = s)
  
  mav.dat <- list(S = S.select,
                  y = y,
                  s = s,
                  L1 = 0.01, L2 = 0.49,
                  U1 = 0.51, U2 = 0.99)
  step1.dat <- list(S = S.select,
                    y = y,
                    s = s,
                    steps = array(c(.05), dim = 1),
                    M = 1)
  step2.dat <- list(S = S.select,
                    y = y,
                    s = s,
                    steps = c(.1, .01),
                    M = 2)
  step2.dat.2 <- list(S = S.select,
                      y = y,
                      s = s,
                      steps = c(.2, .05),
                      M = 2)
  # step3.dat <- list(S = S.select,
  #                   y = y,
  #                   s = s,
  #                   steps = c(.2, .1, .05),
  #                   M = 3)
  
  fit.std <- jags(data = std.dat, parameters.to.save = std.params, inits = init.gen.std,
                  model.file = here("R", "std.meta.txt"), 
                  n.iter = 2000, n.chains = 3, DIC = FALSE)
  copas.mav <- jags(data = mav.dat, inits = init.gen.mav.adj, parameters.to.save = stack.params,
                    model.file = here("R", "copas.jags.mavridis.adj.txt"),
                    n.iter = 2000, n.chains = 3, DIC = FALSE)
  copas.bai <- jags(data = bai.dat, inits = init.gen.bai.adj, parameters.to.save = stack.params,
                    model.file = here("R", "copas.jags.bai.adj.txt"),
                    n.iter = 2000, n.chains = 3, DIC = FALSE)
  
  # copas.sqrt.p <- jags(data = bai.dat, inits = init.gen.logp, parameters.to.save = stack.params,
  #                      model.file = here('R', 'copas.pval.sqrt.txt'),
  #                      n.iter = 10000, n.thin = 4, n.chains = 4, DIC = FALSE)
  # copas.log.p <- jags(data = log.p.dat, inits = init.gen.logp, parameters.to.save = stack.params,
  #                     model.file = here('R', 'copas.pval.log.txt'),
  #                     n.iter = 2000, n.chains = 3, DIC = FALSE)
  
  step.1 <- stan(file = here("R", "stepfunc.stan"), data = step1.dat,
                 iter = 2000, chains = 3)
  step.2 <- stan(file = here("R", "stepfunc.stan"), data = step2.dat,
                 iter = 2000, chains = 3)
  step.3 <- stan(file = here("R", "stepfunc.stan"), data = step2.dat.2,
                 iter = 2000, chains = 3)
  
  model.sums[(6 * j - 5):(6 * j),] <- cbind(c("std", "Mavridis", "Bai", "one.step", "two.step", "three.step"),
                                            rep(j, 6),
                                            rbind(fit.std$BUGSoutput$summary[c(1:3, 7)],
                                                  copas.mav$BUGSoutput$summary[(S.select + 1), c(1:3, 7)], #1-3, 7 are mean, sd, 2.5, 97.5
                                                  copas.bai$BUGSoutput$summary[(S.select + 1), c(1:3, 7)],
                                                  summary(step.1, pars = "theta")$summary[c(1, 3, 4, 8)],
                                                  summary(step.2, pars = "theta")$summary[c(1, 3, 4, 8)],
                                                  summary(step.3, pars = "theta")$summary[c(1, 3, 4, 8)]))
  
  loglik <- list()
  loglik[[1]] <- copas.mav$BUGSoutput$sims.list$loglik
  loglik[[2]] <- copas.bai$BUGSoutput$sims.list$loglik
  #loglik[[3]] <- fit.std$BUGSoutput$sims.list$loglik
  loglik[[3]] <- extract_log_lik(step.1, parameter_name = "loglik")
  loglik[[4]] <- extract_log_lik(step.2, parameter_name = "loglik")
  loglik[[5]] <- extract_log_lik(step.3, parameter_name = "loglik")
  
  r_eff <- lapply(loglik, function(x){
    relative_eff(exp(x), chain_id = rep(1:3, each = 1000)) # 3 chains, each with 1000 draws
  })
  
  loo_list <- lapply(1:length(loglik), function(j){
    loo(loglik[[j]], r_eff = r_eff[[j]], k_threshold = 0.7)
  })
  ### extract pareto_k values
  pareto_values <- vector(length = 5 * S.select)
  for(i in 1:5){
    pareto_values[((i - 1) * S.select + 1):(i * S.select)] <- loo_list[[i]]$diagnostics$pareto_k
  }
  if(sum(pareto_values > 0.7) > 0) {  #### if we can't trust the loo object then stacking isn't good
    next
  }
  weights[j,] <- loo_model_weights(loo_list, method = 'stacking', r_eff_list = r_eff)
  
  num.sims <- round(weights[j,] * 3000)
  sims <- list()
  
  sims[[1]] <- sample(copas.mav$BUGSoutput$sims.list$theta0, num.sims[1], replace = TRUE)
  sims[[2]] <- sample(copas.bai$BUGSoutput$sims.list$theta0, num.sims[2], replace = TRUE)
  # sims[[3]] <- sample(fit.std$BUGSoutput$sims.list$theta0, num.sims[3], replace = TRUE)
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
names(stacked.summary) <- c("est.mean", "sd", "ci.lower", "ci.upper")
names(model.sums) <- c("model", "iteration", "est.mean", "sd", "ci.lower", "ci.upper")
names(weights) <- c("mav", "bai", "one.step", "two.step", "three.step")

sim.results.1 <- list(stacked.summary, model.sums, weights, theta0, funnel.plots)
names(sim.results.1) <- c("stacked", "models", "weights", "theta0", "funnel.plots")
saveRDS(sim.results.1, file = here("R", "Results", "sim.results.1.rds"))



model.sums %>%
  filter(iteration %in% which(!is.na(stacked.summary[,1]))) %>%
  group_by(model) %>%
  summarize(m = mean(as.numeric(est.mean)),
            s = sd(as.numeric(est.mean)),
            ms = mean(as.numeric(sd)),
            cover = sum(ci.lower < 0.4 & ci.upper > 0.4),
            length = mean(as.numeric(ci.upper) - as.numeric(ci.lower)),
            rmse = sqrt(mean((as.numeric(est.mean) - 0.4)^2)))




