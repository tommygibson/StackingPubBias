####### simulation for stacking

library(rstan)
library(R2jags)
library(here)
library(tidyverse)
library(loo)
library(metasens)
library(RobustBayesianCopas)

# theta0 <- 0
# tau <- 0.2
# rho <- 0.9
# gamma0.true <- -0.5
# gamma1.true <- 0.3
# S <- 50

# observed effects (y), latent selection process (z)
# dat <- matrix(nrow = S, ncol = 2)
# 
# set.seed(2021)
# 
# # simulate standard errors and true means for studies
# s <- runif(S, 0.2, 0.8)
# theta <- rnorm(S, mean = theta0, sd = tau)
# 
# for(i in 1:S){
#   Sigma <- matrix(c(s[i] ^ 2, rho * s[i], rho * s[i], 1), byrow = T, nrow = 2)
#   mu <- c(theta[i], gamma0.true + gamma1.true / s[i])
#   
#   dat[i,] <- MASS::mvrnorm(n = 1, mu = mu, Sigma = Sigma)
# }
# 
# # complete data
# dat.full <- cbind.data.frame(dat, s)
# names(dat.full) <- c('y', 'z', 's')
# 
# # observed studies
# dat.select <- dat.full %>%
#   filter(z > 0)
# 
# S.select <- dim(dat.select)[1]
# 
# ggplot(dat.full, aes(x = y, y = s, color = (z > 0))) + 
#   geom_point()
# 
# # we have our data! Now let's set up stacking
# 
# gamma0.list <- seq(-1, 1, 0.25)
# gamma1.list <- seq(0, 0.5, length.out = length(gamma0.list))
# 
# gammas <- expand.grid(gamma0.list, gamma1.list)
# gamma0 <- gammas[,1]
# gamma1 <- gammas[,2]
# 
# K <- length(gamma0)
# 
# loglik <- list()
# model.sum <- matrix(nrow = K, ncol = 2)
# 
# for(k in 1:K){
#   copas.dat <- list(S = S.select,
#                     y = dat.select$y,
#                     s = dat.select$s,
#                     gamma0 = gamma0[k],
#                     gamma1 = gamma1[k])
#   copas.fit <- stan(here("R", "copas_adj.stan"), data = copas.dat,
#                     iter = 5000, thin = 2)
#   
#   model.sum[k,] <- summary(copas.fit, pars = "mu", probs = c(.025, .5, .975))$summary[c(1, 3)]
#   
#   loglik[[k]] <- extract_log_lik(copas.fit)
# }
# 
# r_eff <- lapply(loglik, function(x){
#   relative_eff(exp(x), chain_id = rep(1:4, each = 1250))
# })
# 
# stacking_weights <- loo_model_weights(loglik, method = 'stacking', r_eff_list = r_eff)
# 
# dat.stacking <- cbind.data.frame(as.vector(stacking_weights), gammas, model.sum)
# names(dat.stacking) <- c('weight', 'gamma0', 'gamma1', 'mean', 'sd')
# stacked_mean <- weighted.mean(model.sum[,1], stacking_weights)
# 
# dat.stacking %>%
#   ggplot(aes(x = gamma0, y = gamma1, size = weight, color = mean)) +
#   geom_point() +
#   scale_color_gradient(low = "green", high = "red") +
#   scale_size_continuous(range = c(0, 16)) +
#   labs(title = "Stacking weight and estimated mean") +
#   annotate("label", x = 0.6, y = 1.05, label = paste("Stacked Mean = ", stacked_mean, sep = "")) +
#   theme_bw()
# 
# 
# ##### And then try it with jags (full means full dataset)
# 
# loglik_jags <- list()
# loglik_jags_full <- list()
# model.sum.jags <- matrix(nrow = length(gamma0), ncol = 2)
# model.sum.jags.full <- matrix(nrow = length(gamma0), ncol = 2)
# 
# 
# ### does it give correct results with correct assumption of gamma0 + gamma1?
# 
# copas.dat <- list(S = S.select,
#                   y = dat.select$y,
#                   s = dat.select$s,
#                   gamma0 = gamma0.true,
#                   gamma1 = gamma1.true)
# copas.dat.full <- list(S = S,
#                   y = dat.full$y,
#                   s = dat.full$s,
#                   gamma0 = gamma0.true,
#                   gamma1 = gamma1.true)
# 
# copas.params.jags <- c("theta", "rho", "tau")
# 
# copas.jags.init <- function(){
#   list(
#     z = runif(S.select, 0, 1),
#     theta = rnorm(0, 1),
#     tau = runif(1, 0, 0.5),
#     rho = runif(1, -0.9, 0.9)
#   )
# }
# copas.jags.init.full <- function(){
#   list(
#     z = runif(S, 0, 1),
#     theta = rnorm(0, 1),
#     tau = runif(1, 0, 0.5),
#     rho = runif(1, -0.9, 0.9)
#   )
# }
# 
# copas.fit <- jags(data = copas.dat,
#                   inits = copas.jags.init,
#                   parameters.to.save = copas.params.jags,
#                   model.file = here("R", "copas_jags_adj.txt"),
#                   n.chains = 4, n.iter = 10000, n.thin = 2, DIC = FALSE)
# 
# copas.fit.full <- jags(data = copas.dat.full,
#                   inits = copas.jags.init.full,
#                   parameters.to.save = copas.params.jags,
#                   model.file = here("R", "copas_jags_adj.txt"),
#                   n.chains = 4, n.iter = 10000, n.thin = 2, DIC = FALSE)
# 
# 
# copas.params.jags <- c("theta", "loglik")
# 
# for(k in 1:K){
#   copas.dat <- list(S = S.select,
#                     y = dat.select$y,
#                     s = dat.select$s,
#                     gamma0 = gamma0[k],
#                     gamma1 = gamma1[k])
#   copas.dat.full <- list(S = S,
#                       y = dat.full$y,
#                       s = dat.full$s,
#                       gamma0 = gamma0[k],
#                       gamma1 = gamma1[k])
#   
#   copas.fit <- jags(data = copas.dat,
#                     inits = copas.jags.init,
#                     parameters.to.save = copas.params.jags,
#                     model.file = here("R", "copas_jags_adj.txt"),
#                     n.chains = 4, n.iter = 10000, n.thin = 2, DIC = FALSE)
#   
#   model.sum.jags[k,] <- copas.fit$BUGSoutput$summary[(S.select + 1), c(1, 2)]
#   
#   loglik_jags[[k]] <- as.matrix(copas.fit$BUGSoutput$sims.list$loglik)
#   
#   copas.fit.full <- jags(data = copas.dat.full,
#                          inits = copas.jags.init.full,
#                          parameters.to.save = copas.params.jags,
#                          model.file = here('R', 'copas_jags_adj.txt'),
#                          n.chains = 4, n.iter = 10000, n.thin = 2, DIC = FALSE)
#   
#   model.sum.jags.full[k,] <- copas.fit.full$BUGSoutput$summary[(S + 1), c(1, 2)]
#   
#   loglik_jags_full[[k]] <- as.matrix(copas.fit.full$BUGSoutput$sims.list$loglik)
# }
# 
# r_eff_jags <- lapply(loglik_jags, function(x){
#   relative_eff(exp(x), chain_id = rep(1:4, each = 2500))
# })
# r_eff_jags_full <- lapply(loglik_jags_full, function(x){
#   relative_eff(exp(x), chain_id = rep(1:4, each = 2500))
# })
# weights_jags <- loo_model_weights(loglik_jags, method = 'stacking', r_eff_list = r_eff_jags)
# weights_jags_full <- loo_model_weights(loglik_jags_full, method = 'stacking', r_eff_list = r_eff_jags_full)
# 
# dat.stacking.jags <- cbind.data.frame(gammas, 
#                                       as.vector(weights_jags), model.sum.jags, 
#                                       as.vector(weights_jags_full), model.sum.jags.full)
# names(dat.stacking.jags) <- c('gamma0', 'gamma1', 'weight', 'mean', 'sd', 'weight.full', 'mean.full', 'sd.full')
# stacked_mean <- round(weighted.mean(model.sum.jags[,1], weights_jags), 3)
# stacked_mean_full <- round(weighted.mean(model.sum.jags.full[,1], weights_jags_full), 3)
# 
# dat.stacking.jags %>%
#   ggplot(aes(x = gamma0, y = gamma1, size = weight, color = mean)) +
#   geom_point() +
#   scale_color_gradient(low = "green", high = "red") +
#   scale_size_continuous(range = c(0, 16)) +
#   labs(title = "Stacking weight and estimated mean") +
#   annotate("label", x = 0.6, y = 0.6, label = paste("Stacked Mean = ", stacked_mean, sep = "")) +
#   theme_bw()
# 
# dat.stacking.jags %>%
#   ggplot(aes(x = gamma0, y = gamma1, size = weight.full, color = mean.full)) +
#   geom_point() +
#   scale_color_gradient(low = "green", high = "red") +
#   scale_size_continuous(range = c(0, 16)) +
#   labs(title = "Stacking weight and estimated mean") +
#   annotate("label", x = 0.6, y = 0.6, label = paste("Stacked Mean = ", stacked_mean_full, sep = "")) +
#   theme_bw()


######### Repeated simulation

set.seed(123)

theta0 <- 0
tau <- 0.2
rho <- 0.9
gamma0.true <- -0.5
gamma1.true <- 0.3
S <- 200
M <- 100

gamma0.list <- seq(-1, 0.5, 0.25)
gamma1.list <- seq(0, 0.5, length.out = length(gamma0.list))

gammas <- expand.grid(gamma0.list, gamma1.list)
gamma0 <- gammas[,1]
gamma1 <- gammas[,2]

K <- length(gamma0)

# two columns for two values of theta0
stacked.summary.0 <- matrix(nrow = M, ncol = 4)
bai.summary <- matrix(nrow = M, ncol = 4)
mav.summary <- matrix(nrow = M, ncol = 4)
data.summary.0 <- matrix(nrow = M, ncol = 5)
funnel.plots.0 <- list()
funnel.plots.1 <- list()

summary.func <- function(x){
  c(mean(x), sd(x), quantile(x, c(.025, .975)))
}



# for(j in 1:M){
#   
#   # dat contains simulated y and z
#   dat <- matrix(nrow = S, ncol = 2)
#   # standard errors
#   s <- runif(S, 0.2, 0.8)
#   # true study-level effects
#   theta <- rnorm(S, mean = theta0, sd = tau)
#   
#   # simulate y and z from bivariate normal
#   for(k in 1:S){
#     Sigma <- matrix(c(s[k] ^ 2, rho * s[k], rho * s[k], 1), byrow = T, nrow = 2)
#     mu <- c(theta[k], gamma0.true + gamma1.true / s[k])
#     
#     dat[k,] <- MASS::mvrnorm(n = 1, mu = mu, Sigma = Sigma)
#   }
#   
#   # complete data
#   dat.full <- cbind.data.frame(dat, s)
#   names(dat.full) <- c('y', 'z', 's')
#   
#   # observed studies
#   dat.select <- dat.full %>%
#     filter(z > 0)
#   
#   S.select <- dim(dat.select)[1]
# 
#   data.summary.0[j,] <- c(S.select, mean(dat.select$y), mean(dat.full$y), mean(dat.select$s), mean(dat.full$s))
#   funnel.plots.0[[j]] <- ggplot(dat.full, aes(x = y, y = s, color = (z > 0))) + 
#     geom_point()
#   
#   copas.jags.init <- function(){
#     list(
#       z = runif(S.select, 0, 1),
#       theta = rnorm(0, 1),
#       tau = runif(1, 0, 0.5),
#       rho = runif(1, -0.9, 0.9)
#     )
#   }
#   
#   copas.params.jags <- c("theta", "loglik")
#   
#   loglik <- list()
#   model.means <- vector(length = K)
#   mean.sims <- list()
#   
#   for(k in 1:K){
#     copas.dat <- list(S = S.select,
#                       y = dat.select$y,
#                       s = dat.select$s,
#                       gamma0 = gamma0[k],
#                       gamma1 = gamma1[k])
#     copas.dat.full <- list(S = S,
#                            y = dat.full$y,
#                            s = dat.full$s,
#                            gamma0 = gamma0[k],
#                            gamma1 = gamma1[k])
#     
#     copas.fit <- jags(data = copas.dat,
#                       inits = copas.jags.init,
#                       parameters.to.save = copas.params.jags,
#                       model.file = here("R", "copas_jags_adj.txt"),
#                       n.chains = 2, n.iter = 5000, n.thin = 2, DIC = FALSE)
#     
#     mean.sims[[k]] <- copas.fit$BUGSoutput$sims.list$theta
#     loglik[[k]] <- as.matrix(copas.fit$BUGSoutput$sims.list$loglik)
#   }
#   r_eff <- lapply(loglik, function(x){
#     relative_eff(exp(x), chain_id = rep(1:2, each = 1250))
#   })
#   
#   weights <- loo_model_weights(loglik, method = 'stacking', r_eff_list = r_eff)
#   
#   # number of samples to take from each model
#   num.sims <- round(weights * 2500)
#   sims <- list()
#   for(k in 1:K){
#     sims[[k]] <- sample(mean.sims[[k]], num.sims[k], replace = FALSE)
#   }
#   sims <- unlist(sims)
#   
#   # summary
#   stacked.summary.0[j, ] <- summary.func(sims)
#   
#   ### Bai 2020 analysis
#   bai.fit <- RobustBayesianCopas(y = dat.select$y, s = dat.select$s, re.dist = "normal")
#   
#   bai.summary[j,] <- summary.func(bai.fit$theta.samples)
#   
#   ### Mavridis 2013 analysis
#   mav.dat <- list(S = S.select,
#                   y = dat.select$y,
#                   s = dat.select$s,
#                   L1 = 0.30, L2 = 0.50,
#                   U1 = 0.6, U2 = 0.9)
#   
#   init.gen.mav <- function(){
#     list(
#       z = runif(S.select, 0, 1),
#       theta = rnorm(S.select, 0, 0.25),
#       theta0 = rnorm(1, 0, 0.25),
#       rho = runif(1, -0.5, 0.5),
#       tau = runif(1, 0.1, 0.5),
#       p.low = runif(1, 0.38, 0.42),
#       p.high = runif(1, 0.78, 0.82)
#     )
#   }
#   mav.params <- c("theta0", "rho", "tau", "gamma0", "gamma1")
#   mav.fit <- jags(data = mav.dat, inits = init.gen.mav, parameters.to.save = mav.params,
#                     model.file = here("R", "copas.jags.mavridis.txt"),
#                     n.iter = 10000, n.thin = 3, n.chains = 2, n.burnin = 5000)
#   
#   mav.summary[j,] <- summary.func(mav.fit$BUGSoutput$sims.list$theta0)
#   
#   if(j < M) {
#     rm(loglik, mean.sims, mav.fit, bai.fit, copas.fit)
#   }
#   gc()
#   
# }
# 
# 
# stacked.summary <- cbind.data.frame(stacked.summary.0, data.summary.0)
# names(stacked.summary) <- c("stacked.mean", "stacked.sd", "stacked.lower", "stacked.upper",
#                             "N.select", "select.mean", "full.mean", "select.sd", "full.sd")
# stacked.summary <- list(stacked.summary, c(gamma0.true, gamma1.true), theta0, tau, rho)
# names(stacked.summary) <- c("sim.dat", "gamma", "theta0", "tau", "rho")
# 
# saveRDS(stacked.summary, here("R", "Results", "small.stack.simulation.rds"))
# simulate standard errors and true means for studies

summary.func <- function(x){
  c(mean(x), sd(x), quantile(x, c(.025, .975)))
}
propensity.func <- function(y, s, p){
  p.prop <- ifelse(p < 0.005, 1, 
                   ifelse(p < 0.2, exp(-2 * p), 
                          ifelse(p < 0.5, exp(-5 * p), .1)))
  y.prop <- 1 / (1 + exp(-(3 * y)))
  s.prop <- 1 / (1 + exp(-(2 - 4*s)))
  
  p.w <- 1/3
  y.w <- 1/3
  s.w <- 1/3
  
  return(p.w * p.prop + y.w * y.prop + s.w * s.prop)
}

set.seed(1232)

theta0 <- 0.4
tau <- 0.2
# rho <- 0.7
# gamma0.true <- -0.5
# gamma1.true <- 0.3
S <- 50
M <- 500

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
  
  # simulate y and z from bivariate normal
  # for(k in 1:S){
  #   Sigma <- matrix(c(s[k] ^ 2, rho * s[k], rho * s[k], 1), byrow = T, nrow = 2)
  #   mu <- c(theta[k], gamma0.true + gamma1.true / s[k])
  #   
  #   dat[k,] <- MASS::mvrnorm(n = 1, mu = mu, Sigma = Sigma)
  # }
  
  # simulate y and censor it based on p-values
  
  for(k in 1:S){
    dat[k, 1:2] <- c(rnorm(1, theta[k], s[k]), s[k])
    dat[k, 3] <- 1 - pnorm(abs(dat[k, 1]) / s[k])
  }
  
  # complete data
  dat.full <- as.data.frame(dat)
  # names(dat.full) <- c('y', 'z', 's')
  names(dat.full) <- c('y', 's', 'p')
  
  # observed studies
  # dat.select <- dat.full %>%
  #   filter(z > 0)
  dat.full$propensity <- with(dat.full,
                              propensity.func(y, s, p))
  dat.full$select <- NA
  for(k in 1:S){
    # dat.full$select[k] <- ifelse(dat.full$p[k] > 0.5, rbinom(1, 1, 0.2),
    #                              ifelse(dat.full$p[k] > .05, rbinom(1, 1, 0.6), 1))
    dat.full$select[k] <- rbinom(1, 1, dat.full$propensity[k])
  }
  dat.select <- dat.full %>%
    filter(select == 1)
  S.select <- dim(dat.select)[1]
  y <- dat.select$y
  s <- dat.select$s
  p <- dat.select$p
  
  data.summary[j,] <- c(S.select, mean(dat.select$y), mean(dat.full$y), mean(dat.select$s), mean(dat.full$s))
  funnel.plots[[j]] <- ggplot(dat.full, aes(x = y, y = s, color = (select == 1))) + 
    geom_point()
  
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

  init.gen.logp <- function(){
    list(
      z = runif(S.select, 0, 1),
      theta0 = rnorm(1, 0, 0.25),
      rho = runif(1, -0.5, 0.5),
      tau = runif(1, 0.1, 0.5),
      gamma0 = runif(1, -1, 1),
      gamma1 = runif(1, 0, 1 / max(-log(p)))
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
  
  fit.std <- fit.std <- jags(data = std.dat, parameters.to.save = std.params, inits = init.gen.std,
                             model.file = here("R", "std.meta.txt"), n.chains = 2, n.iter = 5000, DIC = FALSE)
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
  # loglik[[3]] <- copas.sqrt.p$BUGSoutput$sims.list$loglik
  # loglik[[3]] <- copas.log.p$BUGSoutput$sims.list$loglik
  loglik[[3]] <- extract_log_lik(step.1, parameter_name = "loglik")
  loglik[[4]] <- extract_log_lik(step.2, parameter_name = "loglik")
  loglik[[5]] <- extract_log_lik(step.3, parameter_name = "loglik")
  
  r_eff <- lapply(loglik, function(x){
    relative_eff(exp(x), chain_id = rep(1:3, each = 1000)) # 3 chains, each with 1250 draws
  })
  
  loo_list <- lapply(1:length(loglik), function(j){
    loo(loglik[[j]], r_eff = r_eff[[j]], k_threshold = 0.7)
  })
  
  weights[j,] <- loo_model_weights(loo_list, method = 'stacking', r_eff_list = r_eff)
  
  num.sims <- round(weights[j,] * 3000)
  sims <- list()
  
  sims[[1]] <- sample(copas.mav$BUGSoutput$sims.list$theta0, num.sims[1], replace = FALSE)
  sims[[2]] <- sample(copas.bai$BUGSoutput$sims.list$theta0, num.sims[2], replace = FALSE)
  # sims[[3]] <- sample(copas.sqrt.p$BUGSoutput$sims.list$theta0, num.sims[3], replace = FALSE)
  # sims[[3]] <- sample(copas.log.p$BUGSoutput$sims.list$theta0, num.sims[3], replace = FALSE)
  sims[[3]] <- sample(rstan::extract(step.1)$theta, num.sims[3], replace = FALSE)
  sims[[4]] <- sample(rstan::extract(step.2)$theta, num.sims[4], replace = FALSE)
  sims[[5]] <- sample(rstan::extract(step.3)$theta, num.sims[5], replace = FALSE)
  sims <- unlist(sims)
  # summary
  stacked.summary[j, ] <- summary.func(sims)
  
  if(j > 10){
    rm(sims, loo_list, r_eff, loglik)
    gc()
  }
  seeds[[j]] <- .Random.seed
}
rm(sims, loo_list, r_eff, loglik)
gc()
stacked.summary <- as.data.frame(stacked.summary)
model.sums <- as.data.frame(model.sums)
weights <- as.data.frame(weights)
names(stacked.summary) <- c("est.mean", "sd", "ci.lower", "ci.upper")
names(model.sums) <- c("model", "iteration", "est.mean", "sd", "ci.lower", "ci.upper")
names(weights) <- c("mav", "bai", "one.step", "two.step", "three.step")
copas.sim.firsthalf <- list(stacked.summary, model.sums, weights)
names(copas.sim.firsthalf) <- c("stacked", "models", "weights")
saveRDS(copas.sim.firsthalf, here("R", "Results", "copas.sim.firsthalf.rds"))


model.sums %>%
  group_by(model) %>%
  summarize(m = mean(as.numeric(est.mean)),
            s = sd(as.numeric(est.mean)),
            ms = mean(as.numeric(sd)),
            cover = sum(ci.lower < 0.4 & ci.upper > 0.4))
for(j in (M / 2 + 1):M){
  
  # dat contains simulated y and z
  dat <- matrix(nrow = S, ncol = 2)
  # standard errors
  s <- runif(S, 0.2, 0.8)
  # true study-level effects
  theta <- rnorm(S, mean = theta0, sd = tau)
  
  # simulate y and z from bivariate normal
  for(k in 1:S){
    Sigma <- matrix(c(s[k] ^ 2, rho * s[k], rho * s[k], 1), byrow = T, nrow = 2)
    mu <- c(theta[k], gamma0.true + gamma1.true / s[k])
    
    dat[k,] <- MASS::mvrnorm(n = 1, mu = mu, Sigma = Sigma)
  }
  
  # complete data
  dat.full <- cbind.data.frame(dat, s)
  names(dat.full) <- c('y', 'z', 's')
  
  # observed studies
  dat.select <- dat.full %>%
    filter(z > 0)
  
  S.select <- dim(dat.select)[1]
  y <- dat.select$y
  s <- dat.select$s
  p <- 2 * pnorm(-y / s)
  
  data.summary[j,] <- c(S.select, mean(dat.select$y), mean(dat.full$y), mean(dat.select$s), mean(dat.full$s))
  funnel.plots[[j]] <- ggplot(dat.full, aes(x = y, y = s, color = (z > 0))) + 
    geom_point()
  
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
  
  init.gen.logp <- function(){
    list(
      z = runif(S.select, 0, 1),
      theta0 = rnorm(1, 0, 0.25),
      rho = runif(1, -0.5, 0.5),
      tau = runif(1, 0.1, 0.5),
      gamma0 = runif(1, -1, 1),
      gamma1 = runif(1, 0, 1 / max(-log(p)))
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
                  L1 = 0.1, L2 = 0.45,
                  U1 = 0.55, U2 = 0.9)
  
  fit.std <- fit.std <- jags(data = std.dat, parameters.to.save = std.params, inits = init.gen.std,
                             model.file = here("R", "std.meta.txt"), n.chains = 2, n.iter = 5000, DIC = FALSE)
  copas.mav <- jags(data = mav.dat, inits = init.gen.mav.adj, parameters.to.save = stack.params,
                    model.file = here("R", "copas.jags.mavridis.adj.txt"),
                    n.iter = 5000, n.thin = 2, n.chains = 3, DIC = FALSE)
  copas.bai <- jags(data = bai.dat, inits = init.gen.bai.adj, parameters.to.save = stack.params,
                    model.file = here("R", "copas.jags.bai.adj.txt"),
                    n.iter = 5000, n.thin = 2, n.chains = 3, DIC = FALSE)
  
  # copas.sqrt.p <- jags(data = bai.dat, inits = init.gen.logp, parameters.to.save = stack.params,
  #                      model.file = here('R', 'copas.pval.sqrt.txt'),
  #                      n.iter = 10000, n.thin = 4, n.chains = 4, DIC = FALSE)
  copas.log.p <- jags(data = log.p.dat, inits = init.gen.logp, parameters.to.save = stack.params,
                      model.file = here('R', 'copas.pval.log.txt'),
                      n.iter = 5000, n.thin = 2, n.chains = 3, DIC = FALSE)
  
  model.sums[(4 * j - 3):(4 * j),] <- cbind(c("std", "Mavridis", "Bai", "log.p"),
                                            rep(j, 4),
                                            rbind(fit.std$BUGSoutput$summary[c(1:3, 7)],
                                                  copas.mav$BUGSoutput$summary[(S.select + 1), c(1:3, 7)], #1-3, 7 are mean, sd, 2.5, 97.5
                                                  copas.bai$BUGSoutput$summary[(S.select + 1), c(1:3, 7)],
                                                  copas.log.p$BUGSoutput$summary[(S.select + 1), c(1:3, 7)]))
  
  loglik <- list()
  loglik[[1]] <- copas.mav$BUGSoutput$sims.list$loglik
  loglik[[2]] <- copas.bai$BUGSoutput$sims.list$loglik
  # loglik[[3]] <- copas.sqrt.p$BUGSoutput$sims.list$loglik
  loglik[[3]] <- copas.log.p$BUGSoutput$sims.list$loglik
  
  r_eff <- lapply(loglik, function(x){
    relative_eff(exp(x), chain_id = rep(1:3, each = 1250)) # 3 chains, each with 1250 draws
  })
  
  loo_list <- lapply(1:length(loglik), function(j){
    loo(loglik[[j]], r_eff = r_eff[[j]], k_threshold = 0.7)
  })
  
  weights[j,] <- loo_model_weights(loo_list, method = 'stacking', r_eff_list = r_eff)
  
  num.sims <- round(weights[j,] * 3750)
  sims <- list()
  
  sims[[1]] <- sample(copas.mav$BUGSoutput$sims.list$theta0, unlist(num.sims)[1], replace = FALSE)
  sims[[2]] <- sample(copas.bai$BUGSoutput$sims.list$theta0, unlist(num.sims)[2], replace = FALSE)
  # sims[[3]] <- sample(copas.sqrt.p$BUGSoutput$sims.list$theta0, num.sims[3], replace = FALSE)
  sims[[3]] <- sample(copas.log.p$BUGSoutput$sims.list$theta0, unlist(num.sims)[3], replace = FALSE)
  sims <- unlist(sims)
  # summary
  stacked.summary[j, ] <- summary.func(sims)
  
  if(j > 100){
    rm(sims, loo_list, r_eff, loglik)
    gc()
  }
}

params <- list(theta0, tau, gamma0.true, gamma1.true, rho)
names(params) <- c("theta0", "tau", "gamma0", "gamma1", "rho")
copas.sim.full <- list(stacked.summary, model.sums, weights, params)
names(copas.sim.full) <- c("stacked", "models", "weights", "params")
saveRDS(copas.sim.full, here("R", "Results", "copas.sim.full.rds"))

weights <- matrix(nrow = M, ncol = 3)
stacked.summary <- matrix(nrow = M, ncol = 4)
each.model.mean <- matrix(nrow = M, ncol = 3)
data.summary <- matrix(nrow = M, ncol = 5)
funnel.plots <- list()
saveRDS()



