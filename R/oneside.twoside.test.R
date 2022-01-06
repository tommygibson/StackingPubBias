##### testing three-step stan function for publication bias

library(rstan)
library(tidyverse)
library(here)
library(loo)
library(rjags)
library(R2jags)

options(mc.cores = 4)
# simulate some data that matches the assumed pattern
set.seed(421)
S <- 400
theta0 <- 0.2
tau <- 0.01

s <- runif(S, 0.1, 0.8)
theta <- rnorm(S, theta0, tau)
y <- vector(length = S)
p1 <- vector(length = S)
p2 <- vector(length = S)
select1 <- vector(length = S)
select2 <- vector(length = S)

for(i in 1:S){
  y[i] <- rnorm(1, theta[i], sd = s[i])
  p1[i] <- (1 - pnorm(y[i] / s[i]))
  p2[i] <- 2 * (1 - pnorm(abs(y[i]) / s[i]))
  select1[i] <- ifelse(p1[i] < 0.025, 1,
                      #ifelse(p[i] < 0.05, rbinom(1, 1, 0.6), 
                      ifelse(p1[i] < 0.05, rbinom(1, 1, 0.6), rbinom(1, 1, 0.3)))
  select2[i] <- ifelse(p2[i] < 0.05, 1,
                       ifelse(p2[i] < 0.1, rbinom(1, 1, 0.6), rbinom(1, 1, 0.3)))
}

dat <- data.frame(y = y,
                  s = s,
                  p1 = p1,
                  p2 = p2,
                  select1 = select1,
                  select2 = select2)

dat.select1 <- dat %>%
  filter(select1 == 1) %>%
  select(!c(select1, select2))
dat.select2 <- dat %>%
  filter(select2 == 1) %>%
  select(!c(select1, select2))


steps.oneside <- list(c(.025), c(.05, .025), c(.1, .05, .025))
steps.twoside <- list(c(.05), c(.1, .05), c(.2, .1, .05))

steps2.twoside.dat <- list(y = dat.select2$y,
                           s = dat.select2$s,
                           S = dim(dat.select2)[1],
                           steps = steps.twoside[[2]],
                           M = length(steps.twoside[[2]]))

steps2.oneside.dat <- list(y = dat.select1$y,
                           s = dat.select1$s,
                           S = dim(dat.select1)[1],
                           steps = steps.oneside[[2]],
                           M = length(steps.oneside[[2]]))
twoside.wrong.dat <- list(y = dat.select1$y,
                          s = dat.select1$s,
                          S = dim(dat.select1)[1],
                          steps = steps.twoside[[2]],
                          M = length(steps.twoside[[2]]))
oneside.wrong.dat <- list(y = dat.select2$y,
                          s = dat.select2$s,
                          S = dim(dat.select2)[1],
                          steps = steps.oneside[[2]],
                          M = length(steps.oneside[[2]]))
steps2.twoside <- stan(file = here("R", "Models", "step.twoside.stan"), data = steps2.twoside.dat,
                       iter = 2000, chains = 4)
steps2.oneside <- stan(file = here("R", "Models", "step.oneside.stan"), data = steps2.oneside.dat,
                       iter = 2000, chains = 4)
twoside.wrong <- stan(file = here("R", "Models", "step.twoside.stan"), data = twoside.wrong.dat,
                       iter = 2000, chains = 4)
oneside.wrong <- stan(file = here("R", "Models", "step.oneside.stan"), data = oneside.wrong.dat,
                       iter = 2000, chains = 4)
summary(steps2.twoside, pars = c("theta", "tau", "omega"))$summary
summary(steps2.oneside, pars = c("theta", "tau", "omega"))$summary
summary(twoside.wrong, pars = c("theta", "tau", "omega"))$summary
summary(oneside.wrong, pars = c("theta", "tau", "omega"))$summary

loglik_twoside <- list()
loglik_twoside[[1]] <- extract_log_lik(steps2.twoside, parameter_name = "loglik")
loglik_twoside[[2]] <- extract_log_lik(oneside.wrong, parameter_name = "loglik")
loglik_oneside <- list()
loglik_oneside[[1]] <- extract_log_lik(steps2.oneside, parameter_name = "loglik")
loglik_oneside[[2]] <- extract_log_lik(twoside.wrong, parameter_name = "loglik")

r_eff_twoside <- lapply(loglik_twoside, function(x){
  relative_eff(exp(x), chain_id = rep(1:4, each = 1000))
})
loo_twoside <- lapply(1:length(loglik_twoside), function(j){
  loo::loo(loglik_twoside[[j]], r_eff = r_eff_twoside[[j]])
})

r_eff_oneside <- lapply(loglik_oneside, function(x){
  relative_eff(exp(x), chain_id = rep(1:4, each = 1000))
})
loo_oneside <- lapply(1:length(loglik_oneside), function(j){
  loo::loo(loglik_oneside[[j]], r_eff = r_eff_oneside[[j]])
})

r_eff_list <- lapply(loglik_list, function(x){
  relative_eff(exp(x), chain_id = rep(1:2, each = 1000))
})

step.weights <- loo_model_weights(loglik_list, method = 'stacking', r_eff_list = r_eff_list)
