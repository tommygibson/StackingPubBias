##### testing three-step stan function for publication bias

library(rstan)
library(tidyverse)
library(here)
library(loo)
library(rjags)
library(R2jags)


# simulate some data that matches the assumed pattern
set.seed(421)
S <- 200
theta0 <- 0.4
tau <- 0.2

s <- runif(S, 0.2, 0.8)
theta <- rnorm(S, theta0, tau)
y <- vector(length = S)
p <- vector(length = S)
select <- vector(length = S)


for(i in 1:S){
  y[i] <- rnorm(1, theta[i], sd = s[i])
  p[i] <- min(2 * (1 - pnorm(y[i] / s[i])), 1)
  select[i] <- ifelse(p[i] < 0.05, 1,
                      ifelse(p[i] < 0.1, rbinom(1, 1, 0.6), rbinom(1, 1, 0.5)))
}

dat <- data.frame(y = y,
                  s = s,
                  p = p,
                  select = select)

dat.select <- dat %>%
  filter(select == 1)

threestep.dat <- list(y = dat.select$y,
                      s = dat.select$s,
                      S = dim(dat.select)[1])

threestep.fit <- stan(file = here("R", "three_step.stan"), data = threestep.dat,
                      iter = 1000, chains = 4)
twostep.fit <- stan(file = here("R", "two_step.stan"), data = threestep.dat,
                    iter = 1000, chains = 2)



## success!

summary(threestep.fit, pars = c("theta", "tau", "omega"))$summary
threestep.loglik <- extract_log_lik(threestep.fit, parameter_name = "loglik")

# let's stack it with some other models

mav.dat <- list(y = dat.select$y,
                s = dat.select$s,
                S = dim(dat.select)[1],
                L1 = .01, L2 = .50,
                U1 = .50, U2 = 1)

mav.inits <- function(){
  list(
    z = runif(dim(dat.select)[1], 0, 1),
    theta0 = rnorm(1, 0, 0.25),
    rho = runif(1, -0.5, 0.5),
    tau = runif(1, 0.1, 0.5),
    p.low = runif(1, 0.1, .4),
    p.high = runif(1, 0.6, 0.9)
  )
}
mav.params <- c("theta0", "tau")

mav.fit <- jags(mav.dat, mav.inits, mav.params, model.file = here("R", "copas.jags.mavridis.adj.txt"),
                n.iter = 1000, n.chains = 2)
mav.fit$BUGSoutput$summary
loglik_list <- list()
loglik_list[[1]] <- threestep.loglik
loglik_list[[2]] <- mav.fit$BUGSoutput$sims.list$loglik

r_eff <- lapply(loglik_list, function(x){
  relative_eff(exp(x), chain_id = rep(1:2, each = 500))
})

w <- loo_model_weights(loglik_list, method = "stacking", r_eff_list = r_eff)

steps <- list(c(.05), c(.1, .05), c(.2, .1, .05))

steps1.dat <- list(y = dat.select$y,
                   s = dat.select$s,
                   S = dim(dat.select)[1],
                   steps = steps[[1]],
                   M = length(steps[[1]]))
steps2.dat <- list(y = dat.select$y,
                   s = dat.select$s,
                   S = dim(dat.select)[1],
                   steps = steps[[2]],
                   M = length(steps[[2]]))
steps3.dat <- list(y = dat.select$y,
                   s = dat.select$s,
                   S = dim(dat.select)[1],
                   steps = steps[[3]],
                   M = length(steps[[3]]))

steps1.fit <- stan(file = here("R", "stepfunc.stan"), data = steps1.dat,
                   iter = 1000, chains = 4)
steps2.fit <- stan(file = here("R", "stepfunc.stan"), data = steps2.dat,
                   iter = 1000, chains = 4)
steps3.fit <- stan(file = here("R", "stepfunc.stan"), data = steps3.dat,
                   iter = 1000, chains = 4)
summary(steps3.fit, pars = c("theta", "tau", "omega"))$summary


