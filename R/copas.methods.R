######## Experimenting with different copas modeling
### using stan, jags, maximum likelihood

library(rstan)
library(tidyverse)
library(here)
library(metasens)
library(R2jags)
library(optimx)
library(nloptr)

data("Fleiss1993cont", package = "meta")

bad_ex <- read.csv(here("R", "dataset11.csv"))

# m <- metacont(n.psyc, mean.psyc, sd.psyc, n.cont, mean.cont, sd.cont, data = Fleiss1993cont)

# example with some more studies

m1 <- metabin(Ee, Ne, Ec, Nc, data = bad_ex, sm = "OR")

copas.standard <- copas(m1)
summary(copas.standard)

# Fleiss.stan.dat <- list(y = m$TE,
#                         s = m$seTE,
#                         S = 5,
#                         gamma0 = 0,
#                         gamma1 = .5)
copas.stan.dat <- list(y = m1$TE,
                        s = m1$seTE,
                        S = length(m1$TE),
                        gamma0 = -1,
                        gamma1 = 0)

# copas.fit <- stan(file = here("R", "copas.stan"), data = Fleiss.stan.dat,
#                   iter = 50000, thin = 3)
copas.fit2 <- stan(file = here("R", "copas2.stan"), data = copas.stan.dat,
                   iter = 10000, thin = 2)
copas.fit.adj <- stan(file = here("R", "copas_adj.stan"), data = copas.stan.dat,
                      iter = 10000, thin = 2)
copas.stan.summary <- summary(copas.fit.adj, pars = c("mu", "tau", "rho"), probs = c(.025, .5, .975))$summary




####### Let's try copas with jags

copas.jags.dat <- copas.stan.dat

copas.jags.init.gen <- function(){
  list(
    z = runif(S, 0, 1),
    theta = rnorm(0, 1),
    tau = runif(1, 0, 0.5),
    rho_t = runif(1, .01, .99)
    # rho_t = runif(1, -0.5, 0.5)
  )
}

copas.params <- c("theta", "tau", "rho", "loglik")

copas.jags <- jags(data = copas.jags.dat,
                   inits = copas.jags.init.gen,
                   parameters.to.save = copas.params,
                   model.file = here("R", "copas_jags.txt"),
                   n.chains = 3, n.iter = 10000, n.thin = 2, DIC = FALSE)
copas.jags.adj <- jags(data = copas.jags.dat,
                       inits = copas.jags.init.gen,
                       parameters.to.save = copas.params,
                       model.file = here("R", "copas_jags_adj.txt"),
                       n.chains = 3, n.iter = 10000, n.thin = 2, DIC = FALSE)

copas.jags.summary <- copas.jags$BUGSoutput$summary[, c(1, 2, 3, 5, 7)]
copas.jags.adj.summary <- copas.jags.adj$BUGSoutput$summary[, c(1, 2, 3, 5, 7)]

#### Using log-likelihood

LL_copas <- function(theta, gamma0, gamma1, y, s){
  
  S = length(y)
  
  mu <- theta[1]
  tau <- theta[2]
  rho <- theta[3]
  
  tau2 <- tau ^ 2
  
  
  rho_bar <- vector(length = S)
  u <- vector(length = S)
  v <- vector(length = S)
  loglik <- vector(length = S)
  
  for(i in 1:S){
    rho_bar[i] <- rho * s[i] / sqrt(tau2 + s[i]^2)
    u[i] <- gamma0 + gamma1 / s[i]
    v[i] <- (u[i] + rho_bar[i] * (y[i] - mu) / sqrt(tau2 + s[i] ^ 2)) / sqrt(1 - rho_bar[i] ^ 2)
    loglik[i] <- -1 / 2 * log(tau2 + s[i] ^ 2) - (y[i] - mu) ^ 2 / (2 * (tau2 + s[i]^2)) -
      log(pnorm(u[i])) + log(pnorm(v[i]))
  }
  
  return(-sum(loglik))
  
}

LL_copas_adj <- function(theta, gamma0, gamma1, y, s){
  
  S = length(y)
  
  mu <- theta[1]
  tau <- theta[2]
  rho <- theta[3]
  
  tau2 <- tau ^ 2
  
  c2 <- vector(length = S)
  sigma2 <- vector(length = S)
  rho_bar <- vector(length = S)
  u <- vector(length = S)
  v <- vector(length = S)
  loglik <- vector(length = S)
  
  for(i in 1:S){
    u[i] <- gamma0 + gamma1 / s[i]
    c2[i] <- dnorm(u[i]) / pnorm(u[i]) * (u[i] + dnorm(u[i]) / pnorm(u[i]))
    sigma2[i] <- s[i] ^ 2 / (1 - c2[i] * rho)
    
    rho_bar[i] <- rho * sqrt(sigma2[i]) / sqrt(tau2 + sigma2[i])
    v[i] <- (u[i] + rho_bar[i] * (y[i] - mu) / sqrt(tau2 + sigma2[i])) / sqrt(1 - rho_bar[i] ^ 2)
    loglik[i] <- -1 / 2 * log(tau2 + sigma2[i]) - (y[i] - mu) ^ 2 / (2 * (tau2 + sigma2[i])) -
      log(pnorm(u[i])) + log(pnorm(v[i]))
  }
  
  return(-sum(loglik))
  
}

lb <- c(-10, 0, -1)
ub <- c(10, 10, 1)

mu.init <- 0
tau.init <- 1
rho.init <- 0

theta.inits <- c(mu.init, tau.init, rho.init)

copas.opt <- nloptr(x0 = theta.inits, 
                    eval_f = LL_copas, 
                    lb = lb, ub = ub,
                    opts = list("algorithm" = "NLOPT_LN_COBYLA",
                                "xtol_rel" = 1e-8),
                    gamma0 = 0,
                    gamma1 = 0.08,
                    y = m1$TE,
                    s = m1$seTE)

copas.opt.adj <- nloptr(x0 = theta.inits, 
                        eval_f = LL_copas_adj, 
                        lb = lb, ub = ub,
                        opts = list("algorithm" = "NLOPT_LN_COBYLA",
                                    "xtol_rel" = 1e-8),
                        gamma0 = -1,
                        gamma1 = 0,
                        y = m1$TE,
                        s = m1$seTE)


# copas.opt$solution
# print(copas.stan.summary)
# print(copas.jags.summary)

copas.opt.adj$solution
print(copas.stan.summary)
# print(copas.jags.summary)


