######## Experimenting with different copas modeling
### using stan, jags, maximum likelihood

library(rstan)
library(tidyverse)
library(here)
library(metasens)
library(R2jags)
library(optimx)

data("Fleiss1993cont", package = "meta")

m <- metacont(n.psyc, mean.psyc, sd.psyc, n.cont, mean.cont, sd.cont, data = Fleiss1993cont)

copas.standard <- copas(m)

Fleiss.stan.dat <- list(y = m$TE,
                        s = m$seTE,
                        S = 5,
                        gamma0 = 0,
                        gamma1 = .5)

copas.fit <- stan(file = here("R", "copas.stan"), data = Fleiss.stan.dat,
                  iter = 50000, thin = 3)
copas.fit2 <- stan(file = here("R", "copas2.stan"), data = Fleiss.stan.dat,
                   iter = 10000, thin = 2)
copas.summary <- summary(copas.fit2)
copas.summary$summary



####### Let's try copas with jags
S <- length(m$TE)
copas.jags.dat <- list(y = m$TE, s = m$seTE, S = length(m$TE),
                       gamma0 = -.12, gamma1 = .15)

copas.jags.init.gen <- function(){
  list(
    z = runif(S, 0, 1),
    theta = rnorm(0, 1),
    tau = runif(1, 0, 0.5),
    rho = runif(1, -.99, .99)
    # rho_t = runif(1, -0.5, 0.5)
  )
}

copas.params <- c("theta", "tau", "rho")

copas.jags <- jags(data = copas.jags.dat,
                   inits = copas.jags.init.gen,
                   parameters.to.save = copas.params,
                   model.file = here("R", "copas_jags.txt"),
                   n.chains = 3, n.iter = 2000, n.thin = 2)

#### Using log-likelihood

LL_copas <- function(theta, gamma0, gamma1, y, s){
  S = length(y)
  mu <- theta[1]
  tau <- exp(theta[2])
  rho <- tanh(theta[3])
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

mu.init <- 0
tau.init <- 0
rho.init <- atanh(0.1)

theta.inits <- c(mu.init, tau.init, rho.init)

copas.opt <- optimx(theta.inits, LL_copas, 
                    gamma0 = -.65, gamma1 = 0.4,
                    y = copas.standard$TE, s = copas.standard$seTE,
                    control = list(maxit = 500))
