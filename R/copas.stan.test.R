### Testing the stan copas files

library(R2jags)
library(rstan)
library(here)
library(tidyverse)
library(loo)


options(mc.cores = parallel::detectCores())

set.seed(1214)
propensity.func.1 <- function(p){
  p.prop <- ifelse(p < 0.005, 1, 
                   ifelse(p < 0.2, exp(-2 * p), 
                          ifelse(p < 0.5, exp(-4 * p), .1)))
  
  return(p.prop)
}

theta0 <- 0.4
tau <- 0.2
S <- 60

s <- runif(S, 0.2, 0.8)
theta <- rnorm(S, theta0, tau)
y <- vector(length = S)
p <- vector(length = S)
select <- vector(length = S)
for(i in 1:S){
  y[i] <- rnorm(1, theta[i], s[i])
  p[i] <- 1 - pnorm(y[i] / s[i])
  
  select[i] <- rbinom(1, 1, propensity.func.1(p[i]))
}

dat <- tibble(y, s, select)
dat.select <- dat %>%
  filter(select == 1)
S.select <- dim(dat.select)[1]
y <- dat.select$y
s <- dat.select$s

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

bai.dat <- log.p.dat <- std.dat <- list(S = S.select,
                                        y = y,
                                        s = s)

mav.dat <- list(S = S.select,
                y = y,
                s = s,
                L1 = 0.01, L2 = 0.49,
                U1 = 0.51, U2 = 0.99)
bai.params <- c("theta0", "tau", "rho", "gamma0", "gamma1")
mav.params <- c("theta0", "tau", "rho", "gamma0", "gamma1")

bai.stan.fit <- stan(file = here("R", "bai.stan"), data = bai.dat)
bai.stan.fit.2 <- stan(file = here("R", "bai.2.stan"), data = bai.dat)
bai.jags.fit <- jags(data = bai.dat, inits = init.gen.bai.adj, bai.params, model.file = here("R", "copas.jags.bai.adj.txt"),
                     n.chains = 4, n.iter = 20000, DIC = FALSE)
mav.stan.fit <- stan(file = here("R", "mav.stan"), data = mav.dat, iter = 4000)
mav.jags.fit <- jags(data = mav.dat, inits = init.gen.mav.adj, mav.params, model.file = here("R", "copas.jags.mavridis.adj.txt"),
                     n.chains = 4, n.iter = 20000, DIC = FALSE)

summary(bai.stan.fit.2, pars = c("theta", "tau", "rho", "gamma0", "gamma1"))$summary
bai.jags.fit$BUGSoutput$summary
summary(mav.stan.fit, pars = c("theta", "tau", "rho", "gamma0", "gamma1"))$summary
mav.jags.fit$BUGSoutput$summary


