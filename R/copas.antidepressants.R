#### Comparing results from antidepressants dataset

library(RobustBayesianCopas)
library(tidyverse)
library(R2jags)
library(here)
library(loo)
library(nloptr)

anti.select <- filter(antidepressants,
                      !is.na(Standardized_effect_size))

y <- anti.select$Standardized_effect_size
s <- anti.select$Standardized_SE
S <- length(y)


### Data for each model

bai.dat <- list(S = S,
                y = y,
                s = s)
mav.dat <- list(S = S,
                y = y,
                s = s,
                L1 = 0.35, L2 = 0.45,
                U1 = 0.75, U2 = 0.85)
mav.dat.diff <- list(S = S,
                     y = y,
                     s = s,
                     L1 = 0.3, L2 = 0.99,
                     U2 = 1)
init.gen.bai <- function(){
  list(
    z = runif(S, 0, 1),
    theta = rnorm(S, 0, 0.25),
    theta0 = rnorm(1, 0, 0.25),
    rho = runif(1, -0.5, 0.5),
    tau = runif(1, 0.1, 0.5),
    gamma0 = runif(1, -1, 1),
    gamma1 = runif(1, 0, max(s))
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
init.gen.mav <- function(){
  list(
    z = runif(S, 0, 1),
    theta = rnorm(S, 0, 0.25),
    theta0 = rnorm(1, 0, 0.25),
    rho = runif(1, -0.5, 0.5),
    tau = runif(1, 0.1, 0.5),
    p.low = runif(1, 0.38, 0.42),
    p.high = runif(1, 0.78, 0.82)
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
init.gen.logp <- function(){
  list(
    z = runif(S, 0, 1),
    theta0 = rnorm(1, 0, 0.25),
    rho = runif(1, -0.5, 0.5),
    tau = runif(1, 0.1, 0.5),
    gamma0 = runif(1, -1, 1),
    gamma1 = runif(1, 0, 0.05)
  )
}
init.gen.intp <- function(){
  list(
    z = runif(S, 0, 1),
    theta0 = rnorm(1, 0, 0.25),
    rho = runif(1, -0.5, 0.5),
    tau = runif(1, 0.1, 0.5),
    gamma0 = runif(1, -1, 1),
    gamma1 = runif(1, -0.5, 0.5),
    gamma2 = runif(1, -0.5, 0.5)
  )
}
init.gen.std <- function(){
  list(
    theta0 = rnorm(1, 0, 0.25),
    tau = runif(1, 0.1, 0.4)
  )
}

std.params <- c("theta0", "tau")
std.dat <- list(S = S,
                y = y,
                s = s)

fit.std <- jags(data = std.dat, parameters.to.save = std.params, inits = init.gen.std,
                model.file = here("R", "std.meta.txt"), n.chains = 2, n.iter = 5000, DIC = FALSE)
mav.params <- c("theta0", "rho", "tau", "gamma0", "gamma1")
bai.params <- c("theta0", "rho", "tau", "gamma0", 'gamma1')
intp.params <- c("theta0", "rho", "tau", "gamma0", "gamma1", "gamma2")
copas.mav <- jags(data = mav.dat, inits = init.gen.mav, parameters.to.save = mav.params,
                  model.file = here("R", "copas.jags.mavridis.txt"),
                  n.iter = 10000, n.thin = 4, n.chains = 4)
copas.mav.adj <- jags(data = mav.dat, inits = init.gen.mav.adj, parameters.to.save = mav.params,
                      model.file = here("R", "copas.jags.mavridis.adj.txt"),
                      n.iter = 10000, n.thin = 4, n.chains = 4)
copas.mav.diff <- jags(data = mav.dat.diff, inits = init.gen.mav.adj, parameters.to.save = mav.params,
                       model.file = here("R", "mavridis.diff.prior.txt"),
                       n.iter = 10000, n.thin = 4, n.chains = 4)

copas.bai <- jags(data = bai.dat, inits = init.gen.bai, parameters.to.save = bai.params,
                  model.file = here("R", "copas.jags.bai.txt"),
                  n.iter = 10000, n.thin = 2, n.chains = 4)
copas.bai.adj <- jags(data = bai.dat, inits = init.gen.bai.adj, parameters.to.save = bai.params,
                      model.file = here("R", "copas.jags.bai.adj.txt"),
                      n.iter = 100000, n.thin = 2, n.chains = 4)
stack.params <- c('theta0', 'loglik')

copas.mav.stack <- jags(data = mav.dat, inits = init.gen.mav.adj, parameters.to.save = stack.params,
                        model.file = here("R", "copas.jags.mavridis.adj.txt"),
                        n.iter = 5000, n.thin = 2, n.chains = 4, DIC = FALSE)
copas.bai.stack <- jags(data = bai.dat, inits = init.gen.bai.adj, parameters.to.save = stack.params,
                        model.file = here("R", "copas.jags.bai.adj.txt"),
                        n.iter = 5000, n.thin = 2, n.chains = 4, DIC = FALSE)



mav.bai.loglik <- list()
mav.bai.loglik[[1]] <- copas.mav.stack$BUGSoutput$sims.list$loglik
mav.bai.loglik[[2]] <- copas.bai.stack$BUGSoutput$sims.list$loglik
mav.bai.loglik[[3]] <- copas.log.p$BUGSoutput$sims.list$loglik

r_eff_mav_bai <- lapply(mav.bai.loglik, function(x){
  relative_eff(exp(x), chain_id = rep(1:4, each = 1250))
})

loo_list_mav_bai <- lapply(1:length(mav.bai.loglik), function(j){
  loo(mav.bai.loglik[[j]], r_eff = r_eff_mav_bai[[j]], k_threshold = 0.7)
})

mav.bai.weights <- loo_model_weights(loo_list_mav_bai, method = 'stacking', r_eff_list = r_eff_mav_bai)


copas.sqrt.p <- jags(data = bai.dat, inits = init.gen.bai.adj, parameters.to.save = bai.params,
                     model.file = here('R', 'copas.pval.sqrt.txt'),
                     n.iter = 10000, n.thin = 4, n.chains = 4, DIC = FALSE)
copas.log.p <- jags(data = bai.dat, inits = init.gen.logp, parameters.to.save = bai.params,
                     model.file = here('R', 'copas.pval.log.txt'),
                     n.iter = 10000, n.thin = 4, n.chains = 4, DIC = FALSE)

plot(s ~ y)

gamma0.list <- seq(-1.2, -0.4, length.out = 10)
gamma1.list <- seq(0, 0.3, length.out = 10)

gammas <- expand.grid(gamma0.list, gamma1.list)

gamma0 <- gammas[,1]
gamma1 <- gammas[,2]

K <- length(gamma0)

loglik <- list()
mean.sims <- list()
mean.sum <- matrix(nrow = length(gamma0), ncol = 2)

copas.jags.init <- function(){
  list(
    z = runif(S, 0, 1),
    theta = rnorm(0, 1),
    tau = runif(1, 0, 0.5),
    rho = runif(1, -0.99, .99)
  )
}

copas.params.jags <- c("theta", "rho", "loglik")

for(k in 1:K){
  copas.dat <- list(S = S,
                    y = y,
                    s = s,
                    gamma0 = gamma0[k],
                    gamma1 = gamma1[k])
  copas.fit <- jags(data = copas.dat,
                    inits = copas.jags.init,
                    parameters.to.save = copas.params.jags,
                    model.file = here("R", "copas_jags_adj.txt"),
                    n.chains = 4, n.iter = 5000, n.thin = 2, DIC = FALSE)
  
  mean.sims[[k]] <- cbind(copas.fit$BUGSoutput$sims.list$theta, copas.fit$BUGSoutput$sims.list$rho)
  mean.sum[k,] <- copas.fit$BUGSoutput$summary[(S + 1) : (S + 2), 1]
  
  loglik[[k]] <- as.matrix(copas.fit$BUGSoutput$sims.list$loglik)
}

r_eff <- lapply(loglik, function(x){
  relative_eff(exp(x), chain_id = rep(1:4, each = 1250))
})

loo_list <- lapply(1:length(loglik), function(j) {
  loo(loglik[[j]], r_eff = r_eff[[j]], k_threshold = 0.7)
})

weights <- loo_model_weights(loo_list, method = "stacking", r_eff_list = r_eff)

num.sims <- round(weights * 5000)
sims <- list()
for(k in 1:K){
  sims[[k]] <- cbind(sample(mean.sims[[k]][,1], num.sims[k], replace = FALSE),
                     sample(mean.sims[[k]][,2], num.sims[k], replace = FALSE))
                     
}
sims <- do.call(rbind, sims)
summary.func <- function(x){
  c(mean(x), sd(x), quantile(x, c(.025, 0.5, .975)))
}

summary.func(sims[,1])
summary.func(sims[,2])


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
                    gamma0 = -0.5,
                    gamma1 = 0.08,
                    y = y,
                    s = s)

copas.opt.adj <- nloptr(x0 = theta.inits, 
                        eval_f = LL_copas_adj, 
                        lb = lb, ub = ub,
                        opts = list("algorithm" = "NLOPT_LN_COBYLA",
                                    "xtol_rel" = 1e-8),
                        gamma0 = -0.5,
                        gamma1 = 0.08,
                        y = y,
                        s = s)

