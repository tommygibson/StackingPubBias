###### Stacking with copas model??

library(rstan)
library(tidyverse)
library(here)
library(metasens)
library(R2jags)
library(optimx)
library(nloptr)
library(loo)
library(scico)
library(RobustBayesianCopas)
# data from meta-analysis with 37 studies
# "dataset11" from R meta-analysis textbook

row.match <- function(x, y){
  one.index <- function(a) {
    which(apply(y, 1, function(z) all(z == a)))
  }
  apply(x, 1, one.index)
}


bad_ex <- read.csv(here("R", "dataset11.csv"))

m1 <- metabin(Ee, Ne, Ec, Nc, data = bad_ex, sm = "OR")

copas.standard <- copas(m1, gamma0.range = c(-1.5, 1.5), gamma1.range = c(0, 1))

y <- copas.standard$TE
s <- copas.standard$seTE
S <- length(y)
gammas <- expand.grid(seq(-1, 1, .125), seq(0, 1, length = length(seq(-1, 1, .125))))
gamma0 <- gammas[,1]
gamma1 <- gammas[,2]
K <- length(gamma0)

gammas.1 <- expand.grid(seq(-1, 1, .25), seq(0, 1, length = length(seq(-1, 1, .25))))
gammas.2 <- expand.grid(seq(-1, 1, .5), seq(0, 1, length = length(seq(-1, 1, .5))))

# indices for subsets of gamma0 and gamma1
gammas.1.index <- row.match(gammas.1, gammas)
gammas.2.index <- row.match(gammas.2, gammas)

log_lik_list_stan <- list()
mean.summary.stan <- matrix(nrow = length(gamma0), ncol = 2)

for(k in 1:K){
  copas.dat <- list(S = S,
                    y = y,
                    s = s,
                    gamma0 = gamma0[k],
                    gamma1 = gamma1[k])
  copas.fit <- stan(here("R", "copas_adj.stan"), data = copas.dat,
                    iter = 5000, thin = 2)
  
  mean.summary.stan[k,] <- summary(copas.fit, pars = "mu", probs = c(.025, .5, .975))$summary[c(1, 3)]
  
  log_lik_list_stan[[k]] <- extract_log_lik(copas.fit)
}


log_lik_list_jags <- list()
mean.summary.jags <- matrix(nrow = length(gamma0), ncol = 2)

copas.jags.init <- function(){
  list(
    z = runif(S, 0, 1),
    theta = rnorm(0, 1),
    tau = runif(1, 0, 0.5),
    rho_t = runif(1, .01, .99)
  )
}

copas.params.jags <- c("theta", "loglik")

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
  
  mean.summary.jags[k,] <- copas.fit$BUGSoutput$summary[(S + 1), c(1, 2)]
  
  log_lik_list_jags[[k]] <- as.matrix(copas.fit$BUGSoutput$sims.list$loglik)
}


copas_results <- list(log_lik_list_stan, mean.summary.stan, 
                      log_lik_list_jags, mean.summary.jags)
names(copas_results) <- c("stan_loglik_list", "stan_means",
                          "jags_loglik_list", "jags_means")
saveRDS(copas_results, file = here("R", "Results", "copas_results.rds"))

# r_eff for stan and jags
r_eff_stan_0 <- lapply(copas_results$stan_loglik_list, function(x) {
  relative_eff(exp(x), chain_id = rep(1:4, each = 1250))
})
r_eff_stan_1 <- lapply(copas_results$stan_loglik_list[gammas.1.index], function(x) {
  relative_eff(exp(x), chain_id = rep(1:4, each = 1250))
})
r_eff_stan_2 <- lapply(copas_results$stan_loglik_list[gammas.2.index], function(x) {
  relative_eff(exp(x), chain_id = rep(1:4, each = 1250))
})

r_eff_jags_0 <- lapply(copas_results$jags_loglik_list, function(x) {
  relative_eff(exp(x), chain_id = rep(1:4, each = 1250))
})
r_eff_jags_1 <- lapply(copas_results$jags_loglik_list[gammas.1.index], function(x) {
  relative_eff(exp(x), chain_id = rep(1:4, each = 1250))
})
r_eff_jags_2 <- lapply(copas_results$jags_loglik_list[gammas.2.index], function(x) {
  relative_eff(exp(x), chain_id = rep(1:4, each = 1250))
})

# stacking weights
stan_weights_0 <- loo_model_weights(copas_results$stan_loglik_list, "stacking", r_eff_list = r_eff_stan_0)
stan_weights_1 <- loo_model_weights(copas_results$stan_loglik_list[gammas.1.index], "stacking", r_eff_list = r_eff_stan_1)
stan_weights_2 <- loo_model_weights(copas_results$stan_loglik_list[gammas.2.index], "stacking", r_eff_list = r_eff_stan_2)

jags_weights_0 <- loo_model_weights(copas_results$jags_loglik_list, "stacking", r_eff_list = r_eff_jags_0)
jags_weights_1 <- loo_model_weights(copas_results$jags_loglik_list[gammas.1.index], "stacking", r_eff_list = r_eff_jags_1)
jags_weights_2 <- loo_model_weights(copas_results$jags_loglik_list[gammas.2.index], "stacking", r_eff_list = r_eff_jags_2)

#### take subsets posterior draws to see if results are consistent

stan_loglik_sub1 <- list()
stan_loglik_sub2 <- list()
jags_loglik_sub1 <- list()
jags_loglik_sub2 <- list()


for(i in 1:length(copas_results$stan_loglik_list)){
  stan_loglik_sub1[[i]] <- copas_results$stan_loglik_list[[i]][c(1:625, 1251:1875, 2501:3125, 3751:4375),]
  stan_loglik_sub2[[i]] <- copas_results$stan_loglik_list[[i]][c(626:1250, 1876:2500, 3126:3750, 4376:5000),]
  jags_loglik_sub1[[i]] <- copas_results$jags_loglik_list[[i]][c(1:625, 1251:1875, 2501:3125, 3751:4375),]
  jags_loglik_sub2[[i]] <- copas_results$jags_loglik_list[[i]][c(626:1250, 1876:2500, 3126:3750, 4376:5000),]
}

r_eff_stan_0_sub1 <- lapply(stan_loglik_sub1, function(x) {
  relative_eff(exp(x), chain_id = rep(1:4, each = 625))
})
r_eff_stan_0_sub2 <- lapply(stan_loglik_sub2, function(x) {
  relative_eff(exp(x), chain_id = rep(1:4, each = 625))
})

r_eff_jags_0_sub1 <- lapply(jags_loglik_sub1, function(x) {
  relative_eff(exp(x), chain_id = rep(1:4, each = 625))
})
r_eff_jags_0_sub2 <- lapply(jags_loglik_sub2, function(x) {
  relative_eff(exp(x), chain_id = rep(1:4, each = 625))
})

stan_weights_0_sub1 <- loo_model_weights(stan_loglik_sub1, "stacking", r_eff_list = r_eff_stan_0_sub1)
stan_weights_0_sub2 <- loo_model_weights(stan_loglik_sub2, "stacking", r_eff_list = r_eff_stan_0_sub2)
jags_weights_0_sub1 <- loo_model_weights(jags_loglik_sub1, "stacking", r_eff_list = r_eff_jags_0_sub1)
jags_weights_0_sub2 <- loo_model_weights(jags_loglik_sub2, "stacking", r_eff_list = r_eff_jags_0_sub2)

stacking_0 <- cbind.data.frame(gammas, as.vector(stan_weights_0), as.vector(jags_weights_0), 
                               as.vector(stan_weights_0_sub1), as.vector(stan_weights_0_sub2),
                               as.vector(jags_weights_0_sub1), as.vector(jags_weights_0_sub2),
                               mean.summary.stan[,1], mean.summary.jags[,1])

stacking_1 <- cbind.data.frame(gammas.1, as.vector(stan_weights_1), as.vector(jags_weights_1),
                               mean.summary.stan[gammas.1.index, 1], mean.summary.jags[gammas.1.index, 1])

stacking_2 <- cbind.data.frame(gammas.2, as.vector(stan_weights_2), as.vector(jags_weights_2),
                               mean.summary.stan[gammas.2.index, 1], mean.summary.jags[gammas.2.index, 1])
names(stacking_0) <- c("gamma0", "gamma1", "StanWeight", "JAGSWeight", 
                       "Stan_subset1", "Stan_subset2", "JAGS_subset1", "JAGS_subset2",
                       "StanMean", "JAGSMean")

names(stacking_1) <- names(stacking_2) <-
  c("gamma0", "gamma1", "StanWeight", "JAGSWeight", "StanMean", "JAGSMean")

copas.stacking <- list(stacking_0, stacking_1, stacking_2)
names(copas.stacking) <- c("stacking.0", "stacking.1", "stacking.2")
saveRDS(copas.stacking, here("R", "Results", "copas.stacking.rds"))

### calculate loo and then model weights
# stan_loo <- lapply(1:length(log_lik_list_stan), function(j) {
#   loo(log_lik_list_stan[[j]], r_eff = r_eff_stan[[j]])
# })
# 
# stan_jags <- lapply(1:length(log_lik_list_jags), function(j) {
#   loo(log_lik_list_jags[[j]], r_eff = r_eff_jags[[j]])
# })


