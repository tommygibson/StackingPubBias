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
# data from meta-analysis with 37 studies
# "dataset11" from R meta-analysis textbook

bad_ex <- read.csv(here("R", "dataset11.csv"))

m1 <- metabin(Ee, Ne, Ec, Nc, data = bad_ex, sm = "OR")

copas.standard <- copas(m1, gamma0.range = c(-1.5, 1.5), gamma1.range = c(0, 1))

y <- copas.standard$TE
s <- copas.standard$seTE
S <- length(y)
gamma0.list <- seq(-1, 1, .1)
gamma1.list <- seq(0, 1, length.out = length(gamma0.list))
gammas <- expand.grid(gamma0.list, gamma1.list)

gamma0 <- gammas[,1]
gamma1 <- gammas[,2]
K <- length(gamma0)

log_lik_list_stan <- list()
mean.summary.stan <- matrix(nrow = length(gamma0), ncol = 2)
for(k in 1:K){
  copas.dat <- list(S = S,
                    y = y,
                    s = s,
                    gamma0 = gamma0[k],
                    gamma1 = gamma1[k])
  copas.fit <- stan(here("R", "copas_adj.stan"), data = copas.dat,
                    iter = 10000, thin = 2)
  
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
                    n.chains = 4, n.iter = 10000, n.thin = 2, DIC = FALSE)
  
  mean.summary.jags[k,] <- copas.fit$BUGSoutput$summary[(S + 1), c(1, 2)]
  
  log_lik_list_jags[[k]] <- as.matrix(copas.fit$BUGSoutput$sims.list$loglik)
}


r_eff_stan <- lapply(log_lik_list_stan, function(x) {
  relative_eff(exp(x), chain_id = rep(1:4, each = 2500))
})
r_eff_jags <- lapply(log_lik_list_jags, function(x) {
  relative_eff(exp(x), chain_id = rep(1:4, each = 2500))
})

### calculate loo and then model weights
stan_loo <- lapply(1:length(log_lik_list_stan), function(j) {
  loo(log_lik_list_stan[[j]], r_eff = r_eff_stan[[j]])
})

stan_jags <- lapply(1:length(log_lik_list_jags), function(j) {
  loo(log_lik_list_jags[[j]], r_eff = r_eff_jags[[j]])
})

copas_weights_stan <- loo_model_weights(log_lik_list_stan, method = "stacking", r_eff_list = r_eff_stan)
stan_weights_bma <- loo_model_weights(log_lik_list_stan, method = "pseudobma", r_eff_list = r_eff_stan)
copas_weights_jags <- loo_model_weights(log_lik_list_jags, method = "stacking", r_eff_list = r_eff_jags)
jags_weights_bma <- loo_model_weights(log_lik_list_jags, method = "pseudobma", r_eff_list = r_eff_stan)

post_stacking_data_stan <- cbind.data.frame(gammas, as.vector(copas_weights_stan), mean.summary.stan[,1])
post_stacking_data_jags <- cbind.data.frame(gammas, as.vector(copas_weights_jags), mean.summary.jags[,1])

names(post_stacking_data_stan) <- names(post_stacking_data_jags) <- c("gamma0", "gamma1", "weight", "mean")

post_stacking_data_stan %>%
  ggplot(aes(x = gamma0, y = gamma1, size = weight, color = mean)) +
  geom_point() +
  scale_color_gradient(low = "green", high = "red") +
  scale_size_continuous(range = c(0, 16)) +
  labs(title = "Stacking weight and estimated mean by choice of gamma0 and gamma1") +
  theme_bw()

post_stacking_data_jags %>%
  ggplot(aes(x = gamma0, y = gamma1, size = weight, color = mean)) +
  geom_point() +
  scale_color_gradient(low = "green", high = "red") +
  scale_size_continuous(range = c(0, 16)) +
  labs(title = "Stacking weight and estimated mean by choice of gamma0 and gamma1") +
  theme_bw()

weighted.mean(mean.summary.stan[,1], copas_weights_stan)
weighted.mean(mean.summary.jags[,1], copas_weights_jags)

plot(copas_weights_stan ~ copas_weights_jags)
