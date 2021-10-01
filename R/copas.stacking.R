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

copas.standard <- copas(m1)

y <- copas.standard$TE
s <- copas.standard$seTE
S <- length(y)
gamma0.list <- seq(-1.5, 1.5, 0.25)
gamma1.list <- seq(0, 1, length.out = length(gamma0.list))
gammas <- expand.grid(gamma0.list, gamma1.list)

gamma0 <- gammas[,1]
gamma1 <- gammas[,2]
K <- length(gamma0)

log_lik_list <- list()
mean.summary <- matrix(nrow = length(gamma0), ncol = 2)
for(k in 1:K){
  copas.dat <- list(S = S,
                    y = y,
                    s = s,
                    gamma0 = gamma0[k],
                    gamma1 = gamma1[k])
  copas.fit <- stan(here("R", "copas_adj.stan"), data = copas.dat,
                    iter = 5000, thin = 2)
  
  mean.summary[k,] <- summary(copas.fit, pars = "mu", probs = c(.025, .5, .975))$summary[c(1, 3)]
  
  log_lik_list[[k]] <- rstan::extract(copas.fit, pars = c("log_lik"))$log_lik
}

copas_weights <- loo_model_weights(log_lik_list, method = "stacking")

post_stacking_data <- cbind.data.frame(gammas, as.vector(copas_weights), mean.summary[,1])
names(post_stacking_data) <- c("gamma0", "gamma1", "weight", "mean")

post_stacking_data %>%
  ggplot(aes(x = gamma0, y = gamma1, size = weight, color = mean)) +
  geom_point() +
  scale_color_gradient(low = "blue", high = "pink") +
  labs(title = "Stacking weight and estimated mean by choice of gamma0 and gamma1") +
  theme_bw()
