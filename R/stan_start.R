library(rstan)
library(tidyverse)
library(here)
library(metasens)
library(R2jags)

example(stan_model, package = "rstan", run.dontrun = TRUE)

schools_dat <- list(J = 8, 
                    y = c(28,  8, -3,  7, -1,  1, 18, 12),
                    sigma = c(15, 10, 16, 11,  9, 11, 10, 18))

fit <- stan(file = here("R", "schools.stan"), data = schools_dat)

print(fit)
plot(fit)
pairs(fit, pars = c("mu", "tau", "lp__"))

la <- extract(fit, permuted = TRUE) # return a list of arrays 
mu <- la$mu 

a <- extract(fit, permuted = FALSE) 

### use S3 functions on stanfit objects
a2 <- as.array(fit)
m <- as.matrix(fit)
d <- as.data.frame(fit)



######## Let's make a toy example to try a binomial meta-analysis

S <- 10
N <- round(runif(S, 250, 1000))
beta0 <- -2
delta0 <- 2
nu0 <- -1
sigma.beta <- 0.75
sigma.delta <- 0.75
sigma.nu <- 0.75

betai <- rnorm(S, beta0, sigma.beta)
deltai <- rnorm(S, delta0, sigma.delta)
nui <- rnorm(S, nu0, sigma.nu)

pi1 <- 1 / (1 + exp(-(betai + deltai / 2)))
pi0 <- 1 / (1 + exp(-(betai - deltai / 2)))
psi <- 1 / (1 + exp(-nui))

y1 <- y0 <- n1 <- vector(length = S)

for(i in 1:S){
  
  n1[i] <- rbinom(1, N[i], psi[i])
  n0 <- N - n1
  y1[i] <- rbinom(1, n1[i], pi1[i])
  y0[i] <- rbinom(1, n0[i], pi0[i])
  
}


meta.dat <- list(n1 = n1, n0 = n0,
                 y1 = y1, y0 = y0,
                 S = S, N = N, L = 100)

meta.fit <- stan(file = here("R", "simple_meta.stan"), data = meta.dat)

meta.summary <- summary(meta.fit, pars = c("beta0", "delta0", "nu0", "PPVmean", "sensmean"))
print(meta.summary$summary)

##### copas model!
# set.seed(1)
# S <- 20
# mu0 <- 0.5
# tau <- 0.25
# s <- runif(S, 0.2, 0.6)
# gamma0 <- -0.5
# gamma1 <- 0.2
# rho <- 0.6
# 
# mu <- rnorm(S, mu0, tau)
# 
# Sigma <- list()
# y.z <- as.data.frame(matrix(nrow = S, ncol = 2))
# 
# for(i in 1:S){
#   Sigma[[i]] <- matrix(c(s[i] ^ 2, rho * s[i], rho * s[i], 1), nrow = 2)
#   y.z[i,] <- MASS::mvrnorm(1, c(mu[i], gamma0 + gamma1 / s[i]), Sigma[[i]])
# }
# 
# y.z.s <- cbind.data.frame(y.z, s)
# colnames(y.z.s) <- c("y", "z", "s")
# 
# sum(y.z.s$z > 0)
# 
# # truncate at z = 0
# y.z.s <- y.z.s %>%
#   filter(z > 0)
# 
# copas.dat <- list(S = dim(y.z.s)[1], 
#                   y = y.z.s$y, 
#                   s = y.z.s$s,
#                   gamma0 = -0.5,
#                   gamma1 = 0.4)
# copas.dat.wrong <- list(S = dim(y.z.s)[1], 
#                         y = y.z.s$y, 
#                         s = y.z.s$s,
#                         gamma0 = 3,
#                         gamma1 = 1)
# 
# copas.fit.wrong <- stan(file = here("R", "copas.stan"), data = copas.dat,
#                         iter = 15000)
# copas.wrong.summary <- summary(copas.fit.wrong)
# copas.wrong.summary$summary
# copas.summary$summary



data("Fleiss1993cont", package = "meta")

m <- metacont(n.psyc, mean.psyc, sd.psyc, n.cont, mean.cont, sd.cont, data = Fleiss1993cont)

copas.standard <- copas(m)

Fleiss.stan.dat <- list(y = m$TE,
                        s = m$seTE,
                        S = 5,
                        gamma0 = 1,
                        gamma1 = .6)

copas.fit <- stan(file = here("R", "copas.stan"), data = Fleiss.stan.dat,
                  iter = 50000, thin = 3)
copas.fit2 <- stan(file = here("R", "copas2.stan"), data = Fleiss.stan.dat,
                   iter = 100000, thin = 1)
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





