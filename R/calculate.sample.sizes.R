###### Estimating the number of initial studies necessary to end up with the 
###### right number of target studies for each selection function

propensity.func.1 <- function(p){
  p.prop <- ifelse(p < 0.005, 1, 
                   ifelse(p < 0.2, exp(-2 * p), 
                          ifelse(p < 0.5, exp(-4 * p), .1)))
  
  return(p.prop)
}
propensity.func.2 <- function(p){
  p.prop <- ifelse(p < 0.005, 1,
                   ifelse(p < 0.2, exp(-0.5 * p),
                          ifelse(p < 0.5, exp(-1 * p), .5)))
  return(p.prop)
}

# simulation originally run on 12/14/2021
set.seed(4321)

theta0 <- c(.5, .1)
tau <- 0.2


S.big <- 10000

s1 <- runif(S.big, 0.1, 0.8)
s2 <- runif(S.big, 0.1, 0.8)
theta1 <- rnorm(S.big, theta0[1], tau)
theta2 <- rnorm(S.big, theta0[2], tau)
y1 <- vector(length = S.big)
y2 <- vector(length = S.big)
p1 <- p2 <- vector(length = S.big)
select1.1 <- select1.2 <- select2.1 <- select2.2 <- vector(length = S.big)

for(i in 1:S.big){
  y1[i] <- rnorm(1, theta1[i], s1[i])
  y2[i] <- rnorm(1, theta2[i], s2[i])
  
  p1[i] <- 1 - pnorm(y1[i] / s1[i])
  p2[i] <- 1 - pnorm(y2[i] / s2[i])
  
}

prop1.1 <- propensity.func.1(p1)
prop1.2 <- propensity.func.2(p1)
prop2.1 <- propensity.func.1(p2)
prop2.2 <- propensity.func.2(p2)

for(i in 1:S.big){
  select1.1[i] <- rbinom(1, 1, prop1.1[i])
  select1.2[i] <- rbinom(1, 1, prop1.2[i])
  select2.1[i] <- rbinom(1, 1, prop2.1[i])
  select2.2[i] <- rbinom(1, 1, prop2.2[i])
}

props <- apply(cbind(select1.1, select1.2, select2.1, select2.2), 2, mean)

targets <- matrix(c(10, 20, 40, 80), nrow = 1)
sim_1_4_initials <- round(apply(targets, 2, function(x) x / props))

saveRDS(sim_1_4_initials, file = here("R", "sim.1.4.initial.S.rds"))



