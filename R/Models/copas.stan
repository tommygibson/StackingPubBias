// Copas selection model

data {
  int<lower=0> S; // number of studies
  real y[S]; // observed effect size
  real s[S]; // observed standard error
  
  real gamma0; 
  real gamma1;
}

parameters {
  vector[S] mu; // true study-specific effect sizes
  real mu0; // mean of random effects
  real<lower=0> tau; // standard devaition of random effects
  
  vector<lower=0>[S] z; // latent variable, propensity for publication
  real<lower=0, upper=1> rho_t; // correlation between effect sizes and propensity for publication
}

transformed parameters {
  real<lower=-1, upper=1> rho = rho_t * 2 - 1;
}

model {
  // (yi, zi) is truncated bivariate normal, so we sample truncated zi and then yi given zi
  for(i in 1:S){
    z[i] ~ normal(gamma0 + gamma1 / s[i], 1)T[0,];
  
    y[i] ~ normal(mu[i] + rho * s[i] * (z[i] - gamma0 - gamma1 / s[i]), sqrt(s[i] ^ 2 * (1 - rho ^ 2)));
  }
  
  mu ~ normal(mu0, tau); // random effects
  mu0 ~ normal(0, 1); // hyperparameters
  tau ~ normal(0, 1)T[0,];
  rho_t ~ beta(2, 2); // beta prior on the transformed correlation
  
}

