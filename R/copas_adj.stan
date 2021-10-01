// Copas selection model adjusted with sigma2

data {
  int<lower=0> S; // number of studies
  real y[S]; // observed effect size
  real s[S]; // observed standard error
  
  real gamma0; 
  real gamma1;
}

transformed data {
  real <lower=0> s2[S];
  real u[S];
  real <lower=0> c2[S];
  
  for(i in 1:S){
    s2[i] = s[i] ^ 2;
    u[i] = gamma0 + gamma1 / s[i];
    c2[i] = exp(normal_lpdf(u[i] | 0, 1)) / exp(normal_lcdf(u[i] | 0, 1)) * 
    (u[i] + exp(normal_lpdf(u[i] | 0, 1)) / exp(normal_lcdf(u[i] | 0, 1)));
  }
}

parameters {
  // vector[S] mu; // don't measure random effects explicitly 
  real mu; // overall mean
  real<lower=0> tau; // standard devaition of random effects
  real<lower=0, upper=1> rho_t; // transformed correlation to be modeled with beta distribution
}

transformed parameters {
  real <lower=0> tau2 = tau^2; // variance of random effects
  real <lower=-1, upper=1> rho = rho_t * 2 - 1; // inverse transformed correlation
  vector <lower=-1, upper=1>[S] rho_bar;
  vector <lower=0>[S] sigma2;
  vector[S] v;
  for(i in 1:S){
    sigma2[i] = s2[i] / (1 - c2[i] * rho);
    rho_bar[i] = rho * sqrt(sigma2[i]) / sqrt(tau2 + sigma2[i]);
    v[i] = (u[i] + rho_bar[i] * (y[i] - mu) / sqrt(tau2 + sigma2[i])) / sqrt(1 - rho_bar[i] ^ 2);
  }
}

model {

  mu ~ normal(0, 10); // hyperparameters
  tau ~ cauchy(0, 1)T[0,];
  rho_t ~ beta(2, 2); // beta prior on the transformed correlation
  
  
  for(i in 1:S){
    
    target += normal_lcdf(v[i] | 0, 1);
    // target += normal_lpdf(y[i] | mu, (tau2 + s2[i])) ;
    target += -0.5 * log(tau2 + sigma2[i]) - (y[i] - mu) ^ 2 / (2 * (tau2 + sigma2[i])) - normal_lcdf(u[i] | 0, 1);

  }
  
}

generated quantities {
  vector[S] log_lik;
  for(i in 1:S){
    log_lik[i] = normal_lcdf(v[i] | 0, 1)-0.5 * log(tau2 + sigma2[i]) - 
      (y[i] - mu) ^ 2 / (2 * (tau2 + sigma2[i])) - 
      normal_lcdf(u[i] | 0, 1);
  }
}

