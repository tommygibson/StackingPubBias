// Copas selection model adjusted with sigma2

data {
  int<lower=0> S; // number of studies
  real y[S]; // observed effect size
  real<lower=0> s[S]; // observed standard error

}

transformed data {
  real<lower=0> s2[S];
  real<lower=0> s_max = max(s);
  
  for(i in 1:S){
    s2[i] = s[i] ^ 2;
  }
}

parameters {
  // vector[S] theta; // don't measure random effects explicitly 
  real theta; // overall mean
  real<lower=0> tau; // standard devaition of random effects
  real<lower=-1, upper=1> rho; // transformed correlation to be modeled with beta distribution
  real gamma0;
  real<lower=0> gamma1;
}

transformed parameters {
  real <lower=0> tau2 = tau^2; // variance of random effects
  vector <lower=-1, upper=1>[S] rho_bar;
  vector[S] u;
  vector[S] v;
  vector[S] sd_y;
  for(i in 1:S){
    u[i] = gamma0 + gamma1 / s[i];
    rho_bar[i] = rho * sqrt(s2[i]) / sqrt(tau2 + s2[i]);
    v[i] = (u[i] + rho_bar[i] * (y[i] - theta) / sqrt(tau2 + s2[i])) / sqrt(1 - rho_bar[i] ^ 2);
    sd_y[i] = sqrt(tau2 + s2[i]);
  }
}

model {

  theta ~ normal(0, 10); // hyperparameters
  tau ~ cauchy(0, 1)T[0,]; 
  rho ~ uniform(-1, 1); // uniform prior on correlation
  gamma0 ~ uniform(-2, 2);
  gamma1 ~ uniform(0, s_max);
  
  
  
  for(i in 1:S){
    
    target += normal_lcdf(v[i] | 0, 1);
    target += normal_lpdf(y[i] | theta, sd_y[i]) ;
    target += - normal_lcdf(u[i] | 0, 1);

  }
  
}

generated quantities {
  vector[S] log_lik;
  for(i in 1:S){
    log_lik[i] = normal_lcdf(v[i] | 0, 1) +
      normal_lpdf(y[i] | theta, sd_y[i]) -
      normal_lcdf(u[i] | 0, 1);
  }
}

