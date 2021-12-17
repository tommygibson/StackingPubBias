// Copas selection model
// Prior distributions taken from Bai 2020

data {
  int<lower=0> S; // number of studies
  real y[S]; // observed effect size
  real<lower=0> s[S]; // observed standard error

}

transformed data {
  real <lower=0> s2[S];
  real s_max = max(s);
  for(i in 1:S){
    s2[i] = s[i] ^ 2;
  }
}

parameters {
  real theta; // overall mean
  real<lower=0> tau; // standard devaition of random effects
  real<lower=-1, upper=1> rho; // correlation between observed effect y and latent factor z
  real gamma0;   // lower probability of publication
  real<lower=0> gamma1;  // higher probability of publication
  vector<lower=0>[S] z;  // latent factor determining publication probability
}

transformed parameters {
  real <lower=0> tau2 = tau^2; // variance of random effects
  vector[S] mean_y;
  vector<lower=0>[S] sd_y;
  vector<lower=0>[S] tot_var_y;
  vector[S] u;
  vector[S] v;
  vector <lower=-1, upper=1>[S] rho_bar;
  
  for(i in 1:S){
    u[i] = gamma0 + gamma1 / s[i];
    rho_bar[i] = rho * s[i] / sqrt(tau2 + s2[i]);
    v[i] = (u[i] + rho_bar[i] * (y[i] - theta) / sqrt(tau2 + s2[i])) / sqrt(1 - rho_bar[i] ^ 2);
    
    mean_y[i] = theta + rho * s[i] * (z[i] - u[i]);
    sd_y[i] = sqrt(tau2 + s2[i] * (1 - (rho ^ 2)));
    tot_var_y[i] = tau2 + s2[i];
  }
}

model {

  theta ~ normal(0, 10); // hyperparameters
  tau ~ normal(0, 10)T[0,]; // heterogeneity paramter
  rho ~ uniform(-1, 1); // correlation 
  
  gamma0 ~ uniform(-2, 2); // parameters determining mean u of latent factor z
  gamma1 ~ uniform(0, s_max);
  
  for(i in 1:S){
    y[i] ~ normal(mean_y[i], sd_y[i]); // better for gibbs but what the heck
    z[i] ~ normal(u[i], 1) T[0,]; 
  }
  
}

generated quantities {
  vector[S] loglik;
  for(i in 1:S){
    loglik[i] = normal_lcdf(v[i] | 0, 1) + normal_lpdf(y[i] | theta, tot_var_y[i]) - normal_lcdf(u[i] | 0, 1);
  }
}

