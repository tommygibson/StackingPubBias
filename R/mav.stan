// Copas selection model
// Prior distributions taken from Mavridis 2013

data {
  int<lower=0> S; // number of studies
  real y[S]; // observed effect size
  real s[S]; // observed standard error

}

transformed data {
  real <lower=0> s2[S];
  real s_min = min(s);
  real s_max = max(s);
  for(i in 1:S){
    s2[i] = s[i] ^ 2;
  }
}

parameters {
  real theta; // overall mean
  real<lower=0> tau; // standard devaition of random effects
  real<lower=-1, upper=1> rho; // correlation between observed effect y and latent factor z
  real<lower=0, upper=1> p_low;   // lower probability of publication
  real<lower=0, upper=1> p_high;  // higher probability of publication
  vector<lower=0>[S] z;  // latent factor determining publication probability
}

transformed parameters {
  real <lower=0> tau2 = tau^2; // variance of random effects
  real c1 = inv_Phi(p_low);
  real c2 = inv_Phi(p_high);
  real gamma0 = c1 - (c2 - c1) / (1 / s_min + 1 / s_max) / s_max;  // transforming from upper/lower probs to gammas
  real<lower=0> gamma1 = (c2 - c1) / (1 / s_min + 1 / s_max);
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
  tau ~ normal(0, 10)T[0,];
  rho ~ uniform(-1, 1); 
  
  p_low ~ uniform(.01, .5);
  p_high ~ uniform(.5, .99);
  
  for(i in 1:S){
    y[i] ~ normal(mean_y[i], sd_y[i]);
    z[i] ~ normal(u[i], 1)T[0,];
  }
  
}

generated quantities {
  vector[S] loglik;
  for(i in 1:S){
    loglik[i] = normal_lcdf(v[i] | 0, 1) + normal_lpdf(y[i] | theta, tot_var_y[i]) - normal_lcdf(u[i] | 0, 1);
  }
}

