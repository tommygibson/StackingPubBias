// Recreating the model in Maier 2021 RoBMA
// Weighted normal distribution is the key, which we need to define
// three-step selection function, with different publication probabilities
// for p < .05, p \in (.05, .10), p > .10
functions {
  real w(vector omega, real p){
    if(p > 0.10){
      return omega[1];
    }
    else if(p > 0.05){
      return omega[2];
    }
    else{
      return omega[3];
    }
  }
  
  real weightednormal_lpdf(real y, real sigma, real p, real theta, real tau, vector omega) {
    // likelihood definition from vevea & hedges 1995
    int K = num_elements(omega); // dimension of weights (last element is 1)
    real num; // numerator 
    real denom[K]; // denominator will be sum of weighted integrals
    real denom_sum = 0; // 
    real tot_var = sigma ^ 2 + tau ^ 2; // study-specific var plus random effects var
    real tot_sd = sqrt(tot_var);
    real logdensity;
    real crit1 = (1.645 * sigma - theta) / tot_sd; // critical values associated with 
    real crit2 = (1.96 * sigma - theta) / tot_sd;  // p = .10 and p = .05
    
    denom[1] = exp(normal_lcdf(crit1 | 0, 1)); 
    denom[2] = exp(normal_lcdf(crit2 | 0, 1)) - denom[1];
    denom[3] = 1 - exp(normal_lcdf(crit2 | 0, 1));
    
    for(i in 1:K){
      denom_sum = denom_sum + exp(log(denom[i]) + log(omega[i]));
    }
    num = exp(normal_lpdf(y | theta, tot_sd) + log(w(omega, p)));

    
    logdensity = log(num) - log(denom_sum);
    return logdensity;
  }
}

data {
  int<lower=0> S; // number of studies
  real y[S]; // observed effect size
  real<lower = 0> s[S]; // observed standard error
  // real<lower = 0, upper = 1> p[S]; // observed p-value
}

transformed data {
  real dev[S];
  real p[S];
  for(i in 1:S){
    dev[i] = fabs(y[i]) / s[i];
    p[i] = 2 * (1 - exp(normal_lcdf(dev[i] | 0, 1))); // two-sided p-values
  }
}


parameters {
  real theta;
  real<lower = 0> tau;
  simplex[3] pre_omega;
}

transformed parameters {
  vector<lower = 0, upper = 1>[3] omega = cumulative_sum(pre_omega);
}

model {

  theta ~ normal(0, 1); // prior for mean hyperparameter
  tau ~ cauchy(0, 1)T[0,]; // half-cauchy for heterogeneity sd
  // weights \omega are cumulative dirichlet
  pre_omega ~ dirichlet(rep_vector(1, 3));
  
  for(i in 1:S){
    y[i] ~ weightednormal(s[i], p[i], theta, tau, omega);
  }
}

generated quantities {
  vector[S] loglik;
  for(i in 1:S){
    loglik[i] = weightednormal_lpdf(y[i] | s[i], p[i], theta, tau, omega);
  }
}
