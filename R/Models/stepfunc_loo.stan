// Recreating the model in Maier 2021 RoBMA
// Weighted normal distribution is the key, which we need to define
// can handle any set of steps, which are input as data
// This one can also handle loo or k-fold cross-validation
// to calculate elpd

functions {
  real w(vector omega, real p, vector steps){
    int K = num_elements(steps) + 1;
    vector[K] steps1 = append_row(1, sort_desc(steps));
    real weight = omega[K];
    
    for(k in 1:(K - 1)){
      if((p <= steps1[k]) && (p > steps1[k + 1])){
        weight = omega[k];
        break;
      }
      else continue;
    }
    return weight;
  }
  
  // real crit(real steps){
  //   return fabs(inv_Phi(steps / 2));
  // }
  
  real weightednormal_lpdf(real y, real sigma, real p, real theta, real tau, vector omega, vector steps) {
    // likelihood definition from vevea & hedges 1995
    int K = num_elements(omega); // number of weights (last element is 1)
    int T = num_elements(steps); // number of steps
    
    real num; // numerator 
    real denom = 0;
    real denoms[K]; // denominator will be sum of weighted integrals
    real denom_sum = 0;
    real tot_var = sigma ^ 2 + tau ^ 2; // study-specific var plus random effects var
    real tot_sd = sqrt(tot_var);
    real logdensity;
    vector[T] crits; // critical values

    for(t in 1:T){
      crits[t] = (fabs(inv_Phi(steps[t] / 2)) * sigma - theta) / tot_sd; // define std normal quantiles associated with steps
    }

    // first denominator integral
    denoms[1] = exp(normal_lcdf(crits[1] | 0, 1));
    denom_sum = denoms[1];
    if(K > 2) {
      for(k in 2:(K - 1)){ 
        // middle denominator integrals
        denoms[k] = exp(normal_lcdf(crits[k] | 0, 1)) - denom_sum;
        denom_sum = denom_sum + denoms[k];
      }
      denoms[K] = 1 - denom_sum;
    }
    else {
      // last denominator integral
      denoms[K] = 1 - denom_sum;
    }
    // sum the weighted integrals
    for(i in 1:K){
      denom = denom + exp(log(denoms[i]) + log(omega[i]));
    }
    num = exp(normal_lpdf(y | theta, tot_sd) + log(w(omega, p, steps)));

    
    logdensity = log(num) - log(denom);
    return logdensity;
  }
}

data {
  int<lower = 0> S; // number of studies
  int<lower = 0> H; // number of holdout studies (gonna be 1 for loo)
  int<lower = 0> M; // number of steps
  real y[S]; // observed effect size
  real<lower = 0> s[S]; // observed standard error
  real y_h[H]; // observed effects for holdouts
  real<lower = 0> s_h[H]; // observed std error for holdouts
  
  vector<lower = 0, upper = 1>[M] steps;
}

transformed data {
  real dev[S];
  real p[S];
  real dev_h[H];
  real<lower = 0, upper = 1> p_h[H];
  vector[M] steps1 = sort_desc(steps);
  for(i in 1:S){
    dev[i] = fabs(y[i]) / s[i];
    p[i] = 2 * (1 - exp(normal_lcdf(dev[i] | 0, 1))); // two-sided p-values
  }
  for(i in 1:H){
    dev_h[i] = fabs(y_h[i]) / s_h[i];
    p_h[i] = 2 * (1 - exp(normal_lcdf(dev_h[i] | 0, 1))); // two-sided p-values
  }
}

// overall mean, heterogeneity, non-cumulative-sum weights
parameters {
  real theta;
  real<lower = 0> tau;
  simplex[M + 1] pre_omega;
}

// transform from pre-omega to omega with cumulative sum
transformed parameters {
  vector<lower = 0, upper = 1>[M + 1] omega = cumulative_sum(pre_omega);
}

model {

  theta ~ normal(0, 10); // prior for mean hyperparameter
  tau ~ cauchy(0, 1)T[0,]; // half-cauchy for heterogeneity sd
  // weights \omega are cumulative dirichlet
  // auxiliary parameter pre_omega is dirichlet
  pre_omega ~ dirichlet(rep_vector(1, (M + 1)));
  
  for(i in 1:S){
    y[i] ~ weightednormal(s[i], p[i], theta, tau, omega, steps1);
  }
}

// We only calculate loglik for holdouts here
generated quantities {
  vector[H] loglik_h;
  for(i in 1:H){
    loglik_h[i] = weightednormal_lpdf(y_h[i] | s_h[i], p_h[i], theta, tau, omega, steps1);
  }
}
