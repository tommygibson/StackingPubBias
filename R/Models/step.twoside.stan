// Recreating the model in Maier 2021 RoBMA
// Weighted normal distribution is the key, which we need to define
// can handle any set of steps, which are input as data
// two-sided p-values!

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
      crits[t] = (fabs(inv_Phi(steps[t] / 2)) * sigma); // define quantiles associated with steps (p-value cutoffs)
    }

    // first denominator integral
    // two-sided, so middle chunk of normal distrubution
    denoms[1] = exp(normal_lcdf(crits[1] | theta, tot_sd)) - exp(normal_lcdf(-crits[1] | theta, tot_sd));
    denom_sum = denoms[1];
    if(K > 2) {
      for(k in 2:(K - 1)){ 
        // middle denominator integrals
        denoms[k] = exp(normal_lcdf(crits[k] | theta, tot_sd)) - exp(normal_lcdf(-crits[k] | theta, tot_sd)) - denom_sum;
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
  int<lower = 0> M; // number of steps
  real y[S]; // observed effect size
  real<lower = 0> s[S]; // observed standard error
  vector<lower = 0, upper = 1>[M] steps;
}

transformed data {
  real dev[S];
  real p[S];
  vector[M] steps1 = sort_desc(steps);
  for(i in 1:S){
    dev[i] = fabs(y[i]) / s[i];
    p[i] = 2 * (1 - exp(normal_lcdf(dev[i] | 0, 1))); // two-sided p-values
  }
}


parameters {
  real theta;
  real<lower = 0> tau;
  simplex[M + 1] pre_omega;
}

transformed parameters {
  vector<lower = 0, upper = 1>[M + 1] omega = cumulative_sum(pre_omega);
}

model {

  theta ~ normal(0, 1); // prior for mean hyperparameter
  tau ~ cauchy(0, 1)T[0,]; // half-cauchy for heterogeneity sd
  // weights \omega are cumulative dirichlet
  pre_omega ~ dirichlet(rep_vector(1, (M + 1)));
  
  for(i in 1:S){
    y[i] ~ weightednormal(s[i], p[i], theta, tau, omega, steps1);
  }
}

generated quantities {
  vector[S] loglik;
  for(i in 1:S){
    loglik[i] = weightednormal_lpdf(y[i] | s[i], p[i], theta, tau, omega, steps1);
  }
}
