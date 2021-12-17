
data {
  int<lower=0> S; // number of studies
  int y1[S]; // events in risk group
  int y0[S]; // events in nonrisk group
  int n1[S]; // sample size in risk group
  int n0[S]; // sample size in non-risk group
  int N[S]; // total sample size per study
  int L; // number of draws from posterior predicive per mcmc iteration
}

parameters {
  vector[S] beta; // random intercept
  vector[S] delta; // random effect log(OR)
  vector[S] nu; // random effect logit(P(risk factor))
  real beta0; // mean hyperparamters
  real delta0;
  real nu0;
  real<lower=0> sigma_beta; // standard dev hyperparameters
  real<lower=0> sigma_delta;
  real<lower=0> sigma_nu;

}
transformed parameters {
  vector<lower=0>[S] pi1; // P(event | risk factor)
  vector<lower=0>[S] pi0; // P(event | no risk factor)
  vector<lower=0>[S] psi; // P(risk factor)
  
  pi1 = inv_logit(beta + delta / 2);
  pi0 = inv_logit(beta - delta / 2);
  psi = inv_logit(nu);
}

model {
  y1 ~ binomial(n1, pi1); // likelihood
  y0 ~ binomial(n0, pi0);
  n1 ~ binomial(N, psi);
  
  beta ~ normal(beta0, sigma_beta); // random effects
  delta ~ normal(delta0, sigma_delta);
  nu ~ normal(nu0, sigma_nu);
  
  beta0 ~ normal(0, 5); // mean hyperparameters
  delta0 ~ normal(0, 5);
  nu0 ~ normal(0, 5);
  
  sigma_beta ~ cauchy(0, 1); // standard deviation hyperparameters
  sigma_delta ~ cauchy(0, 1);
  sigma_nu ~ cauchy(0, 1);
  
}

generated quantities{
  real betanew[L];
  real deltanew[L];
  real nunew[L];
  real pi1new[L];
  real pi0new[L];
  real psinew[L];
  
  real sensnew[L];
  real PPVmean;
  real sensmean;
  for(l in 1:L){
    
    betanew[l] = normal_rng(beta0, sigma_beta);
    deltanew[l] = normal_rng(delta0, sigma_delta);
    nunew[l] = normal_rng(nu0, sigma_nu);

    pi1new[l] = inv_logit(betanew[l] + deltanew[l] / 2);
    pi0new[l] = inv_logit(betanew[l] - deltanew[l] / 2);    
    psinew[l] = inv_logit(nunew[l]);
    
    sensnew[l] = pi1new[l] * psinew[l] / (pi1new[l] * psinew[l] + pi0new[l] * (1 - psinew[l]));
    
  }
  PPVmean = mean(pi1new);
  sensmean = mean(sensnew); 
}

