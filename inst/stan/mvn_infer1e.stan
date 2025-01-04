data{
  int Nobs;// number of observations
  int NP;// number of variables
  int NV;// number of variates
  matrix[Nobs,NP] X; //covariate data
  matrix[Nobs,NV] Y; //outcomes
  real beta_prior_sd;//prior for Betas, def=5
  real tau_prior_sd; //prior for tau,   def=2.5
  real lkj_prior_scale;//prior for cor, def=2
}
transformed data{

}
parameters{
  cholesky_factor_corr[NV] L_Omega; //cholesky factor for variance
  vector<lower=0, upper=pi() / 2>[NV] tau_unif;  // prior scale  
  matrix[NP,NV] Beta;           //regression coefs
}
transformed parameters{
  matrix[Nobs,NV] mu = X * Beta; //mean responses
  vector<lower=0>[NV] tau = tau_prior_sd * tan(tau_unif);//gives Cauchy
  matrix[NV, NV] omega = diag_pre_multiply(tau, L_Omega);// * z matrix to get var Sigma  
}
model{
  //priors
  L_Omega ~ lkj_corr_cholesky(lkj_prior_scale);
  to_vector(Beta) ~ normal(0, beta_prior_sd);

  //likelihood
  for(n in 1:Nobs){
    Y[n] ~ multi_normal_cholesky(mu[n],omega);
  }

}
generated quantities{

}
