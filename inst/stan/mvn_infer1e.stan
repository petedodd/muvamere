data{
  int Nobs;// number of observations
  int NP;// number of variables
  int NV;// number of variates
  matrix[Nobs,NP] X; //covariate data
  matrix[Nobs,NV] Y; //outcomes
  real beta_prior_sd;//prior for Betas, def=5
  real tau_prior_sd; //prior for tau_prior_sd, def=2.5
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
  vector<lower=0>[NV] tau = tau_prior_sd * tan(tau_unif);//gives Cauchy
  matrix[NV, NV] omega = diag_pre_multiply(tau, L_Omega);// * z matrix to get var Sigma
  matrix[NV,Nobs] resid = (Y-X * Beta)';
  matrix[NV,Nobs] z = mdivide_left_tri_low(omega,resid); //this is std normal
}
model{
  //priors
  L_Omega ~ lkj_corr_cholesky(lkj_prior_scale);
  to_vector(Beta) ~ normal(0, beta_prior_sd);

  // likelihood
  to_vector(z) ~ std_normal();//Y = X*beta + (diag(tau)*L_Omega * z')'
}
generated quantities{

}
