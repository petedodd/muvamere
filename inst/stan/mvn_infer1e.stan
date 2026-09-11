data{
  int Nobs;// number of observations
  int NP;// number of variables
  int NV;// number of variates
  matrix[Nobs,NP] X; //covariate data
  matrix[Nobs,NV] Y; //outcomes
  real beta_prior_sd;//prior for Betas, def=5
  real tau_prior_sd; //prior for tau,   def=2.5
  real lkj_prior_scale;//prior for cor, def=2
  int NewObs;// number of new observations to predict for; 0 for none
  matrix[NewObs,NP] Z; //covariate data to predict/simulate on (NP columns, not NV --
  //this dimension was wrong in the original mvn_infer1be.stan prototype, only masked
  //there because its one test happened to use NP==NV)
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
  matrix[NV, NV] L = diag_pre_multiply(tau, L_Omega);// * z matrix to get var Sigma
  matrix[Nobs,NV] Yresid = Y-mu;  
  matrix[NV,Nobs] z = mdivide_left_tri_low(L,Yresid'); //should be std normal: same as X2 = Z * L' in cholsim.stan  
}
model{
  //priors
  L_Omega ~ lkj_corr_cholesky(lkj_prior_scale);
  to_vector(Beta) ~ normal(0, beta_prior_sd);

  //likelihood
  to_vector(z) ~ std_normal();
  target += -Nobs * sum(log(diagonal(L))); //jacobian
}
generated quantities{
  cov_matrix[NV] Sigma = L * L'; //L = diag_pre_multiply(tau, L_Omega) is already Sigma's Cholesky factor
  matrix[NewObs,NV] Ynew; //declared unconditionally (top-level) so it's always a real
  //output, even when NewObs==0 (an empty matrix) -- a variable declared only inside an
  //if-block is LOCAL to that block and would not appear as an output at all, which is
  //the likely reason the original prototype's version of this was flagged unreliable
  {
    matrix[NewObs,NV] munew = Z * Beta; //mean responses for the new/target covariate data
    for(n in 1:NewObs){
      Ynew[n] = multi_normal_cholesky_rng(munew[n]', L)';
    }
  }
}
