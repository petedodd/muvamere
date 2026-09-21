/* Hierarchical multivariate regression with horseshoe-sparsity shrinkage on
   the global correlation structure only.

   also now revised to use block-study implementation of likelihood
*/
data{
  int Nrecords; //number of records/patients
  int Nstudies; //number of studies
  int study[Nrecords]; //which study does each record correspond to?
  int NP;//number of variables
  int NV;//number of variates
  matrix[Nrecords,NP] X; //covariate data
  matrix[Nrecords,NV] Y; //outcomes
  real betaM_prior_sd;//prior for Betas
  real betaS_prior_sd;//prior for Betas
  real tauM_prior_sd; //prior for tau
  real tauS_prior_sd; //prior for tau
  real lkj_local_prior_scale;//prior for local cor
  // Omega_global is horseshoe-shrunk below, not LKJ, so there is no lkj_global_prior_scale here
  real rhoA;//beta parameter for rho
  real rhoB;//beta parameter for rho
}
transformed data{
  int ns[Nstudies] = rep_array(0, Nstudies);  //records per study
  matrix[NP,Nrecords] Xt = X';
  matrix[NV,Nrecords] Yt = Y';
  for(n in 1:Nrecords){
    if(study[n] < 1 || study[n] > Nstudies) reject("study id out of range at record ", n);
    if(n > 1 && study[n] < study[n-1]) reject("records must be sorted by study; violated at record ", n);
    ns[study[n]] += 1;
  }
  for(i in 1:Nstudies){
    if(ns[i] == 0) reject("study ", i, " has no records (ids must be 1..Nstudies)");
  }
  int loc[NV,NV];//location key for upper triangle of Omega_global (horseshoe)
  int k=1;
  for(i in 1:(NV-1)){
    for(j in (i+1):NV){
      loc[i,j] = k;
      k = k + 1;
    }
  }
}
parameters{
  corr_matrix[NV] Omega_local[Nstudies];  //local correlation, plain LKJ
  corr_matrix[NV] Omega_global;           //global correlation, horseshoe-shrunk
  vector<lower=0>[choose(NV,2)] lam;      //horseshoe local scales for Omega_global
  real<lower=0> T;                        //horseshoe global scale for Omega_global
  vector<lower=0>[NV] tau[Nstudies];      //scale
  matrix[NP,NV] Betas[Nstudies];          //regression coefs
  matrix[NP,NV] BetaM;                    //global regression coef means
  matrix<lower=0>[NP,NV] BetaS;           //global regression coef SDs
  vector<lower=0>[NV] taum;               //global cor scale means
  vector<lower=0>[NV] sigt;               //global cor scale SDs
  real<lower=0,upper=1> rho;              //local-global cor interpolant
}
transformed parameters{
  matrix[NV,NV] Sigs[Nstudies];
  for(i in 1:Nstudies){
    Sigs[i] = quad_form_diag(rho * Omega_global + (1-rho) * Omega_local[i], tau[i]);
  }
}
model{

  //top-level priors
  //horseshoe shrinkage on Omega_global's off-diagonal upper triangle
  T ~ cauchy(0,1);
  lam ~ cauchy(0,1);
  for(i in 1:(NV-1)){   //upper
    for(j in (i+1):NV){ //triangle
      Omega_global[i,j] ~ normal(0, lam[loc[i,j]] * T);
    }
  }
  rho ~ beta(rhoA, rhoB); //interpolant
  taum ~ normal(0, tauM_prior_sd);
  sigt ~ normal(0, tauS_prior_sd);
  to_vector(BetaM) ~ normal(0, betaM_prior_sd);
  to_vector(BetaS) ~ normal(0, betaS_prior_sd);

  //study-level likelihood
  for(i in 1:Nstudies){
    tau[i] ~ normal(taum,sigt);
    Omega_local[i] ~ lkj_corr(lkj_local_prior_scale);
    to_vector(Betas[i]) ~ normal(to_vector(BetaM),to_vector(BetaS));
  }
  //tau[i] ~ normal(taum,sigt) is truncated at 0 (tau has lower=0); its normaliser depends on taum, sigt
  target += -Nstudies * sum(log(Phi(taum ./ sigt)));

  //individual-level likelihood, one block per study
  {
    int pos = 1;
    for(i in 1:Nstudies){
      int a = pos;
      int b = pos + ns[i] - 1;
      matrix[NV,NV] L = cholesky_decompose(Sigs[i]);
      //residuals for the whole block, NV x ns[i]: Y_i' - Betas[i]' * X_i'
      matrix[NV,ns[i]] Z = mdivide_left_tri_low(L, Yt[:, a:b] - Betas[i]' * Xt[:, a:b]);
      // = sum_n log MVN(Y_n | mu_n, L L'), including the -0.5*log(2 pi) constants
      target += -0.5 * dot_self(to_vector(Z))
                - ns[i] * sum(log(diagonal(L)))
                - 0.5 * ns[i] * NV * log(2 * pi());
      pos = b + 1;
    }
  }

}
generated quantities{

}
