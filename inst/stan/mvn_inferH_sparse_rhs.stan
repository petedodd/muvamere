/* hierarchical multivariate regression with a non-centered regularized horseshoe
   on the global correlation structure only.

   Regularized horseshoe: Piironen & Vehtari (2017), Electronic Journal of
   Statistics 11(2):5018-5051, doi:10.1214/17-EJS1337SI (arXiv:1707.01694).
     Omega_global[i,j] = z_ij * T * lam_tilde_ij,   z_ij ~ N(0,1)
     lam_tilde^2 = c^2 lam^2 / (c^2 + T^2 lam^2),   lam ~ half-Cauchy(0,1)
     c = slab_scale * sqrt(caux),  caux ~ inv_gamma(slab_df/2, slab_df/2)  (Student-t slab)
     T ~ half-Cauchy(0, T_scale),  T_scale = p0/(D_R - p0)/sqrt(Nrecords),  D_R = choose(NV,2)
   p0 is a prior guess for the number of non-zero (effectively unshrunk) correlations.
   The paper's sigma/sqrt(n) is replaced here by 1/sqrt(Nrecords), the standard error of
   a correlation near zero: a heuristic transplant, not derived in the paper.

   non-centered: Omega_global is built as I + off-diagonal entries and proposals for which it is
   not positive definite are rejected (cholesky_decompose throws). That restricts the prior to the
   set of valid correlation matrices, i.e. the same target as putting the same normal prior on the
   entries of a corr_matrix parameter (the centered form, which samples badly at NV~10).

   This model needs valid initial values: with Stan's default random inits (T ~ 1) Omega_global is
   not positive definite and initialisation fails, and sparse starts (T small, z ~ N(0,1)) can fall
   into a spurious all-zero mode when the true correlations are dense. mvn_infer_mlm_sparse() supplies
   data-informed starts and a start-agreement check; call the model directly only if you do the same.

   Betas are also non-centered (Betas[i] = BetaM + BetaS .* Bz[i]); tau is centered.
   The likelihood is evaluated study-by-study (records must be sorted by study). */
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
  real rhoA;//beta parameter for rho
  real rhoB;//beta parameter for rho
  real<lower=0> p0;         //prior guess for the number of non-zero correlations, 0 < p0 < choose(NV,2)
  real<lower=0> slab_scale; //scale of the Student-t slab for large correlations
  real<lower=0> slab_df;    //degrees of freedom of the slab
}
transformed data{
  int ns[Nstudies] = rep_array(0, Nstudies);  //records per study
  matrix[NP,Nrecords] Xt = X';
  matrix[NV,Nrecords] Yt = Y';
  int D_R = choose(NV,2);   //number of correlations
  real T_scale;
  int loc[NV,NV];//location key for upper triangle of Omega_global
  int k=1;
  if(NV < 2) reject("NV must be at least 2");
  if(p0 >= D_R) reject("p0 must be < choose(NV,2) = ", D_R);
  T_scale = p0 / (D_R - p0) / sqrt(Nrecords);
  for(n in 1:Nrecords){
    if(study[n] < 1 || study[n] > Nstudies) reject("study id out of range at record ", n);
    if(n > 1 && study[n] < study[n-1]) reject("records must be sorted by study; violated at record ", n);
    ns[study[n]] += 1;
  }
  for(i in 1:Nstudies){
    if(ns[i] == 0) reject("study ", i, " has no records (ids must be 1..Nstudies)");
  }
  for(i in 1:(NV-1)){
    for(j in (i+1):NV){
      loc[i,j] = k;
      k = k + 1;
    }
  }
}
parameters{
  corr_matrix[NV] Omega_local[Nstudies];  //local correlation, plain LKJ
  vector[D_R] zg;                         //non-centered global correlations (see header)
  vector<lower=0>[D_R] lam;               //horseshoe local scales
  real<lower=0> T;                        //horseshoe global scale
  real<lower=0> caux;                     //slab: c = slab_scale*sqrt(caux)
  vector<lower=0>[NV] tau[Nstudies];      //scale
  matrix[NP,NV] Bz[Nstudies];             //non-centered regression-coef deviations
  matrix[NP,NV] BetaM;                    //global regression coef means
  matrix<lower=0>[NP,NV] BetaS;           //global regression coef SDs
  vector<lower=0>[NV] taum;               //global cor scale means
  vector<lower=0>[NV] sigt;               //global cor scale SDs
  real<lower=0,upper=1> rho;              //local-global cor interpolant
}
transformed parameters{
  real<lower=0> c = slab_scale * sqrt(caux);
  vector<lower=0>[D_R] lam_tilde = sqrt( c^2 * square(lam) ./ (c^2 + T^2 * square(lam)) );
  matrix[NV,NV] Omega_global = diag_matrix(rep_vector(1.0, NV));
  matrix[NV,NV] Sigs[Nstudies];
  matrix[NP,NV] Betas[Nstudies];          //Betas[i] = BetaM + BetaS .* Bz[i]
  for(i in 1:(NV-1)){
    for(j in (i+1):NV){
      real v = zg[loc[i,j]] * T * lam_tilde[loc[i,j]];
      Omega_global[i,j] = v;
      Omega_global[j,i] = v;
    }
  }
  {
    // reject proposals where Omega_global is not positive definite (see header)
    matrix[NV,NV] Lg = cholesky_decompose(Omega_global);
  }
  for(i in 1:Nstudies){
    Betas[i] = BetaM + BetaS .* Bz[i];
    Sigs[i] = quad_form_diag(rho * Omega_global + (1-rho) * Omega_local[i], tau[i]);
  }
}
model{

  //top-level priors
  //regularized horseshoe on Omega_global's off-diagonal upper triangle, non-centered
  T ~ cauchy(0, T_scale);
  lam ~ cauchy(0,1);
  caux ~ inv_gamma(0.5 * slab_df, 0.5 * slab_df);
  zg ~ std_normal();
  rho ~ beta(rhoA, rhoB); //interpolant
  taum ~ normal(0, tauM_prior_sd);
  sigt ~ normal(0, tauS_prior_sd);
  to_vector(BetaM) ~ normal(0, betaM_prior_sd);
  to_vector(BetaS) ~ normal(0, betaS_prior_sd);

  //study-level likelihood
  for(i in 1:Nstudies){
    tau[i] ~ normal(taum,sigt);
    Omega_local[i] ~ lkj_corr(lkj_local_prior_scale);
    to_vector(Bz[i]) ~ std_normal();
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
