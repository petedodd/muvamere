/* Hierarchical multivariate regression: non-centered regularized horseshoe on the global
   correlation (as mvn_inferH_sparse_rhs.stan), but the rho-blend + free-LKJ Omega_local[i] is
   replaced by a direct, non-centered, estimated-strength deviation of each study's correlation
   from Omega_global:
     Omega_s[i][j,k] = Omega_global[j,k] + kappa * eta[i][loc[j,k]],   eta[i] ~ N(0,1) iid
     kappa ~ half-normal(0, kappa_prior_scale)
   with Omega_s[i] rejected if not positive definite (same mechanism as Omega_global below).

   The rho-blend model (mvn_inferH_sparse_rhs.stan) has a genuine non-identifiability issue
   because Omega_local[i] was a free per-study corr_matrix, for any (rho, Omega_global) there is a
   valid Omega_local[i] that reproduces study i's data-implied correlation exactly, at zero prior
   cost (the LKJ prior on Omega_local[i] doesn't reference Omega_global at all). This created a
   second posterior mode where Omega_global collapses to ~0 (favoured by the horseshoe) and each
   study's correlation is "absorbed" independently into Omega_local[i], with rho left unidentified
   at its prior mean. That mode is invisible to per-study fit checks (Sigs[i] is fine either way)
   but corrupts exactly what mvn_generate_AP() uses to simulate an unobserved cohort (rho,
   Omega_global), since a new cohort's correlation is a fresh draw blended via those two
   quantities, not anything from the fitted Omega_local[i]'s.

   Betas are non-centered (Betas[i] = BetaM + BetaS .* Bz[i]); tau is centered (with the
   truncation-normaliser term, as in the shipped hierarchical models). Likelihood is evaluated
   study-by-study (records must be sorted by study). */
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
  real<lower=0> p0;         //prior guess for the number of non-zero correlations, 0 < p0 < choose(NV,2)
  real<lower=0> slab_scale; //scale of the Student-t slab for large correlations
  real<lower=0> slab_df;    //degrees of freedom of the slab
  real<lower=0> kappa_prior_scale; //scale of the half-normal prior on kappa (local-deviation strength)
}
transformed data{
  int ns[Nstudies] = rep_array(0, Nstudies);  //records per study
  matrix[NP,Nrecords] Xt = X';
  matrix[NV,Nrecords] Yt = Y';
  int D_R = choose(NV,2);   //number of correlations
  real T_scale;
  int loc[NV,NV];//location key for upper triangle of Omega_global / each study's deviation
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
  vector[D_R] zg;                         //non-centered global correlations (see header)
  vector<lower=0>[D_R] lam;               //horseshoe local scales for Omega_global
  real<lower=0> T;                        //horseshoe global scale for Omega_global
  real<lower=0> caux;                     //slab: c = slab_scale*sqrt(caux)
  vector[D_R] eta[Nstudies];              //non-centered per-study deviations from Omega_global
  real<lower=0> kappa;                    //deviation strength (shared across studies and pairs)
  vector<lower=0>[NV] tau[Nstudies];      //scale
  matrix[NP,NV] Bz[Nstudies];             //non-centered regression-coef deviations
  matrix[NP,NV] BetaM;                    //global regression coef means
  matrix<lower=0>[NP,NV] BetaS;           //global regression coef SDs
  vector<lower=0>[NV] taum;               //global cor scale means
  vector<lower=0>[NV] sigt;               //global cor scale SDs
}
transformed parameters{
  real<lower=0> c = slab_scale * sqrt(caux);
  vector<lower=0>[D_R] lam_tilde = sqrt( c^2 * square(lam) ./ (c^2 + T^2 * square(lam)) );
  matrix[NV,NV] Omega_global = diag_matrix(rep_vector(1.0, NV));
  matrix[NV,NV] Omega_s[Nstudies];
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
    // reject proposals where Omega_global is not positive definite
    matrix[NV,NV] Lg = cholesky_decompose(Omega_global);
  }
  for(s in 1:Nstudies){
    Omega_s[s] = diag_matrix(rep_vector(1.0, NV));
    for(i in 1:(NV-1)){
      for(j in (i+1):NV){
        real v = Omega_global[i,j] + kappa * eta[s][loc[i,j]];
        Omega_s[s][i,j] = v;
        Omega_s[s][j,i] = v;
      }
    }
    {
      // reject proposals where Omega_s[s] is not positive definite
      matrix[NV,NV] Ls = cholesky_decompose(Omega_s[s]);
    }
    Betas[s] = BetaM + BetaS .* Bz[s];
    Sigs[s] = quad_form_diag(Omega_s[s], tau[s]);
  }
}
model{

  //top-level priors
  //regularized horseshoe on Omega_global's off-diagonal upper triangle, non-centered
  T ~ cauchy(0, T_scale);
  lam ~ cauchy(0,1);
  caux ~ inv_gamma(0.5 * slab_df, 0.5 * slab_df);
  zg ~ std_normal();
  kappa ~ normal(0, kappa_prior_scale); //half-normal (kappa has lower=0)
  taum ~ normal(0, tauM_prior_sd);
  sigt ~ normal(0, tauS_prior_sd);
  to_vector(BetaM) ~ normal(0, betaM_prior_sd);
  to_vector(BetaS) ~ normal(0, betaS_prior_sd);

  //study-level likelihood
  for(i in 1:Nstudies){
    tau[i] ~ normal(taum,sigt);
    to_vector(eta[i]) ~ std_normal();
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
