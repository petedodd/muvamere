/* Hierarchical multivariate regression, deviation-penalty ("rhs_kappa")
   model: non-centered regularized horseshoe on the global correlation
   Omega_global, and each study's correlation deviating from it with an
   estimated strength kappa:
     Omega_s[i][j,k] = rr[i][loc[j,k]],  rr[i] ~ normal(Omega_global, kappa)
     kappa ~ half-normal(0, kappa_prior_scale)
   with Omega_s[i] rejected if not positive definite (as for Omega_global).
   Scales use a non-centered log-normal hierarchy:
     log tau[i] = ltaum + lsig .* tz[i],  tz[i] ~ N(0,1).
   Betas are non-centered (Betas[i] = BetaM + BetaS .* Bz[i]).

   History: this replaced the rho-blend + free per-study LKJ
   local layer (still available as prior = "rhs"), which had a spurious
   all-zero Omega_global mode. The first version of this model used
   NON-centered deviations (Omega_global + kappa * eta[i]) and a centered,
   truncated-normal tau hierarchy; a one-change-at-a-time ablation
   (apopwork/apop_study/bench/target/kdiv/) located its divergences in the
   horseshoe block (coupled to every study through the non-centered
   deviations) and its slow mixing in the tau hierarchy. Centered
   deviations + log-normal non-centered tau fixed both: fewer divergences
   and 1.6-2.7x faster across sparse/medium/dense truths, with identical
   accuracy; only with weak data (small studies, large kappa) do centered
   deviations sometimes diverge more, while still mixing better.

   Likelihood is evaluated study-by-study (records must be sorted by
   study). */
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
  real tauM_prior_sd; //prior SD for the mean of log(tau)
  real tauS_prior_sd; //prior SD (half-normal) for the SD of log(tau)
  //prior guess for the number of non-zero correlations, 0 < p0 < choose(NV,2)
  real<lower=0> p0;
  real<lower=0> slab_scale; //Student-t slab scale for large correlations
  real<lower=0> slab_df;    //degrees of freedom of the slab
  //scale of the half-normal prior on kappa (local-deviation strength)
  real<lower=0> kappa_prior_scale;
}
transformed data{
  int ns[Nstudies] = rep_array(0, Nstudies);  //records per study
  matrix[NP,Nrecords] Xt = X';
  matrix[NV,Nrecords] Yt = Y';
  int D_R = choose(NV,2);   //number of correlations
  real T_scale;
  //location key for the upper triangle of Omega_global / study correlations
  int loc[NV,NV];
  int k=1;
  if(NV < 2) reject("NV must be at least 2");
  if(p0 >= D_R) reject("p0 must be < choose(NV,2) = ", D_R);
  T_scale = p0 / (D_R - p0) / sqrt(Nrecords);
  for(n in 1:Nrecords){
    if(study[n] < 1 || study[n] > Nstudies)
      reject("study id out of range at record ", n);
    if(n > 1 && study[n] < study[n-1])
      reject("records must be sorted by study; violated at record ", n);
    ns[study[n]] += 1;
  }
  for(i in 1:Nstudies){
    if(ns[i] == 0)
      reject("study ", i, " has no records (ids must be 1..Nstudies)");
  }
  for(i in 1:(NV-1)){
    for(j in (i+1):NV){
      loc[i,j] = k;
      k = k + 1;
    }
  }
}
parameters{
  vector[D_R] zg;                  //non-centered global correlations
  vector<lower=0>[D_R] lam;        //horseshoe local scales
  real<lower=0> T;                 //horseshoe global scale
  real<lower=0> caux;              //slab: c = slab_scale*sqrt(caux)
  vector[D_R] rr[Nstudies];        //centered per-study correlations
  real<lower=0> kappa;             //deviation strength (all studies/pairs)
  vector[NV] tz[Nstudies];         //non-centered log-scale deviations
  matrix[NP,NV] Bz[Nstudies];      //non-centered regression-coef deviations
  matrix[NP,NV] BetaM;             //global regression coef means
  matrix<lower=0>[NP,NV] BetaS;    //global regression coef SDs
  vector[NV] ltaum;                //global log-scale means
  vector<lower=0>[NV] lsig;        //global log-scale SDs
}
transformed parameters{
  real<lower=0> c = slab_scale * sqrt(caux);
  vector<lower=0>[D_R] lam_tilde =
    sqrt(c^2 * square(lam) ./ (c^2 + T^2 * square(lam)));
  matrix[NV,NV] Omega_global = diag_matrix(rep_vector(1.0, NV));
  matrix[NV,NV] Omega_s[Nstudies];
  matrix[NV,NV] Sigs[Nstudies];
  vector[NV] tau[Nstudies];
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
        real v = rr[s][loc[i,j]];
        Omega_s[s][i,j] = v;
        Omega_s[s][j,i] = v;
      }
    }
    {
      // reject proposals where Omega_s[s] is not positive definite
      matrix[NV,NV] Ls = cholesky_decompose(Omega_s[s]);
    }
    Betas[s] = BetaM + BetaS .* Bz[s];
    tau[s] = exp(ltaum + lsig .* tz[s]);
    Sigs[s] = quad_form_diag(Omega_s[s], tau[s]);
  }
}
model{
  vector[D_R] og;
  for(i in 1:(NV-1)) for(j in (i+1):NV) og[loc[i,j]] = Omega_global[i,j];

  //top-level priors
  //regularized horseshoe on Omega_global's upper triangle, non-centered
  T ~ cauchy(0, T_scale);
  lam ~ cauchy(0,1);
  caux ~ inv_gamma(0.5 * slab_df, 0.5 * slab_df);
  zg ~ std_normal();
  kappa ~ normal(0, kappa_prior_scale); //half-normal (kappa has lower=0)
  ltaum ~ normal(0, tauM_prior_sd);
  lsig ~ normal(0, tauS_prior_sd); //half-normal (lsig has lower=0)
  to_vector(BetaM) ~ normal(0, betaM_prior_sd);
  to_vector(BetaS) ~ normal(0, betaS_prior_sd);

  //study-level likelihood
  for(i in 1:Nstudies){
    tz[i] ~ std_normal();
    rr[i] ~ normal(og, kappa);
    to_vector(Bz[i]) ~ std_normal();
  }

  //individual-level likelihood, one block per study
  {
    int pos = 1;
    for(i in 1:Nstudies){
      int a = pos;
      int b = pos + ns[i] - 1;
      matrix[NV,NV] L = cholesky_decompose(Sigs[i]);
      //residuals for the whole block, NV x ns[i]: Y_i' - Betas[i]' * X_i'
      matrix[NV,ns[i]] Z =
        mdivide_left_tri_low(L, Yt[:, a:b] - Betas[i]' * Xt[:, a:b]);
      // = sum_n log MVN(Y_n | mu_n, L L'), incl. the -0.5*log(2 pi) terms
      target += -0.5 * dot_self(to_vector(Z))
                - ns[i] * sum(log(diagonal(L)))
                - 0.5 * ns[i] * NV * log(2 * pi());
      pos = b + 1;
    }
  }

}
generated quantities{

}
