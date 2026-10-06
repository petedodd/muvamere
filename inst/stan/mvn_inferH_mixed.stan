/* Hierarchical multivariate regression for mixed continuous + binary
   variates: the deviation-penalty model of mvn_inferH_sparse_kappa.stan
   (regularized horseshoe on Omega_global; centered per-study correlations
   rr[s] ~ normal(Omega_global, kappa); non-centered Betas; log-normal scale
   hierarchy) with a multivariate-PROBIT observation layer for the binaries.

   Latent Y*_n ~ MVN(X_n Betas[s], Sigma_s), Sigma_s = D_s Omega_s D_s,
   D_s = diag(tau[s], 1, ..., 1): continuous variates first (observed
   directly), then binary ones with B_nj = 1[Y*_nj > 0] and latent scale
   fixed at 1 (probit identification), so every correlation in Omega_s -
   binary-binary (tetrachoric) and binary-continuous (biserial) included -
   is identified.

   Exact likelihood by GHK (Stan User's Guide, "Multivariate probit"): with
   L = chol(Sigma_s), continuous residuals are whitened (z = L^-1 (y - mu));
   then each binary in turn gets a truncation bound from the z's already
   drawn and one GHK uniform u in (0,1) per record: log-lik += log P(bound),
   z_j = truncated draw. A missing binary (coded -1) is an untruncated step
   (z_j = inv_Phi(u), no log-lik term): exact marginalisation under MAR.
   Missing continuous values are parameters (ymiss).

   Ordinal variates (nlev[j] = K > 2 levels, coded 0..K-1) share the
   binary block: Y_nj = k iff c_k < Y*_nj <= c_{k+1}, with c_0 = -inf,
   c_1 = 0 (so the intercept has the binary meaning, P(Y >= 1) =
   Phi(mu)), c_K = +inf and K-2 free increasing cutpoints c_2..c_{K-1}
   SHARED across studies (study intercepts carry the shifts). The GHK step
   is then a two-sided truncation; binaries (K = 2) use the original
   one-sided code unchanged.

   The u (and ymiss) are per-record quantities that, with the other
   parameters, reveal individual records: the R wrapper does not save them
   by default, and mvn_make_generator() never keeps them.

   Records must be sorted by study. At least one binary or ordinal variate
   is required (use mvn_inferH_sparse_kappa.stan for all-continuous data).
*/
functions{
  // GHK step for latent truncated to (lo, hi] (has_lo/has_hi = 0 for an
  // infinite bound), standardised as z = (lo or hi)/Ljj - pre. Returns
  // [log P(interval), z draw]. A two-sided interval gives z = the
  // u-quantile of the truncated normal (continuous in the parameters),
  // computed from the upper tail when a > 0 for accuracy.
  vector ghk_interval(real pre, real Ljj, int has_lo, real lo, int has_hi,
                      real hi, real u){
    vector[2] out;
    real p;
    real a = has_lo ? lo / Ljj - pre : 0;
    real b = has_hi ? hi / Ljj - pre : 0;
    if(!has_lo){
      p = Phi(b);
      out[2] = inv_Phi(p * u);
    } else if(!has_hi){
      p = Phi(-a);
      out[2] = -inv_Phi(p * u);
    } else if(a > 0){
      // same u-quantile as the branch below (continuous in a), computed
      // from the upper tail: P(Z <= z) = u  <=>  P(-Z < -z) = 1 - u
      real pb = Phi(-b);
      p = Phi(-a) - pb;
      out[2] = -inv_Phi(pb + p * (1 - u));
    } else {
      real pa = Phi(a);
      p = Phi(b) - pa;
      out[2] = inv_Phi(pa + p * u);
    }
    out[1] = log(p);
    return out;
  }
}
data{
  int Nrecords;
  int Nstudies;
  int study[Nrecords];
  int NP;
  int<lower=0> NC;                 //continuous variates (first)
  int<lower=1> NB;                 //binary/ordinal variates (after the cts)
  int<lower=2> nlev[NB];           //levels of each (2 = binary)
  matrix[Nrecords,NP] X;
  matrix[Nrecords,NC] Yc;          //continuous outcomes (missing: any value)
  int<lower=-1> Yb[Nrecords,NB];   //level 0..nlev-1, or -1 = missing
  int<lower=0> Nmiss_c;            //number of missing continuous cells
  int miss_c[Nmiss_c,2];           //(record, column) of each missing cell
  real betaM_prior_sd;
  real betaS_prior_sd;
  real tauM_prior_sd;              //prior SD for the mean of log(tau)
  real tauS_prior_sd;              //prior SD (half-normal) for SD of log(tau)
  //prior guess for the number of non-zero correlations, 0 < p0 < D_R
  real<lower=0> p0;
  real<lower=0> slab_scale;        //Student-t slab scale (large correlations)
  real<lower=0> slab_df;           //degrees of freedom of the slab
  //scale of the half-normal prior on kappa (study-deviation strength)
  real<lower=0> kappa_prior_scale;
  //prior SD (half-normal) of ordinal cutpoint increments (probit scale)
  real<lower=0> cut_prior_sd;
}
transformed data{
  int NV = NC + NB;
  int ns[Nstudies] = rep_array(0, Nstudies);
  int start[Nstudies];
  int D_R = choose(NV,2);
  real T_scale;
  int loc[NV,NV];
  int k = 1;
  int NCUT = 0;                    //free cutpoints, all ordinal variates
  int cs[NB];                      //offset of variate j's cutpoints
  for(j in 1:NB){
    cs[j] = NCUT;
    NCUT += nlev[j] - 2;
  }
  for(n in 1:Nrecords) for(j in 1:NB) if(Yb[n,j] >= nlev[j])
    reject("Yb[", n, ",", j, "] exceeds nlev - 1");
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
  start[1] = 1;
  for(i in 2:Nstudies) start[i] = start[i-1] + ns[i-1];
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
  vector[D_R] rr[Nstudies];        //CENTERED per-study correlations
  real<lower=0> kappa;             //deviation strength (all studies/pairs)
  vector[NC] tz[Nstudies];         //non-centered log-scale deviations
  matrix[NP,NV] Bz[Nstudies];      //non-centered regression-coef deviations
  matrix[NP,NV] BetaM;             //global regression coef means
  matrix<lower=0>[NP,NV] BetaS;    //global regression coef SDs
  vector[NC] ltaum;                //global log-scale means (continuous)
  vector<lower=0>[NC] lsig;        //global log-scale SDs (continuous)
  vector<lower=0,upper=1>[NB] u[Nrecords]; //GHK uniforms (per record)
  vector[Nmiss_c] ymiss;           //missing continuous values (per record)
  vector<lower=0>[NCUT] cinc;      //ordinal cutpoint increments
}
transformed parameters{
  real<lower=0> c = slab_scale * sqrt(caux);
  vector<lower=0>[D_R] lam_tilde =
    sqrt(c^2 * square(lam) ./ (c^2 + T^2 * square(lam)));
  matrix[NV,NV] Omega_global = diag_matrix(rep_vector(1.0, NV));
  matrix[NV,NV] Omega_s[Nstudies];
  matrix[NP,NV] Betas[Nstudies];
  vector[NC] tau[Nstudies];
  vector[NCUT] cuts;               //c_2..c_{K-1} per ordinal variate
  for(j in 1:NB){
    for(m in 1:(nlev[j] - 2)){
      cuts[cs[j] + m] = (m == 1 ? 0 : cuts[cs[j] + m - 1]) + cinc[cs[j] + m];
    }
  }
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
    Betas[s] = BetaM + BetaS .* Bz[s];
    tau[s] = exp(ltaum + lsig .* tz[s]);
  }
}
model{
  matrix[Nrecords,NC] Ycf = Yc;
  vector[D_R] og;
  for(m in 1:Nmiss_c) Ycf[miss_c[m,1], miss_c[m,2]] = ymiss[m];
  for(i in 1:(NV-1)) for(j in (i+1):NV) og[loc[i,j]] = Omega_global[i,j];

  //regularized horseshoe on Omega_global's upper triangle, non-centered
  T ~ cauchy(0, T_scale);
  lam ~ cauchy(0,1);
  caux ~ inv_gamma(0.5 * slab_df, 0.5 * slab_df);
  zg ~ std_normal();
  kappa ~ normal(0, kappa_prior_scale); //half-normal (kappa has lower=0)
  ltaum ~ normal(0, tauM_prior_sd);
  lsig ~ normal(0, tauS_prior_sd);      //half-normal (lsig has lower=0)
  to_vector(BetaM) ~ normal(0, betaM_prior_sd);
  to_vector(BetaS) ~ normal(0, betaS_prior_sd);
  cinc ~ normal(0, cut_prior_sd);       //half-normal (cinc has lower=0)
  for(i in 1:Nstudies){
    tz[i] ~ std_normal();
    rr[i] ~ normal(og, kappa);
    to_vector(Bz[i]) ~ std_normal();
  }

  //likelihood, one block per study
  for(s in 1:Nstudies){
    int a = start[s];
    int b = a + ns[s] - 1;
    matrix[NV,NV] L = cholesky_decompose(
      quad_form_diag(Omega_s[s], append_row(tau[s], rep_vector(1.0, NB))));
    matrix[ns[s],NV] Mu = X[a:b] * Betas[s];
    //continuous block, whole study at once: Zc is NC x ns
    matrix[NC,ns[s]] Zc;
    if(NC > 0){
      Zc = mdivide_left_tri_low(L[1:NC,1:NC], (Ycf[a:b] - Mu[:,1:NC])');
      target += -0.5 * dot_self(to_vector(Zc))
                - ns[s] * sum(log(diagonal(L[1:NC,1:NC])));
    }
    //binary block, GHK record by record
    for(r in 1:ns[s]){
      int n = a + r - 1;
      vector[NV] z;
      if(NC > 0) z[1:NC] = Zc[:,r];
      for(j in 1:NB){
        int jj = NC + j;
        real pre = Mu[r,jj];
        if(jj > 1) pre += L[jj,1:(jj-1)] * z[1:(jj-1)];
        pre /= L[jj,jj];
        if(Yb[n,j] >= 0 && nlev[j] > 2){ //ordinal: (c_y, c_{y+1}]
          int y = Yb[n,j];
          int top = y == nlev[j] - 1;
          real lo = y >= 2 ? cuts[cs[j] + y - 1] : 0;
          real hi = (y >= 1 && !top) ? cuts[cs[j] + y] : 0;
          vector[2] g = ghk_interval(pre, L[jj,jj], y >= 1, lo, !top, hi,
                                     u[n][j]);
          z[jj] = g[2];
          target += g[1];
        } else if(Yb[n,j] == 1){   //latent > 0  <=>  z_jj > -pre
          real p = Phi(pre);
          z[jj] = -inv_Phi(p * u[n][j]);
          target += log(p);
        } else if(Yb[n,j] == 0){   //latent <= 0
          real p = Phi(-pre);
          z[jj] = inv_Phi(p * u[n][j]);
          target += log(p);
        } else {                   //missing: untruncated
          z[jj] = inv_Phi(u[n][j]);
        }
      }
    }
  }
}
