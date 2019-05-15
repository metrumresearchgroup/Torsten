functions {

/**
 * turn a slicing variable for a ragged array
 * S = {5, 6, 11}
 * into
 * Si = {0, 5, 5 + 6, 5+6+11} + 1
 * such that we can index the ragged array A as
 * A[Si[u] : Si[u+1]-1]
 * for the uth unit
 */
int[] make_slice_index(int[] S) {
  int Si[size(S)+1];
  int cv = 1;
  Si[1] = cv;
  for(i in 1:size(S)) {
    cv = cv + S[i];
    Si[i+1] = cv;
  }
  return(Si);
}

/**
 * calculate the absolute value of a - b in log-space with log(a)
 * and log(b) given. Does so by squaring and taking the root, i.e.
 *   la = log(a)
 *   lb = log(b)
 * sqrt( (a - b)^2 ) = sqrt( a^2 - 2 * a * b + b^2 )
 * <=> 0.5 * log_diff_exp(log_sum_exp(2*la, 2*lb), log(2) + la + lb)
 */
real log_diff_exp_abs(real la, real lb) {
  return(0.5 * log_diff_exp(log_sum_exp(2*la, 2*lb), log(2) + la + lb));
}



real gamma2_overdisp_lpdf(vector y, vector mu, real sigma, real kappa) {
  int N = rows(y);
  vector[N] muSq = square(mu);
  vector[N] sigmaSq = square(sigma) + muSq / kappa;
  return gamma_lpdf(y| muSq ./ sigmaSq, mu ./ sigmaSq);
}

real gamma2_overdisp_array_lpdf(real[] y, vector mu, real sigma, real kappa) {
  int N = size(y);
  vector[N] muSq = square(mu);
  vector[N] sigmaSq = square(sigma) + muSq / kappa;
  return gamma_lpdf(y| muSq ./ sigmaSq, mu ./ sigmaSq);
}

real gamma2_overdisp_cdf(real y, real mu, real sigma, real kappa) {
  real muSq = square(mu);
  real sigmaSq = square(sigma) + muSq / kappa;
  return gamma_cdf(y, muSq ./ sigmaSq, mu ./ sigmaSq);
}

real gamma2_overdisp_rng(real mu, real sigma, real kappa) {
  real muSq = square(mu);
  real sigmaSq = square(sigma) + muSq / kappa;
  return gamma_rng(muSq ./ sigmaSq, mu ./ sigmaSq);
}

  vector pk_1cmt_oral(vector tad, real ldose, real lka, real lCL, real lV) {
    real lk = lCL - lV;
    real k = exp(lk);
    real ka = exp(lka);
    real lkdelta = log_diff_exp_abs(lka, lk);
    real lscale = ldose - lV + lka - lkdelta;
    int N = num_elements(tad);
    vector[N] pk;

    for(i in 1:N) {
      pk[i] = lscale + log_diff_exp_abs(-k * tad[i], -ka * tad[i]);
    }

    return(pk);
  }

  // this version expects that the lag time is never greater than the
  // first tad time.
  vector pk_1cmt_oral_tlagMax(vector tad, real ldose, real llag, real lka, real lCL, real lV) {
    int k = 0;
    int N = num_elements(tad);
    real lag = exp(llag);
    vector[N] pk;

    if(lag > tad[1])
      reject("Tlag must not exceed first observation");
 
    return pk_1cmt_oral(tad - lag, ldose, lka, lCL, lV);
  }

// this version skips according to llag the first elements of tad
// (better use this if llag and tad is data only)
  vector pk_1cmt_oral_tlag(vector tad, real ldose, real llag, real lka, real lCL, real lV) {
    int k = 0;
    int N = num_elements(tad);
    real lag = exp(llag);
    vector[N] pk;
 
    // find out which elements have to be skipped due to lag
    while(tad[k+1] < lag && k+1 != N) {
      k = k + 1;
      pk[k] = -25;
    }

    pk[k+1:N] = pk_1cmt_oral(tail(tad, N-k) - lag, ldose, lka, lCL, lV);
    return(pk);
  }

// work around the problem that we do not want things to land on the AD stack
  real pk_1cmt_oral_t(real tad, real ldose, real lka, real lCL, real lV) {
    real lk = lCL - lV;
    real k = exp(lk);
    real ka = exp(lka);
    real lkdelta = log_diff_exp_abs(lka, lk);
    real lscale = ldose - lV + lka - lkdelta;
    return lscale + log_diff_exp_abs(-k * tad, -ka * tad);
  }
  real pk_1cmt_oral_tlag_t(real t, real ldose, real llag, real lka, real lCL, real lV) {
    real lag = exp(llag);
    if (t < lag)
      return -25;
    return pk_1cmt_oral_t(t-lag , ldose, lka, lCL, lV);
  }

  real[] turnover_kin_inhib_2(real t, real[] R, real[] theta, real[] x_r, int[] x_i) {
    //real ldose = x_r[1];
    //real llag = x_r[2];
    //real lka = x_r[3];
    //real lCl = x_r[4];
    //real lV = x_r[5];
    real lconc = pk_1cmt_oral_tlag_t(t, x_r[1], x_r[2], x_r[3], x_r[4], x_r[5]);
    real lkout = -theta[2];
    real lkin = theta[1] + lkout;
    real lEC50 = theta[3];
    real lS = log_inv_logit(lconc - lEC50);
    return { exp(lkin + log1m_exp(lS)) - R[1] * exp(lkout) };
  }
  
  vector lpdf_subject(vector phi, vector eta, real[] x_r, int[] x_i) {
    vector[3] theta = phi[1:3];
    vector[3] sigma_eta = phi[4:6];
    real sigma_y = phi[7];
    real kappa = phi[8];
    real eta_cov[3] = to_array_1d(theta + eta .* sigma_eta);
    int start = x_i[1];
    int end = x_i[2]-1;
    int Ndv = x_i[3];
    int N = end-start+1;
    real mu_temp[N,1] = pmx_integrate_ode_rk45(turnover_kin_inhib_2,
                                           { exp(eta_cov[1]) }, -1E-4,
                                           x_r[7+start-1 : 7+end-1], // time
                                           eta_cov,
                                           x_r[1:5], // PK parameters
                                           x_i[1:0],
                                           1E-5, 1E-3, 500
                                           );
    vector[N] mu = to_vector(mu_temp[:,1]);
    return [ gamma2_overdisp_array_lpdf(x_r[7+Ndv+start-1 : 7+Ndv+end-1]| mu, sigma_y, kappa * x_r[6] ) ]';
  }
}
data {
  int<lower=1> J;
  vector<lower=0>[J] dose;
  int<lower=1> N[J];
  vector[sum(N)] dv;
  real<lower=0> time[sum(N)];
  matrix[4,J] Eta_est;
  real<lower=0> ref;
}
transformed data {
  int Ndv = sum(N);
  int sidx[J+1] = make_slice_index(N);
  int x_i[J,3];
  real x_r[J,6 + 2*Ndv];
  real t0 = -1E-4;
  real refSq = square(ref);
  vector[3] zero = rep_vector(0.0, 3);
  matrix[3,3] Identity = diag_matrix(rep_vector(1.0,3));

  for(j in 1:J) {
    x_r[j,1] = log(dose[j]);
    x_r[j,2:5] = to_array_1d(Eta_est[:,j]);
    x_r[j,6] = refSq;
    x_r[j,7:7+Ndv-1] = to_array_1d(time);
    x_r[j,7+Ndv:7+2*Ndv-1] = to_array_1d(dv);
    
    x_i[j,1] = sidx[j];
    x_i[j,2] = sidx[j+1];
    x_i[j,3] = Ndv;
  }
}
parameters {
  // log for 1-R0=kin/kout, 2-1/kout, 3-EC50
  vector[3] theta;

  //matrix[3,J] Eta;
  vector[3] Eta[J];
  vector<lower=0>[3]  sigma_eta;

  real<lower=0> sigma_y;
  real<lower=0> kappa;
}
transformed parameters {
}
model {
  vector[8] phi;
  phi[1:3] = theta;
  phi[4:6] = sigma_eta;
  phi[7] = sigma_y;
  phi[8] = kappa;
  
  theta[1] ~ normal(log(80),   log(10)/1.96);
  theta[2] ~ normal(log(30),   log(10)/1.96);
  theta[3] ~ normal(log( 2.5), log(10)/1.96);
  
  sigma_eta ~ normal(0, 0.5);
  kappa ~ gamma(0.2, 0.2);

  Eta ~ multi_normal_cholesky(zero, Identity);

  sigma_y ~ normal(0, 10);

  target += map_rect(lpdf_subject, phi, Eta, x_r, x_i);
}
generated quantities {
  real R0 = exp(theta[1]);
  real kin = exp(theta[1] - theta[2]);
  real kout = exp(-theta[2]);
  real EC50 = exp(theta[3]);
  real sigma_y_ref = sqrt(square(sigma_y) + 1.0/kappa);
}
