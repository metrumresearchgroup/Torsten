functions {
  real[] dz_dt(real t,       // time
               real[] z,     // system state {prey, predator}
               real[] theta, // parameters
               real[] x_r,   // unused data
               int[] x_i) {
    real u = z[1];
    real v = z[2];

    real alpha = theta[1];
    real beta = theta[2];
    real gamma = theta[3];
    real delta = theta[4];

    real du_dt = (alpha - beta * v) * u;
    real dv_dt = (-gamma + delta * u) * v;
    return { du_dt, dv_dt };
  }
}
data {
  int<lower = 0> N_subj;      // number of subjects
  int<lower = 0> N;           // number of measurement times
  real ts_0[N];                 // measurement times > 0
  real y0_0[2];     // initial measured populations
  real<lower = 0> y_0[N, 2];    // measured populations
}
transformed data {
  int len[N_subj] = rep_array(N, N_subj);
  real y0[N_subj, 2] = rep_array(y0_0, N_subj);
  real y[N_subj, N, 2] = rep_array(y_0, N_subj);
  real ts[N_subj * N];
  for (i in 1:N_subj) {
    ts[((i-1)*N + 1) : (i*N)] = ts_0;
  }
}

parameters {
  real<lower = 0> theta[N_subj, 4];   // { alpha, beta, gamma, delta }
  real<lower = 0> z_init[N_subj, 2];  // initial population
  real<lower = 0> sigma[N_subj, 2];   // measurement errors
}
transformed parameters {
  matrix[2, N_subj * N] z;
  z = pmx_integrate_ode_group_rk45(dz_dt, z_init, 0, len, ts, theta, rep_array(rep_array(0.0, 0), N_subj), rep_array(rep_array(0, 0),N_subj));
}
model {
  for (isub in 1:N_subj) {
    theta[isub, {1, 3}] ~ normal(1, 0.5);
    theta[isub, {2, 4}] ~ normal(0.05, 0.05);
    sigma[isub] ~ lognormal(-1, 1);
    z_init[isub] ~ lognormal(10, 1);
    for (k in 1:2) {
      y0[isub, k] ~ lognormal(log(z_init[isub, k]), sigma[isub, k]);
      y[isub, , k] ~ lognormal(log(z[k, ((isub-1)*N + 1):(isub*N)]), sigma[isub, k]);
    }
  }
}
