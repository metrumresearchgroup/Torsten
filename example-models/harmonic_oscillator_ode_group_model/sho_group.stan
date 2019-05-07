functions {
  real[] sho(real t, real[] y, real[] theta, real[] x, int[] x_int) {
    real dydt[2];
    dydt[1] =  y[2];
    dydt[2] = -y[1] - theta[1] * y[2];
    return dydt;
  }
}

data {
  int n;                        /* ode system size */
  int N_subj;                   /* nb of subjects */
  int nt;                       /* nb of time steps */
  real ts[N_subj * nt];
  real y1_obs[nt * N_subj];         /* measure rsults */
  real y2_obs[nt * N_subj];         /* measure rsults */
}

transformed data {
  int len[N_subj] = rep_array(nt, N_subj);
  real y0[N_subj, n] = rep_array({0.0, 1.0}, N_subj);
}

parameters {
  real<lower=0> sigma[N_subj, n];
  real<lower=0> theta[N_subj, 1];
}

transformed parameters {
  matrix[n, N_subj * nt] y;
  y = pmx_integrate_ode_group_bdf(sho, y0, 0.0, len, ts, theta, rep_array(rep_array(0.0, 0),N_subj), rep_array(rep_array(0, 0),N_subj));
}

model {
  for (i in 1:N_subj) {
    theta[i, 1] ~ normal(3, 0.5);
    sigma[i] ~ normal(0, 1);
    y1_obs[((i-1)*nt+1):(i*nt)] ~ normal(y[1, ((i-1)*nt+1):(i*nt)], sigma[i, 1]);
    y2_obs[((i-1)*nt+1):(i*nt)] ~ normal(y[2, ((i-1)*nt+1):(i*nt)], sigma[i, 2]);
  }
}
