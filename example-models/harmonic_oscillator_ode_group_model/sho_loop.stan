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
  real<lower=0> mu;
  real<lower=0> tau;
  real<lower=0> sigma[N_subj, n];
  real<lower=0> theta[N_subj, 1];
}

transformed parameters {
  real y[N_subj, nt, n];
  for (i in 1:N_subj) {
    y[i] = integrate_ode_bdf(sho, y0[i], 0.0, ts[1:nt], theta[i], rep_array(0.0, 0), rep_array(0, 0));
  }
}

model {
  mu ~ normal(3, 0.5);
  tau ~ normal(0.5, 0.1);
  for (i in 1:N_subj) {
    theta[i, 1] ~ normal(mu, tau);
    sigma[i] ~ cauchy(0, 1);
    y1_obs[((i-1)*nt+1):(i*nt)] ~ normal(y[i, 1:nt, 1], sigma[i, 1]);
    y2_obs[((i-1)*nt+1):(i*nt)] ~ normal(y[i, 1:nt, 2], sigma[i, 2]);
  }
}
