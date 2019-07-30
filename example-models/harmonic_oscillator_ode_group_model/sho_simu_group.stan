functions {
  real[] sho(real t, real[] y, real[] theta, real[] x, int[] x_int) {
    real dydt[3];
    dydt[1] =  y[2];
    dydt[2] = -y[1] - theta[1] * y[2];
    return dydt;
  }
}

data {
  int n;
  int np;
  int nt;
  real t0;
  real y0[n];

  /* For subject i, len[i] is the length of its data size in ts */
  int len[np];
  real ts[np * nt];
  real x_r[1];
  int x_i[1];
  /* real y1_obs[np * nt]; */
  real<lower=0> theta_prior[1];
  real<lower=0> theta_sigma_prior[1];
}

transformed data {
  real y0_p[np, n];
  real x_r_p[np, 1];
  int x_i_p[np, 1];
  for (i in 1:np) {
    y0_p[i] = y0;
    x_r_p[i] = x_r;
    x_i_p[i] = x_i;
  }
}

parameters {
  real<lower=0> sigma;
  /* real<lower=0> theta[3]; */
}

/* transformed parameters { */
/*   matrix[n, np * nt] yp; */
/*   yp = pmx_integrate_ode_group_bdf(sho, y0_p, t0, len, ts, theta_p, x_r_p, x_i_p); */
/* } */

model {
  /* for (i in 1:np) { */
  /*   theta_p[i] ~ normal( theta_prior,  theta_sigma_prior ); */
  /* } */
  /* y1_obs ~ normal(yp[1], sigma); */
  sigma ~ normal(0.0,1.0);
}

generated quantities {
  matrix[n, np * nt] yp;
  real<lower=0> theta_p[np, 1];
  real y1[np * nt];
  real y2[np * nt];
  for (i in 1:np) {
    theta_p[i] = theta_prior;
  }
  yp = pmx_integrate_ode_group_adams(sho, y0_p, t0, len, ts, theta_p, x_r_p, x_i_p);
  y1 = normal_rng(yp[1], 1.0);
  y2 = normal_rng(yp[2], 1.0);
}
