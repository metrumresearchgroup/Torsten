functions {
  real[] ode(real t,
             real[] y,
             real[] theta,
             real[] x,
             int[] x_int) {
    real dydt[2];
    return dydt;
  }
}

data {
  int n;
  int np;
  int nt;
  real t0;
  real y0[n, np];
  int len[np];
  real ts[np * nt];
  real theta[n, np];
  real x_r[n, np];
  int x_i[n, np];
}

transformed data {
  matrix[n, np * nt] y;

  y = pmx_integrate_ode_group_adams (ode, y0, t0, len, ts, theta, x_r, x_i);
  y = pmx_integrate_ode_group_bdf   (ode, y0, t0, len, ts, theta, x_r, x_i);
  y = pmx_integrate_ode_group_rk45  (ode, y0, t0, len, ts, theta, x_r, x_i);
}

parameters {
  real x_p;
  real theta_p[n, np];  
}

transformed parameters {
  matrix[n, np * nt] yp;
  real y0_p[n, np];

  yp = pmx_integrate_ode_group_adams (ode, y0_p, t0, len, ts, theta_p, x_r, x_i);
  yp = pmx_integrate_ode_group_bdf   (ode, y0_p, t0, len, ts, theta_p, x_r, x_i);
  yp = pmx_integrate_ode_group_rk45  (ode, y0_p, t0, len, ts, theta_p, x_r, x_i);
}

model {
  x_p ~ normal(0,1);
}
