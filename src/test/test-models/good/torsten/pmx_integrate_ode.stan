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
  real y0[n];
  real ts[np * nt];
  real theta[n];
  real x_r[n];
  int x_i[n];
}

transformed data {
  real y[n, np * nt];

  y = pmx_integrate_ode_adams (ode, y0, t0, ts, theta, x_r, x_i);
  y = pmx_integrate_ode_bdf   (ode, y0, t0, ts, theta, x_r, x_i);
}

parameters {
  real x_p;
  real theta_p[n];  
  real ts_v[np * nt];
}

transformed parameters {
  real yp[n, np * nt];
  real y0_p[n];
  yp = pmx_integrate_ode_adams (ode, y0_p, t0, ts, theta_p, x_r, x_i);
  yp = pmx_integrate_ode_bdf   (ode, y0_p, t0, ts, theta_p, x_r, x_i);
  yp = pmx_integrate_ode_adams (ode, y0_p, t0, ts_v, theta_p, x_r, x_i);
  yp = pmx_integrate_ode_bdf   (ode, y0_p, t0, ts_v, theta_p, x_r, x_i);
}

model {
  x_p ~ normal(0,1);
}
