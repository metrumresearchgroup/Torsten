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
  real rtol = 1.0e-5;
  real atol = 1.0e-8;
  int max_num_steps = 10000;

  y = pmx_integrate_ode_adams (ode, y0, t0, ts, theta, x_r, x_i, rtol, atol, max_num_steps);
  y = pmx_integrate_ode_bdf   (ode, y0, t0, ts, theta, x_r, x_i, rtol, atol, max_num_steps);
  y = pmx_integrate_ode_rk45  (ode, y0, t0, ts, theta, x_r, x_i, rtol, atol, max_num_steps);
}

parameters {
  real x_p;
  real theta_p[n];  
  real ts_v[np * nt];
}

transformed parameters {
  real yp[n, np * nt];
  real y0_p[n];

  yp = pmx_integrate_ode_adams (ode, y0_p, t0, ts, theta_p, x_r, x_i, rtol, atol, max_num_steps);
  yp = pmx_integrate_ode_bdf   (ode, y0_p, t0, ts, theta_p, x_r, x_i, rtol, atol, max_num_steps);
  yp = pmx_integrate_ode_rk45  (ode, y0_p, t0, ts, theta_p, x_r, x_i, rtol, atol, max_num_steps);

  yp = pmx_integrate_ode_adams (ode, y0_p, t0, ts_v, theta_p, x_r, x_i, rtol, atol, max_num_steps);
  yp = pmx_integrate_ode_bdf   (ode, y0_p, t0, ts_v, theta_p, x_r, x_i, rtol, atol, max_num_steps);
  yp = pmx_integrate_ode_rk45  (ode, y0_p, t0, ts_v, theta_p, x_r, x_i, rtol, atol, max_num_steps);
}

model {
  x_p ~ normal(0,1);
}
