functions {
  vector ode(real t,
             vector y,
             vector theta,
             array[] real x,
             array[] int x_int) {
    vector[2] dydt;
    return dydt;
  }
}

data {
  int n;
  int np;
  int nt;
  real t0;
  vector[n] y0;
  array[nt] real ts;
  vector[np] theta;
  array[np] real x_r;
  array[np] int x_i;
}

transformed data {
  array[nt] vector[n] y;
  real rtol = 1.0e-5;
  real atol = 1.0e-8;
  int max_num_steps = 10000;

  y = pmx_ode_adams_ctrl (ode, y0, t0, ts, rtol, atol, max_num_steps, theta, x_r, x_i);
  y = pmx_ode_bdf_ctrl   (ode, y0, t0, ts, rtol, atol, max_num_steps, theta, x_r, x_i);
  y = pmx_ode_rk45_ctrl  (ode, y0, t0, ts, rtol, atol, max_num_steps, theta, x_r, x_i);
  y = pmx_ode_ckrk_ctrl  (ode, y0, t0, ts, rtol, atol, max_num_steps, theta, x_r, x_i);
}

parameters {
  real x_p;
  vector[np] theta_p;
  array[nt] real ts_v;
}

transformed parameters {
  array[nt] vector[n] yp;
  vector[n] y0_p;

  yp = pmx_ode_adams_ctrl (ode, y0_p, t0, ts,   rtol, atol, max_num_steps, theta_p, x_r, x_i);
  yp = pmx_ode_bdf_ctrl   (ode, y0_p, t0, ts,   rtol, atol, max_num_steps, theta_p, x_r, x_i);
  yp = pmx_ode_rk45_ctrl  (ode, y0_p, t0, ts,   rtol, atol, max_num_steps, theta_p, x_r, x_i);

  yp = pmx_ode_adams_ctrl (ode, y0_p, t0, ts_v, rtol, atol, max_num_steps, theta_p, x_r, x_i);
  yp = pmx_ode_bdf_ctrl   (ode, y0_p, t0, ts_v, rtol, atol, max_num_steps, theta_p, x_r, x_i);
  yp = pmx_ode_rk45_ctrl  (ode, y0_p, t0, ts_v, rtol, atol, max_num_steps, theta_p, x_r, x_i);
}

model {
  x_p ~ normal(0,1);
}
