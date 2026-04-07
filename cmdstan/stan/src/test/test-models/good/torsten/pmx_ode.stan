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

  y = pmx_ode_adams (ode, y0, t0, ts, theta, x_r, x_i);
  y = pmx_ode_bdf   (ode, y0, t0, ts, theta, x_r, x_i);
  y = pmx_ode_rk45  (ode, y0, t0, ts, theta, x_r, x_i);
  y = pmx_ode_ckrk  (ode, y0, t0, ts, theta, x_r, x_i);
}

parameters {
  real x_p;
  vector[np] theta_p;
  array[nt] real ts_v;
}

transformed parameters {
  array[nt] vector[n] yp;
  vector[n] y0_p;
  yp = pmx_ode_adams (ode, y0_p, t0, ts, theta_p, x_r, x_i);
  yp = pmx_ode_bdf   (ode, y0_p, t0, ts, theta_p, x_r, x_i);
  yp = pmx_ode_adams (ode, y0_p, t0, ts_v, theta_p, x_r, x_i);
  yp = pmx_ode_bdf   (ode, y0_p, t0, ts_v, theta_p, x_r, x_i);
  yp = pmx_ode_rk45  (ode, y0_p, t0, ts, theta_p, x_r, x_i);
  yp = pmx_ode_ckrk  (ode, y0_p, t0, ts, theta_p, x_r, x_i);
  yp = pmx_ode_rk45  (ode, y0_p, t0, ts_v, theta_p, x_r, x_i);
  yp = pmx_ode_ckrk  (ode, y0_p, t0, ts_v, theta_p, x_r, x_i);
}

model {
  x_p ~ normal(0,1);
}
