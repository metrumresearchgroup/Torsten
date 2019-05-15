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
  real t0[n];
  real y0[n];
  int len[np];
  real ts[np * nt];
  real theta[n];
  real x_r[n];
  int x_i[n];
}

transformed data {
  real rtol = 1.0e-5;
  real atol = 1.0e-8;
  int max_num_steps = 10000;
}

parameters {
  real x_p;
  real theta_p[n];  
}

transformed parameters {
  matrix[n, np * nt] yp;
  real y0_p[n];

  yp = pmx_integrate_ode_rk45(ode, y0_p, t0, ts, theta_p, x_r, x_i, rtol, atol, max_num_steps);
}

model {
  x_p ~ normal(0,1);
}
