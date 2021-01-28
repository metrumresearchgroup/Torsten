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
  int len[np];
  real ts[np * nt];
  real theta[n];
  real x_r[n];
  int x_i[n];
}

parameters {
  real x_p;
  real theta_p;
}

transformed parameters {
  matrix[n, np * nt] yp;
  real y0_p[n];

  yp = pmx_integrate_ode_adams(ode, y0_p, t0, ts, theta_p, x_r, x_i);
}

model {
  x_p ~ normal(0,1);
}
