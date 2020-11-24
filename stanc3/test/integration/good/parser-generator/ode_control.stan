functions {
  real[] sho(real t,
             real[] y, 
             real[] theta,
             real[] x,
             int[] x_int) {
    real dydt[2];
    dydt[1] = y[2];
    dydt[2] = -y[1] - theta[1] * y[2];
    return dydt;
  }
}
data {
  int<lower=1> T;
  real y0[2];
  real t0;
  real ts[T];
  real theta[1];
}
transformed data {
  real x[0];
  int x_int[0];
}
model {
}
generated quantities {
  real y_hat[T,3];
  y_hat = integrate_ode_rk45(sho, y0, t0, ts, theta, x, x_int,
                              1e-10, 1e-10, 1e6);
  y_hat = integrate_ode_bdf(sho, y0, t0, ts, theta, x, x_int,
                             1e-10, 1e-10, 1e6);
  y_hat = integrate_ode_adams(sho, y0, t0, ts, theta, x, x_int,
                             1e-10, 1e-10, 1e6);

  // add measurement error
  for (t in 1:T) {
    y_hat[t,1] = y_hat[t,1] + normal_rng(0,0.1);
    y_hat[t,2] = y_hat[t,2] + normal_rng(0,0.1);
    y_hat[t,3] = y_hat[t,3] + normal_rng(0,0.1);
  }
}
