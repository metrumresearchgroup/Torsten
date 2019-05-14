functions {
  real foo(real t,
           real[] theta,
           real[] x,
           int[] x_int) {
    return 1.0;
  }
}
data {
  int<lower=1> T;
  real t0;
  real t1;
  real theta[1];
}
transformed data {
  real x[0];
  int x_int[0];
}
model {
}
generated quantities {
  real y_hat;
  y_hat = univariate_integral_rk45(foo, t0, t1, theta, x, x_int );
}
