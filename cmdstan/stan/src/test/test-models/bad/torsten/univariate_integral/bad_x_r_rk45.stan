functions {
  real foo(real t,
           real[] theta,     // parameters
           real[] x,         // data
           int[] x_int) {    // integer data
    real res;
    res = 1.0;
    return res;
  }
}
data {
  real t0;
  real t1;
  real x[1, 2];   
  int x_int[0];
  real y;
}
parameters {
  real theta[2];
  real<lower=0> sigma;
}
transformed parameters {
  real y_hat;
  y_hat = univariate_integral_rk45(foo,     // system
                                   t0,      // initial time
                                   t1,      // final time
                                   theta,   // parameters
                                   x,       // data
                                   x_int);  // integer data
}
model {
  y ~ normal(y_hat, sigma);  // independent normal noise
}
