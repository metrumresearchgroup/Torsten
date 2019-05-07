

functions {
  real fun_ord0(real t, real[] theta, real[] x_r, int[] x_i) {
    real res;
    res = 2.0;
    return res;
  }
  real fun_ord1(real t, real[] theta, real[] x_r, int[] x_i) {
    real k = 1.2;
    real res;
    res = k * t;
    return res;
  }
  real fun_ord2(real t, real[] theta, real[] x_r, int[] x_i) {
    real a = 2.3;
    real b = 2.0;
    real c = 1.5;
    real res;
    res = a + b * t + c * t * t;
    return res;
  }
}

data {
  real t0;
  real t1;
  real dtheta[2];
  real x_r[0];
  int x_i[0];
}

transformed data {
  real d_univar_integral;

  d_univar_integral = univariate_integral_rk45(fun_ord0, t0, t1, dtheta, x_r, x_i);
  d_univar_integral = univariate_integral_rk45(fun_ord1, t0, t1, dtheta, x_r, x_i);
  d_univar_integral = univariate_integral_rk45(fun_ord2, t0, t1, dtheta, x_r, x_i);

  d_univar_integral = univariate_integral_bdf(fun_ord0, t0, t1, dtheta, x_r, x_i);
  d_univar_integral = univariate_integral_bdf(fun_ord1, t0, t1, dtheta, x_r, x_i);
  d_univar_integral = univariate_integral_bdf(fun_ord2, t0, t1, dtheta, x_r, x_i);
}

parameters {
  real y_p;
  real t0_p;
  real t1_p;
  real ptheta[2];
}

transformed parameters {
  real p_univar_integral;

  p_univar_integral = univariate_integral_rk45(fun_ord0, t0_p, t1_p, ptheta, x_r, x_i);
  p_univar_integral = univariate_integral_rk45(fun_ord1, t0_p, t1_p, ptheta, x_r, x_i);
  p_univar_integral = univariate_integral_rk45(fun_ord2, t0_p, t1_p, ptheta, x_r, x_i);

  p_univar_integral = univariate_integral_bdf(fun_ord0, t0_p, t1_p, ptheta, x_r, x_i);
  p_univar_integral = univariate_integral_bdf(fun_ord1, t0_p, t1_p, ptheta, x_r, x_i);
  p_univar_integral = univariate_integral_bdf(fun_ord2, t0_p, t1_p, ptheta, x_r, x_i);
}

model {
	y_p ~ normal(0,1);
}
