+++
title = "Univariate integral of a quadratic function"
author = ["Yi Zhang"]
date = 2021-06-25T00:00:00-07:00
lastmod = 2021-06-28T19:35:23-07:00
draft = false
weight = 2007
+++

integral of a quadratic function.
This example shows how to use `univariate_integral_rk45` to calculate the
integral of a quadratic function.

```stan
functions {
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
  real univar_integral;
  univar_integral = univariate_integral_rk45(func, t0, t1, dtheta,
                          x_r, x_i);
}
/* ... */
```
