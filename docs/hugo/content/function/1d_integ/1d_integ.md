+++
title = "Univariate integral"
author = ["Yi Zhang"]
date = 2021-06-25T00:00:00-07:00
lastmod = 2021-06-29T11:51:49-07:00
draft = false
weight = 2010
+++

```stan
real univariate_integral_rk45(f, t0, t1, theta, x_r, x_i)
```

```stan
real univariate_integral_bdf(f, t0, t1, theta, x_r, x_i)
```

Based on the ODE solver capability in Stan, Torsten provides functions
calculating the integral of a univariate function. The integrand function \\(f\\) must follow the signature

```stan
  real f(real t, real[] theta, real[] x_r, int[] x_i) {
    /* ... */
}
```
