+++
title = "Two-compartment model solved by numerical integrator for single patient"
author = ["Yi Zhang"]
date = 2021-06-25T00:00:00-07:00
lastmod = 2021-06-30T11:38:31-07:00
draft = false
weight = 2003
+++

Using `pmx_solve_rk45`, the following example fits a two-compartment model
with first order absorption. User-defined function
`ode_rhs` describes the RHS of the ODEs.

```stan
functions{
  vector ode_rhs(real t, vector x, real[] parms, real[] x_r, int[] x_i){
    real CL = parms[1];
    real Q = parms[2];
    real V1 = parms[3];
    real V2 = parms[4];
    real ka = parms[5];

    real k10 = CL / V1;
    real k12 = Q / V1;
    real k21 = Q / V2;

    vector[3] y;

    y[1] = -ka*x[1];
    y[2] = ka*x[1] - (k10 + k12)*x[2] + k21*x[3];
    y[3] = k12*x[2] - k21*x[3];

    return y;
  }
}
```

We omit `data` and
`model` block as they are identical to [Two-compartment model for single patient]({{< relref "twocpt" >}}) Example.

```stan
transformed data {
  row_vector[nObs] logCObs = log(cObs);
  int nTheta = 5;   // number of parameters
  int nCmt = 3;   // number of compartments
}

parameters {
  real<lower = 0> CL;
  real<lower = 0> Q;
  real<lower = 0> V1;
  real<lower = 0> V2;
  real<lower = 0> ka;
  real<lower = 0> sigma;
}

transformed parameters {
  real theta[nTheta];
  row_vector<lower = 0>[nt] cHat;
  row_vector<lower = 0>[nObs] cHatObs;
  matrix<lower = 0>[3, nt] x;

  theta[1] = CL;
  theta[2] = Q;
  theta[3] = V1;
  theta[4] = V2;
  theta[5] = ka;

  x = pmx_solve_rk45(ode_rhs, 3, time, amt, rate, ii, evid, cmt, addl, ss, theta, 1e-5, 1e-8, 1e5);

  cHat = x[2, ] ./ V1;

  for(i in 1:nObs){
    cHatObs[i] = cHat[iObs[i]];  // predictions for observed data records
  }
}

model{
```
