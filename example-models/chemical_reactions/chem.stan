functions{
  vector reaction(real t, vector x, array[] real p, array[] real r,
                  array[] int i){
    vector[3] dxdt;
    real p1 = p[1];
    real p2 = p[2];
    real p3 = p[3];
    dxdt[1] = -p1*x[1] + p2*x[2]*x[3];
    dxdt[2] =  p1*x[1] - p2*x[2]*x[3] - p3*(x[2])^2;
    dxdt[3] =  p3*(x[2])^2;
    return dxdt;
  }
}

data {
  int<lower=1> nsub;
  array[nsub] int<lower=1> len;  
  int<lower=1> ntot;  
  array[ntot] real ts;
  array[ntot] real obs;
}

transformed data {
  array[nsub] int i1;
  array[nsub] int i2;
  real t0 = 0.0;
  array[0] real xr;
  array[0] int xi;
  array[3] real theta = {0.04, 1.0e4, 3.0e7};
  i1[1] = 1;
  i2[1] = len[1];
  for (i in 2:nsub) {
    i1[i] = i2[i-1] + 1;
    i2[i] = i1[i] + len[i] - 1;
  }
}

parameters {
  /*  p1=0.04, p2=1e4, and p3=3e7 */
  real<lower = 0> y0_mu;
  array[nsub] real<lower = 0> y0_1;
  real<lower = 0> sigma;
}

transformed parameters {
  vector[3] y0;
  array[ntot] vector[3] x;
  array[ntot] real x3;
  for (i in 1:nsub) {
    y0[1] = y0_1[i];
    y0[2] = 0.0;
    y0[3] = 0.0;
    x[i1[i]:i2[i], ] = pmx_ode_bdf_ctrl(reaction, y0, t0, ts[i1[i]:i2[i]], 1.e-4, 1.e-8, 10000, theta, xr, xi);
  }
  x3 = x[ , 3];
}

model {
  y0_mu ~ lognormal(log(2.0), 0.5);
  for (i in 1:nsub) {
    y0_1[i] ~ lognormal(y0_mu, 0.5);    
  }
  sigma ~ cauchy(0, 1); 
  obs ~ lognormal(log(x3), sigma);
}
