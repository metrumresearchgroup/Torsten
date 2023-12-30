functions {
  vector airy(real t, vector y) {
  return [ y[2], -t * y[1] ]';
}
}

data {
 int n;
 array[n] real ts;
 vector[n] yobs;
 real rtol;
 real atol;
}

parameters {
 real<lower=0> y0;
 real<lower=0> sigma;
}

transformed parameters {
 vector[2] yinit = [y0, 0.0]';
 array[n] vector[2] y = pmx_ode_adams_ctrl(airy, yinit, 0.0, ts, rtol, atol, 1000000);
}

model {
  y0 ~ cauchy(0, 0.5);
  sigma ~ normal(0, 0.01);
  for (i in 1:n) {
  yobs[i] ~ normal(y[i][1], sigma);
}
}
