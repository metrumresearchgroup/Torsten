data{
  int nObs;
  array[nObs] real xObs;
  array[nObs] real yObs;
  int nx;
  int nPred;
  array[nPred] real xPred;
}

transformed data{
  real xmin = min(xObs);
  real xmax = max(xObs);
}

parameters{
  array[nx] real y;
  real<lower = 0> sigma;
  simplex[nx - 1] xSimplex;
}

transformed parameters{
  array[nObs] real yHat;
  array[nx] real x;

  x[1] = xmin;
  x[nx] = xmax;
  for(i in 2:(nx-1))
    x[i] = x[i-1] + xSimplex[i-1] * (xmax - xmin);

  for(i in 1:nObs) {
    yHat[i] = pmx_ln_interpolate(xObs[i], x, y);
  }

  yHat = pmx_ln_interpolate(xObs, x, y);
}

model{
}
