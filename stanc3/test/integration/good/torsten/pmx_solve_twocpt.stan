data{
  int<lower = 1> nt;  // number of events
  int<lower = 1> nObs;  // number of observation
  array[nObs] int<lower = 1> iObs;  // index of observation
  
  // NONMEM data
  array[nt] int<lower = 1> cmt;
  array[nt] int evid;
  array[nt] int addl;
  array[nt] int ss;
  array[nt] real amt;
  array[nt] real time;
  array[nt] real rate;
  array[nt] real ii;
}

parameters{
  real<lower = 0> CL;
  real<lower = 0> Q;
  real<lower = 0> V1;
  real<lower = 0> V2;
  real<lower = 0> ka;
  real<lower = 0> sigma;
}

transformed parameters{
  array[5] real theta;  // ODE parameters
  row_vector<lower = 0>[nt] cHat;
  matrix<lower = 0>[3, nt] x;

  theta[1] = CL;
  theta[2] = Q;
  theta[3] = V1;
  theta[4] = V2;
  theta[5] = ka;

  x = pmx_solve_twocpt(time, amt, rate, ii, evid, cmt, addl, ss, theta);

  cHat = x[2, :] ./ V1; // we're interested in the amount in the second compartment
}

model{
}
