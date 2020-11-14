data{
  int<lower = 1> nt;  // number of events
  int<lower = 1> cmt[nt];
  int evid[nt];
  int addl[nt];
  int ss[nt];
  real amt[nt];
  real time[nt];
  real rate[nt];
  real ii[nt];
}

transformed data {
  int nTheta = 5;  // number of ODE parameters in Two Compartment Model
  int nCmt = 3;  // number of compartments in model
  real biovar[nCmt] = {1.0, 1.0, 1.0};
  real tlag[nCmt] = {0.0, 0.0, 0.0};
}

parameters {
  real<lower = 0> CL;
  real<lower = 0> Q;
  real<lower = 0> V1;
  real<lower = 0> V2;
  real<lower = 0> ka;
}

transformed parameters{
  real theta[nTheta] = {CL, Q, V1, V2, ka};
  row_vector<lower = 0>[nt] cHat;
  matrix<lower = 0>[nCmt, nt] x;

  // PKModelTwoCpt takes in the NONMEM data, followed by the parameter
  // arrays abd returns a matrix with the predicted amount in each 
  // compartment at each event.
  x = pmx_solve_twocpt(time, amt, rate, ii, evid, cmt, addl, ss,
                       theta, biovar, tlag);

  cHat = x[2, :] ./ V1; // we're interested in the amount in the second compartment
}

model{
  // informative prior
  CL ~ lognormal(log(10), 0.25);
  Q ~ lognormal(log(15), 0.5);
  V1 ~ lognormal(log(35), 0.25);
  V2 ~ lognormal(log(105), 0.5);
  ka ~ lognormal(log(2.5), 1);
}
