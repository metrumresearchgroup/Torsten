functions {
  // define ODE system for two compartmnt model
  vector twoCptModelODE(real t, vector x, array[] real parms,
                              array[] real rate,
                              // in this example, rate is treated as data
                              array[] int dummy) {
    // Parameters
    real CL = parms[1];
    real Q = parms[2];
    real V1 = parms[3];
    real V2 = parms[4];
    real ka = parms[5];
    
    // Re-parametrization
    real k10 = CL / V1;
    real k12 = Q / V1;
    real k21 = Q / V2;
    
    // Return object (derivative)
    vector[3] y; // 1 element per compartment of
    // the model
    
    // PK component of the ODE system
    y[1] = -ka * x[1];
    y[2] = ka * x[1] - (k10 + k12) * x[2] + k21 * x[3];
    y[3] = k12 * x[2] - k21 * x[3];
    
    return y;
  }
}
data {
  int<lower=1> np; /* population size */
  int<lower=1> nt; // number of events
  int<lower=1> nObs; // number of observations
  array[nObs] int<lower=1> iObs; // index of observation
  
  // NONMEM data
  array[np * nt] int<lower=1> cmt;
  array[np * nt] int evid;
  array[np * nt] int addl;
  array[np * nt] int ss;
  array[np * nt] real amt;
  array[np * nt] real time;
  array[np * nt] real rate;
  array[np * nt] real ii;
  
  array[np * nObs] real<lower=0> cObs; // observed concentration (dependent variable)
}
transformed data {
  array[np * nObs] real logCObs;
  array[np] int<lower=1> len;
  array[np] int<lower=1> len_theta;
  array[np] int<lower=1> len_biovar;
  array[np] int<lower=1> len_tlag;
  
  int nTheta = 5; // number of parameters
  int nCmt = 3; // number of compartments
  array[np * nt, nCmt] real biovar;
  array[np * nt, nCmt] real tlag;
  
  logCObs = log(cObs);
  
  for (id in 1 : np) {
    for (j in 1 : nt) {
      for (i in 1 : nCmt) {
        biovar[(id - 1) * nt + j, i] = 1;
        tlag[(id - 1) * nt + j, i] = 0;
      }
    }
    len[id] = nt;
    len_theta[id] = nt;
    len_biovar[id] = nt;
    len_tlag[id] = nt;
  }
}
parameters {
  array[np] real<lower=0> CL;
  array[np] real<lower=0> Q;
  array[np] real<lower=0> V1;
  array[np] real<lower=0> V2;
  array[np] real<lower=0> ka;
  array[np] real<lower=0> sigma;
}
transformed parameters {
  array[np * nt, nTheta] real theta;
  array[np] vector<lower=0>[nt] cHat;
  array[np * nObs] real<lower=0> cHatObs;
  matrix[3, nt * np] x;
  
  for (id in 1 : np) {
    for (it in 1 : nt) {
      theta[(id - 1) * nt + it, 1] = CL[id];
      theta[(id - 1) * nt + it, 2] = Q[id];
      theta[(id - 1) * nt + it, 3] = V1[id];
      theta[(id - 1) * nt + it, 4] = V2[id];
      theta[(id - 1) * nt + it, 5] = ka[id];
    }
  }
  
  x = pmx_solve_group_bdf(twoCptModelODE, 3, len,
                          time, amt, rate, ii, evid, cmt, addl, ss,
                          theta, biovar, tlag);
  
  for (id in 1 : np) {
    for (j in 1 : nt) {
      cHat[id][j] = x[2, (id - 1) * nt + j] ./ V1[id];
    }
  }
  
  for (id in 1 : np) {
    for (i in 1 : nObs) {
      cHatObs[(id - 1) * nObs + i] = cHat[id][iObs[i]]; // predictions for observed data records
    }
  }
}
model {
  // informative prior
  for (id in 1 : np) {
    CL[id] ~ lognormal(log(10), 0.25);
    Q[id] ~ lognormal(log(15), 0.5);
    V1[id] ~ lognormal(log(35), 0.25);
    V2[id] ~ lognormal(log(105), 0.5);
    ka[id] ~ lognormal(log(2.5), 1);
    sigma[id] ~ cauchy(0, 1);
    
    for (i in 1 : nObs) {
      logCObs[(id - 1) * nObs + i] ~ normal(log(cHatObs[(id - 1) * nObs + i]),
                                            sigma[id]);
    }
  }
}
