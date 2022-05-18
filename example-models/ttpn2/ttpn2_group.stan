functions {
  vector oneCptPNODE(real t, vector x, array[] real parms,
                     array[] real rdummy, array[] int idummy) {
    vector[3] dxdt;
    real CL = parms[1];
    real V = parms[2];
    real ke0 = parms[3];
    real alpha = parms[4];
    real beta = parms[5];
    real Edrug;
    real hazard;
    
    dxdt[1] = -(CL / V) * x[1];
    dxdt[2] = ke0 * (x[1] / V - x[2]);
    Edrug = alpha * x[2];
    if (t == 0) {
      hazard = 0;
    } else {
      hazard = beta * Edrug ^ beta * t ^ (beta - 1);
    }
    dxdt[3] = hazard;
    
    return dxdt;
  }
}
data {
  int<lower=1> nId;
  int<lower=1> nt;
  int<lower=1> nPNObs;
  int<lower=1> nPNCens;
  array[nPNObs] int<lower=1> iPNObs;
  array[nPNCens] int<lower=1> iPNCens;
  array[nt] real<lower=0> amt;
  array[nt] real<lower=0> rate;
  array[nt] real<lower=0> ii;
  array[nt] int<lower=0> addl;
  array[nt] int<lower=1> cmt;
  array[nt] int<lower=0> evid;
  array[nId] int<lower=1> start;
  array[nId] int<lower=1> end;
  array[nt] real<lower=0> time;
  array[nId] real<lower=0> CL;
  array[nId] real<lower=0> V;
}
transformed data {
  array[nId] int<lower=0> len;
  array[nt] int<lower=0> ss = rep_array(0, nt);
  int<lower=1> nCmt = 3;
  array[nCmt] real F = rep_array(1.0, nCmt);
  array[nCmt] real tLag = rep_array(0.0, nCmt);
  
  int ntPred = 253;
  real dt = 12;
  array[ntPred] real tPred;
  array[ntPred] real ratePred = rep_array(0.0, ntPred);
  array[ntPred] real iiPred;
  array[ntPred] int addlPred;
  array[ntPred] int cmtPred = rep_array(1, ntPred);
  array[ntPred] int evidPred;
  array[ntPred] int ssPred = rep_array(0, ntPred);
  
  iiPred[1] = 21 * 24;
  iiPred[2 : ntPred] = rep_array(0.0, ntPred - 1);
  addlPred[1] = 5;
  addlPred[2 : ntPred] = rep_array(0, ntPred - 1);
  evidPred[1] = 1;
  evidPred[2 : ntPred] = rep_array(0, ntPred - 1);
  for (i in 1 : ntPred) 
    tPred[i] = dt * (i - 1);
  
  for (i in 1 : nId) 
    len[i] = end[i] - start[i] + 1;
}
parameters {
  real<lower=0> ke0;
  real<lower=0> alpha;
  real<lower=0> beta;
}
transformed parameters {
  vector<lower=0>[nPNObs] survObs;
  row_vector<lower=0>[nPNObs] EdrugObs;
  vector<lower=0>[nPNObs] hazardObs;
  vector<lower=0>[nPNCens] survCens;
  matrix<lower=0>[3, nt] x;
  array[nId, 5] real<lower=0> parms;
  
  for (j in 1 : nId) {
    parms[j,  : ] = {CL[j], V[j], ke0, alpha, beta};
  }
  
  x = pmx_solve_group_rk45(oneCptPNODE, 3, len, time, amt, rate, ii, evid,
                           cmt, addl, ss, parms, rep_array(F, nId),
                           rep_array(tLag, nId), 1e-6, 1e-6, 1e8);
  
  for (i in 1 : nPNObs) 
    survObs[i] = fmax(machine_precision(), exp(-x[3, iPNObs[i]]));
  EdrugObs = alpha * x[2, iPNObs];
  for (i in 1 : nPNObs) 
    hazardObs[i] = fmax(machine_precision(),
                        beta * EdrugObs[i] ^ beta
                        * time[iPNObs[i]] ^ (beta - 1));
  for (i in 1 : nPNCens) 
    survCens[i] = fmax(machine_precision(), exp(-x[3, iPNCens[i]]));
}
model {
  ke0 ~ normal(0, 0.0005);
  alpha ~ normal(0, 0.000003);
  beta ~ normal(0, 1.5);
  
  target += log(hazardObs .* survObs); // observed PN event log likelihood
  target += log(survCens); // censored PN event log likelihood
}
