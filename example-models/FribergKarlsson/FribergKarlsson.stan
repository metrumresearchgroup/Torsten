functions{
    vector twoCptNeutModelODE(real t, vector x, array[] real parms,
                              array[] real rdummy, array[] int idummy){
    real k10;
    real k12;
    real k21;
    real CL;
    real Q;
    real V1;
    real V2;
    real ka;
    real mtt;
    real circ0;
    real gamma;
    real alpha;
    real ktr;
    vector[8] dxdt;
    real conc;
    real EDrug;
    real transit1;
    real transit2;
    real transit3;
    real circ;
    real prol;

    CL = parms[1];
    Q = parms[2];
    V1 = parms[3];
    V2 = parms[4];
    ka = parms[5];
    mtt = parms[6];	
    circ0 = parms[7];
    gamma = parms[8];
    alpha = parms[9];

    k10 = CL / V1;
    k12 = Q / V1;
    k21 = Q / V2;

    ktr = 4 / mtt;
  
    dxdt[1] = -ka * x[1];
    dxdt[2] = ka * x[1] - (k10 + k12) * x[2] + k21 * x[3];
    dxdt[3] = k12 * x[2] - k21 * x[3];
    conc = x[2]/V1;
    EDrug = alpha * conc;
    // x[4], x[5], x[6], x[7] and x[8] are differences from circ0.
    prol = x[4] + circ0;
    transit1 = x[5] + circ0;
    transit2 = x[6] + circ0;
    transit3 = x[7] + circ0;
    circ = fmax(machine_precision(), x[8] + circ0); // Device for implementing a modeled 
                                                    // initial condition
    dxdt[4] = ktr * prol * ((1 - EDrug) * ((circ0 / circ)^gamma) - 1);
    dxdt[5] = ktr * (prol - transit1);
    dxdt[6] = ktr * (transit1 - transit2);
    dxdt[7] = ktr * (transit2 - transit3);
    dxdt[8] = ktr * (transit3 - circ);

    return dxdt;
  }
}

data{
  int<lower = 1> nt;
  int<lower = 1> nObsPK;
  int<lower = 1> nObsPD;
  array[nObsPK] int<lower = 1> iObsPK;
  array[nObsPD] int<lower = 1> iObsPD;
  array[nt] real<lower = 0> amt;
  array[nt] int<lower = 1> cmt;
  array[nt] int<lower = 0> evid;
  array[nt] real<lower = 0> time;
  array[nt] real<lower = 0> ii;
  array[nt] int<lower = 0> addl;
  array[nt] int<lower = 0> ss;
  array[nt] real rate;
  vector<lower = 0>[nObsPK] cObs;
  vector<lower = 0>[nObsPD] neutObs;
  
  // data for population model
  int<lower = 1> nSubjects;
  array[nSubjects] int<lower = 1> start;
  array[nSubjects] int<lower = 1> end;
  array[nSubjects] real<lower = 0> weight;
  
  real<lower = 0> circ0HatPrior;
  real<lower = 0> circ0HatPriorCV;
  real<lower = 0> mttHatPrior;
  real<lower = 0> mttHatPriorCV;
  real<lower = 0> gammaPrior;
  real<lower = 0> gammaPriorCV;
  real<lower = 0> alphaHatPrior;
  real<lower = 0> alphaHatPriorCV;
}

transformed data{
  vector[nObsPK] logCObs;
  vector[nObsPD] logNeutObs;
//  int idummy[0];
//  real rdummy[0];

  int nTheta;
  int nIIV;

  array[nSubjects] int len;

  logCObs = log(cObs);
  logNeutObs = log(neutObs);
  
  nIIV = 7; // parameters with IIV
  nTheta = 9; // number of parameters

  for(i in 1:nSubjects){
    len[i] = end[i] - start[i] + 1;
  }
}

parameters{
  real<lower = 0> CLHat;
  real<lower = 0> QHat;
  real<lower = 0> V1Hat;
  real<lower = 0> V2Hat;
  real<lower = (CLHat / V1Hat + QHat / V1Hat + QHat / V2Hat +
		sqrt((CLHat / V1Hat + QHat / V1Hat + QHat / V2Hat)^2 -
		     4 * CLHat / V1Hat * QHat / V2Hat)) / 2> kaHat; // ka > lambda_1

  real<lower = 0> mttHat;
  real<lower = 0> circ0Hat;
  real<lower = 0> alphaHat;
  real<lower = 0> gamma;
  real<lower = 0> sigma;
  real<lower = 0> sigmaNeut;
  
  // IIV parameters
  cholesky_factor_corr[nIIV] L;
  vector<lower = 0>[nIIV] omega;
  matrix[nIIV, nSubjects] etaStd;
}

transformed parameters{
  row_vector[nt] cHat;
  vector[nObsPK] cHatObs;
  row_vector[nt] neutHat;
  vector[nObsPD] neutHatObs;
  matrix[8, nt] x;
  array[nSubjects, nTheta] real<lower = 0> parms; // The [1] indicates the parameters are constant
  
  // variables for Matt's trick
  vector<lower = 0>[nIIV] thetaHat;
  matrix<lower = 0>[nSubjects, nIIV] thetaM; 

  // Matt's trick to use unit scale
  thetaHat[1] = CLHat; 
  thetaHat[2] = QHat;
  thetaHat[3] = V1Hat;
  thetaHat[4] = V2Hat;
  thetaHat[5] = mttHat;
  thetaHat[6] = circ0Hat;
  thetaHat[7] = alphaHat;
  thetaM = (rep_matrix(thetaHat, nSubjects) .* 
             exp(diag_pre_multiply(omega, L * etaStd)))';
  
  for(i in 1:nSubjects) {
    parms[i, 1] = thetaM[i, 1] * (weight[i] / 70)^0.75; // CL
    parms[i, 2] = thetaM[i, 2] * (weight[i] / 70)^0.75; // Q
    parms[i, 3] = thetaM[i, 3] * (weight[i] / 70); // V1
    parms[i, 4] = thetaM[i, 4] * (weight[i] / 70); // V2
    parms[i, 5] = kaHat; // ka
    parms[i, 6] = thetaM[i, 5]; // mtt
    parms[i, 7] = thetaM[i, 6]; // circ0
    parms[i, 8] = gamma;
    parms[i, 9] = thetaM[i, 7]; // alpha
  }
                             
  /* group solver */
  x = pmx_solve_group_rk45(twoCptNeutModelODE, 8, len,
                           time, amt, rate, ii, evid, cmt, addl, ss,
                           parms, 
                           1e-6, 1e-6, 500);

  for(i in 1:nSubjects) {
    cHat[start[i]:end[i]] = x[2, start[i]:end[i]] / parms[i, 3]; // divide by V1
    neutHat[start[i]:end[i]] = x[8, start[i]:end[i]] + parms[i, 7]; // Add baseline
  }

  for(i in 1:nObsPK) cHatObs[i] = cHat[iObsPK[i]];
  for(i in 1:nObsPD) neutHatObs[i] = neutHat[iObsPD[i]];
}

model{
  // Priors
  CLHat ~ lognormal(log(10.0), 0.2);
  QHat ~ lognormal(log(15), 0.1);
  V1Hat ~ lognormal(log(35), 0.1);
  V2Hat ~ lognormal(log(100), 0.2);
  kaHat ~ lognormal(log(2.0), 0.2);
  sigma ~ cauchy(0, 0.5);

  mttHat ~ lognormal(log(mttHatPrior), mttHatPriorCV);
  circ0Hat ~ lognormal(log(circ0HatPrior), circ0HatPriorCV);
  alphaHat ~ lognormal(log(alphaHatPrior), alphaHatPriorCV);
  gamma ~ lognormal(log(gammaPrior), gammaPriorCV);
  sigmaNeut ~ cauchy(0, 0.5);

  // Parameters for Matt's trick
  L ~ lkj_corr_cholesky(1);
  to_vector(etaStd) ~ normal(0, 1);
  omega ~ cauchy(0, 0.5);

  // observed data likelihood
  logCObs ~ normal(log(cHatObs), sigma);
  logNeutObs ~ normal(log(neutHatObs), sigmaNeut);
}

generated quantities {
  matrix[8, nt] xPred;
  array[nSubjects, nTheta] real<lower = 0> parmsPred; // [1] indicates the parameters are constant
  row_vector[nt] cHatPred;
  row_vector[nt] neutHatPred;
  vector<lower = 0>[nObsPK] cHatObsCond;
  vector<lower = 0>[nObsPK] cHatObsPred;
  vector<lower = 0>[nObsPD] neutHatObsCond;
  vector<lower = 0>[nObsPD] neutHatObsPred;

  // Variables for IIV  
  matrix[nIIV, nSubjects] etaStdPred;
  matrix<lower = 0>[nSubjects, nIIV] thetaPredM;
  corr_matrix[nIIV] rho;
  
  rho = L * L';
  for(i in 1:nSubjects) {
    for(j in 1:nIIV) {
      etaStdPred[j, i] = normal_rng(0, 1);
    }
  }
  thetaPredM = (rep_matrix(thetaHat, nSubjects) .* 
                exp(diag_pre_multiply(omega, L * etaStdPred)))';
                
  for(i in 1:nSubjects) {
    parmsPred[i, 1] = thetaPredM[i, 1] * (weight[i] / 70)^0.75; // CL
    parmsPred[i, 2] = thetaPredM[i, 2] * (weight[i] / 70)^0.75; // Q
    parmsPred[i, 3] = thetaPredM[i, 3] * (weight[i] / 70); // V1
    parmsPred[i, 4] = thetaPredM[i, 4] * (weight[i] / 70); // V2
    parmsPred[i, 5] = kaHat; // ka
    parmsPred[i, 6] = thetaPredM[i, 5]; // mtt
    parmsPred[i, 7] = thetaPredM[i, 6]; // circ0
    parmsPred[i, 8] = gamma; // gamma
    parmsPred[i, 9] = thetaPredM[i, 7]; // alpha
  }

  xPred = pmx_solve_group_rk45(twoCptNeutModelODE, 8, len,
                         time, amt, rate, ii, evid, cmt, addl, ss,
                         parmsPred, 
                         1e-6, 1e-6, 500);

  for(i in 1:nSubjects) {    
    cHatPred[start[i]:end[i]] = xPred[2, start[i]:end[i]] / parmsPred[i, 3]; // divide by V1
    neutHatPred[start[i]:end[i]] = xPred[8, start[i]:end[i]] + parmsPred[i, 7]; // Add baseline
  }

  // predictions for observed data records
  for(i in 1:nObsPK) cHatObsPred[i] = cHatPred[iObsPK[i]];
  for(i in 1:nObsPD) neutHatObsPred[i] = neutHatPred[iObsPD[i]];
  
  for(i in 1:nObsPK) {
    cHatObsCond[i] = exp(normal_rng(log(fmax(machine_precision(), cHatObs[i])), sigma));
    cHatObsPred[i] = exp(normal_rng(log(fmax(machine_precision(), cHatObsPred[i])), sigma));
  }
  
  for(i in 1:nObsPD) {
    neutHatObsCond[i] = exp(normal_rng(log(fmax(machine_precision(), neutHatObs[i])), sigmaNeut));
    neutHatObsPred[i] = exp(normal_rng(log(fmax(machine_precision(), neutHatObsPred[i])), sigmaNeut));
  }
}
