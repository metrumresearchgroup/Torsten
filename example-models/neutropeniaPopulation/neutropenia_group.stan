functions{
  vector twoCptNeutModelODE(real t, vector x, array[] real parms,
                            array[] real rdummy, array[] int idummy){
    real CL = parms[1];
    real Q = parms[2];
    real V1 = parms[3];
    real V2 = parms[4];
    real ka = parms[5];
    real mtt = parms[6];
    real circ0 = parms[7];
    real gamma = parms[8];
    real alpha = parms[9];
    real k10 = CL / V1;
    real k12 = Q / V1;
    real k21 = Q / V2;
    real ktr = 4 / mtt;
    vector[8] dxdt;
    real conc = x[2]/V1;
    real EDrug = fmin(1.0, alpha * conc);
    real prol = x[4] + circ0;
    real transit1 = x[5] + circ0;
    real transit2 = x[6] + circ0;
    real transit3 = x[7] + circ0;
    real circ = fmax(machine_precision(), x[8] + circ0);

    //    print("parms ", parms);
    dxdt[1] = -ka * x[1];
    dxdt[2] = ka * x[1] - (k10 + k12) * x[2] + k21 * x[3];
    dxdt[3] = k12 * x[2] - k21 * x[3];
    // x[4], x[5], x[6], x[7] and x[8] are differences from circ0.
    // This is a device for implementing a modeled initial condition
    dxdt[4] = ktr * prol * ((1 - EDrug) * ((circ0 / circ)^gamma) - 1);
    dxdt[5] = ktr * (prol - transit1);
    dxdt[6] = ktr * (transit1 - transit2);
    dxdt[7] = ktr * (transit2 - transit3);
    dxdt[8] = ktr * (transit3 - circ);

    //    print("dxdt = ", dxdt);
    return dxdt;
  }
}

data{
  int<lower = 1> nId;
  int<lower = 1> nt;
  int<lower = 1> nObsPK;
  int<lower = 1> nObsPD;
  array[nObsPK] int<lower = 1> iObsPK;
  array[nObsPD] int<lower = 1> iObsPD;
  array[nt] real<lower = 0> amt;
  array[nt] int<lower = 1> cmt;
  array[nt] int<lower = 0> evid;
  array[nId] int<lower = 1> start;
  array[nId] int<lower = 1> end;
  array[nId] real<lower = 0> weight;
  array[nt] real<lower = 0> time;
  vector<lower = 0>[nObsPK] cObs;
  vector<lower = 0>[nObsPD] neutObs;
  real<lower = 0> CLPrior;
  real<lower = 0> QPrior;
  real<lower = 0> V1Prior;
  real<lower = 0> V2Prior;
  real<lower = 0> kaPrior;
  real<lower = 0> CLPriorCV;
  real<lower = 0> QPriorCV;
  real<lower = 0> V1PriorCV;
  real<lower = 0> V2PriorCV;
  real<lower = 0> kaPriorCV;
  real<lower = 0> circ0Prior;
  real<lower = 0> circ0PriorCV;
  real<lower = 0> mttPrior;
  real<lower = 0> mttPriorCV;
  real<lower = 0> gammaPrior;
  real<lower = 0> gammaPriorCV;
  real<lower = 0> alphaPrior;
  real<lower = 0> alphaPriorCV;
}

transformed data{
  array[nt] real<lower = 0> rate=rep_array(0.0, nt);
  array[nt] real<lower = 0> ii=rep_array(0.0, nt);
  array[nt] int<lower = 0> addl=rep_array(0, nt);
  array[nt] int<lower = 0> ss=rep_array(0, nt);
  vector[nObsPK] logCObs = log(cObs);
  vector[nObsPD] logNeutObs = log(neutObs);
  int<lower = 1> nRandom = 8;
  int<lower = 1> nCmt = 8;
  array[nId] int len;

  for (j in 1:nId) {
    len[j] = end[j] - start[j] + 1;
  }
}

parameters{
  real<lower = 0> CLHat;
  real<lower = 0> QHat;
  real<lower = 0> V1Hat;
  real<lower = 0> V2Hat;
  real<lower = 0> kaHat;
  real<lower = 0> mttHat;
  real<lower = 0> circ0Hat;
  real<lower = 0> alphaHat;
  real<lower = 0, upper = 1> gamma;
  vector<lower = 0, upper = 1>[nRandom] omega;
  cholesky_factor_corr[nRandom] L;
  real<lower = 0> sigma;
  real<lower = 0> sigmaNeut;
  matrix[nRandom, nId] etaStd;
}

transformed parameters{
  vector<lower = 0>[nRandom] thetaHat =
    to_vector({CLHat, QHat, V1Hat, V2Hat, kaHat, mttHat, circ0Hat, alphaHat});
  array[nId] real<lower = 0> CL;
  array[nId] real<lower = 0> Q;
  array[nId] real<lower = 0> V1;
  array[nId] real<lower = 0> V2;
  array[nId] real<lower = 0> ka;
  array[nId] real<lower = 0> mtt;
  array[nId] real<lower = 0> circ0;
  array[nId] real<lower = 0> alpha;
  matrix<lower = 0>[nId, nRandom] theta;

  vector[nt] cHat;
  vector[nObsPK] cHatObs;
  vector[nt] neutHat;
  vector[nObsPD] neutHatObs;
  matrix[8, nt] x;
  array[nId, 9] real<lower = 0> parms;

  theta = (rep_matrix(thetaHat, nId) .* 
	   exp(diag_pre_multiply(omega, L * etaStd)))';

  for(j in 1:nId){
    CL[j] = theta[j, 1] * (weight[j] / 70)^0.75;
    Q[j] = theta[j, 2] * (weight[j] / 70)^0.75;
    V1[j] = theta[j, 3] * weight[j] / 70;
    V2[j] = theta[j, 4] * weight[j] / 70;
    ka[j] = theta[j, 5];
    mtt[j] = theta[j, 6];
    circ0[j] = theta[j, 7];
    alpha[j] = theta[j, 8];

    parms[j, ] = {CL[j], Q[j], V1[j], V2[j], ka[j], mtt[j], circ0[j], gamma, alpha[j]};
  }

  x = pmx_solve_group_bdf(twoCptNeutModelODE, 8, len,
                          time, amt, rate, ii, evid, cmt, addl, ss,
                          parms, 1e-6, 1e-6, 100000);

  for(j in 1:nId) {
    cHat[start[j]:end[j]] = x[2, start[j]:end[j]]' / V1[j];
    neutHat[start[j]:end[j]] = x[8, start[j]:end[j]]' + circ0[j];
  }

  cHatObs = cHat[iObsPK]; // predictions for observed data records  }
  neutHatObs = neutHat[iObsPD]; // predictions for observed data records
}

model{
  CLHat ~ lognormal(log(CLPrior), CLPriorCV);
  QHat ~ lognormal(log(QPrior), QPriorCV);
  V1Hat ~ lognormal(log(V1Prior), V1PriorCV);
  V2Hat ~ lognormal(log(V2Prior), V2PriorCV);
  kaHat ~ lognormal(log(kaPrior), kaPriorCV);
  sigma ~ cauchy(0, 1);
  mttHat ~ lognormal(log(mttPrior), mttPriorCV);
  circ0Hat ~ lognormal(log(circ0Prior), circ0PriorCV);
  alphaHat ~ lognormal(log(alphaPrior), alphaPriorCV);
  gamma ~ lognormal(log(gammaPrior), gammaPriorCV);
  sigmaNeut ~ cauchy(0, 1);
  omega ~ normal(0, 0.5);
  L ~ lkj_corr_cholesky(1);

  // Inter-individual variability
  to_vector(etaStd) ~ normal(0, 1);

  logCObs ~ normal(log(cHatObs), sigma); // observed data likelihood
  logNeutObs ~ normal(log(neutHatObs), sigmaNeut);
}

generated quantities {
  array[nId] real<lower = 0> CLPred;
  array[nId] real<lower = 0> QPred;
  array[nId] real<lower = 0> V1Pred;
  array[nId] real<lower = 0> V2Pred;
  array[nId] real<lower = 0> kaPred;
  array[nId] real<lower = 0> mttPred;
  array[nId] real<lower = 0> circ0Pred;
  array[nId] real<lower = 0> alphaPred;
  vector[nt] cHatPred;
  vector[nt] neutHatPred;
  array[nt] real cObsCond;
  array[nt] real neutObsCond;
  array[nt] real cObsPred;
  array[nt] real neutObsPred;
  matrix<lower = 0>[nId, nRandom] thetaPred;
  matrix[nCmt, nt] xPred;
  array[nId, 9] real<lower = 0> parmsPred;
  corr_matrix[nRandom] rho;
  matrix[nRandom, nId] etaStdPred;

  rho = L * L';
  for(j in 1:nId) 
    for(i in 1:nRandom)
      etaStdPred[i, j] = normal_rng(0, 1);

  thetaPred = (rep_matrix(thetaHat, nId) .* exp(diag_pre_multiply(omega, L * etaStdPred)))';

  for(j in 1:nId){
    CLPred[j] = thetaPred[j, 1] * (weight[j] / 70)^0.75;
    QPred[j] = thetaPred[j, 2] * (weight[j] / 70)^0.75;
    V1Pred[j] = thetaPred[j, 3] * weight[j] / 70;
    V2Pred[j] = thetaPred[j, 4] * weight[j] / 70;
    kaPred[j] = thetaPred[j, 5];
    mttPred[j] = thetaPred[j, 6];
    circ0Pred[j] = thetaPred[j, 7];
    alphaPred[j] = thetaPred[j, 8];

    parmsPred[j, ] = {CLPred[j], QPred[j], V1Pred[j], V2Pred[j], kaPred[j], mttPred[j],
      circ0Pred[j], gamma, alphaPred[j]};
  }

  xPred = pmx_solve_group_bdf(twoCptNeutModelODE, 8, len, time, amt, rate, ii, evid, cmt, addl, ss, parmsPred, 1e-6, 1e-6, 1e8);

  for(j in 1:nId){
    cHatPred[start[j]:end[j]] = xPred[2, start[j]:end[j]]' / V1Pred[j];
    neutHatPred[start[j]:end[j]] = xPred[8, start[j]:end[j]]' + circ0Pred[j];
  }

  for(i in 1:nt){
    if(time[i] == 0){
      cObsCond[i] = 0;
      cObsPred[i] = 0;
    }else{
      cObsCond[i] = exp(normal_rng(log(fmax(machine_precision(), cHat[i])),
				   sigma));
      cObsPred[i] = exp(normal_rng(log(fmax(machine_precision(),
					    cHatPred[i])), sigma));
    }
    neutObsCond[i] = exp(normal_rng(log(fmax(machine_precision(),
					     neutHat[i])), sigmaNeut));
    neutObsPred[i] = exp(normal_rng(log(fmax(machine_precision(),
					     neutHatPred[i])), sigmaNeut));
  }
}
