data{
  int<lower = 1> nId;
  int<lower = 1> nt;
  int<lower = 1> nObsPK;
  int<lower = 1> nObsPD;
  int<lower = 1> iObsPK[nObsPK];
  int<lower = 1> iObsPD[nObsPD];
  real<lower = 0> amt[nt];
  int<lower = 1> cmt[nt];
  int<lower = 0> evid[nt];
  real<lower = 0> ii[nt];
  int<lower = 0> addl[nt];
  real<lower = 0> rate[nt];
  int<lower = 0> ss[nt];
  int<lower = 1> start[nId];
  int<lower = 1> end[nId];
  real<lower = 0> weight[nId];
  real<lower = 0> time[nt];
  vector<lower = 0>[nObsPK] cObs;
  vector<lower = 0>[nObsPD] respObs;
  real<lower = 0> CLPrior;
  real<lower = 0> VPrior;
  real<lower = 0> kaPrior;
  real<lower = 0> CLPriorCV;
  real<lower = 0> VPriorCV;
  real<lower = 0> kaPriorCV;
  real<lower = 0> ke0Prior;
  real<lower = 0> ke0PriorCV;
  real<lower = 0> E0Prior;
  real<lower = 0> E0PriorCV;
  real<lower = 0> EmaxPrior;
  real<lower = 0> EmaxPriorCV;
  real<lower = 0> EC50Prior;
  real<lower = 0> EC50PriorCV;
  real<lower = 0> gammaPrior;
  real<lower = 0> gammaPriorCV;
}

transformed data{
  vector[nObsPK] logCObs = log(cObs);
  vector[nObsPD] logRespObs = log(respObs);
  int<lower = 1> nRandom = 6;
  int<lower = 1> nCmt = 3;
  real F[nCmt] = rep_array(1.0, nCmt);
  real tLag[nCmt] = rep_array(0.0, nCmt);
}

parameters{
  real<lower = 0> CLHat;
  real<lower = 0> VHat;
  real<lower = CLHat / VHat> kaHat;
  real<lower = 0> ke0Hat;
  real<lower = 0> E0Hat;
  real<lower = 0> EmaxHat;
  real<lower = 0> EC50;
  real<lower = 0> gamma;
  vector<lower = 0>[nRandom] omega;
  cholesky_factor_corr[nRandom] L;
  real<lower = 0> sigma;
  real<lower = 0> sigmaPD;
  matrix[nRandom, nId] etaStd;
}

transformed parameters{
  vector<lower = 0>[nRandom] thetaHat;
  matrix<lower = 0>[nId, nRandom] theta;
  real<lower = 0> CL[nId];
  real<lower = 0> V[nId];
  real<lower = 0> ka[nId];
  real<lower = 0> ke0[nId];
  real<lower = 0> Emax[nId];
  real<lower = 0> E0[nId];
  vector<lower = 0>[nt] cHat;
  vector<lower = 0>[nt] ceHat;
  vector<lower = 0>[nt] respHat;
  vector<lower = 0>[nObsPK] cHatObs;
  vector<lower = 0>[nObsPD] respHatObs;
  matrix[nCmt, nCmt] K;
  matrix[nCmt, nt] x;

  thetaHat = to_vector({CLHat, VHat, kaHat, ke0Hat, EmaxHat, E0Hat});

  theta = (rep_matrix(thetaHat, nId) .* 
	   exp(diag_pre_multiply(omega, L * etaStd)))';

  for(j in 1:nId){
    CL[j] = theta[j, 1] * (weight[j] / 70)^0.75;
    V[j] = theta[j, 2] * weight[j] / 70;
    ka[j] = theta[j, 3];
    ke0[j] = theta[j, 4];
    Emax[j] = theta[j, 5];
    E0[j] = theta[j, 6];

    K = rep_matrix(0, nCmt, nCmt);
    
    K[1, 1] = -ka[j];
    K[2, 1] = ka[j];
    K[2, 2] = -CL[j] / V[j];
    K[3, 2] = ke0[j];
    K[3, 3] = -ke0[j];

    x[, start[j]:end[j] ] = pmx_solve_linode(time[start[j]:end[j]], 
                                             amt[start[j]:end[j]],
                                             rate[start[j]:end[j]],
                                             ii[start[j]:end[j]],
                                             evid[start[j]:end[j]],
                                             cmt[start[j]:end[j]],
                                             addl[start[j]:end[j]],
                                             ss[start[j]:end[j]],
                                             K, F, tLag);

    cHat[start[j]:end[j]] = x[2, start[j]:end[j]]' / V[j];
    ceHat[start[j]:end[j]] = x[3, start[j]:end[j]]' / V[j];
    for(i in start[j]:end[j])
      respHat[i] = E0[j] + Emax[j] * ceHat[i]^gamma / (EC50^gamma + ceHat[i]^gamma);
  }

  cHatObs = cHat[iObsPK]; // predictions for observed data records
  respHatObs = respHat[iObsPD]; // predictions for observed data records
}

model{
  CLHat ~ lognormal(log(CLPrior), CLPriorCV);
  VHat ~ lognormal(log(VPrior), VPriorCV);
  kaHat ~ lognormal(log(kaPrior), kaPriorCV);
  sigma ~ cauchy(0, 1);
  ke0Hat ~ lognormal(log(ke0Prior), ke0PriorCV);
  E0Hat ~ lognormal(log(E0Prior), E0PriorCV);
  EmaxHat ~ lognormal(log(EmaxPrior), EmaxPriorCV);
  EC50 ~ lognormal(log(EC50Prior), EC50PriorCV);
  gamma ~ lognormal(log(gammaPrior), gammaPriorCV);
  sigmaPD ~ cauchy(0, 1);
  omega ~ cauchy(0, 1);
  L ~ lkj_corr_cholesky(1);

  // Inter-individual variability
  to_vector(etaStd) ~ normal(0, 1);

  logCObs ~ normal(log(cHatObs), sigma); // observed data likelihood
  logRespObs ~ normal(log(respHatObs), sigmaPD);
}

generated quantities{
  real<lower = 0> CLPred[nId];
  real<lower = 0> VPred[nId];
  real<lower = 0> kaPred[nId];
  real<lower = 0> ke0Pred[nId];
  real<lower = 0> EmaxPred[nId];
  real<lower = 0> E0Pred[nId];
  vector[nt] cHatPred;
  vector[nt] ceHatPred;
  vector[nt] respHatPred;
  real cObsCond[nt];
  real respObsCond[nt];
  real cObsPred[nt];
  real respObsPred[nt];
  matrix<lower = 0>[nId, nRandom] thetaPred;
  matrix[nCmt, nt] xPred;
  matrix[nCmt, nCmt] KPred;
  corr_matrix[nRandom] rho;
  matrix[nRandom, nId] etaStdPred;

  rho = multiply_lower_tri_self_transpose(L);
  for(j in 1:nId) 
    for(i in 1:nRandom)
      etaStdPred[i, j] = normal_rng(0, 1);

  thetaPred = (rep_matrix(thetaHat, nId) .* 
	       exp(diag_pre_multiply(omega, L * etaStdPred)))';

  for(j in 1:nId){
    CLPred[j] = thetaPred[j, 1] * (weight[j] / 70)^0.75;
    VPred[j] = thetaPred[j, 2] * weight[j] / 70;
    kaPred[j] = thetaPred[j, 3];
    ke0Pred[j] = thetaPred[j, 4];
    EmaxPred[j] = thetaPred[j, 5];
    E0Pred[j] = thetaPred[j, 6];

    KPred = rep_matrix(0, nCmt, nCmt);
    
    KPred[1, 1] = -kaPred[j];
    KPred[2, 1] = kaPred[j];
    KPred[2, 2] = -CLPred[j] / VPred[j];
    KPred[3, 2] = ke0Pred[j];
    KPred[3, 3] = -ke0Pred[j];

    xPred[, start[j]:end[j] ] = pmx_solve_linode(time[start[j]:end[j]], 
                                                 amt[start[j]:end[j]],
                                                 rate[start[j]:end[j]],
                                                 ii[start[j]:end[j]],
                                                 evid[start[j]:end[j]],
                                                 cmt[start[j]:end[j]],
                                                 addl[start[j]:end[j]],
                                                 ss[start[j]:end[j]],
                                                 KPred, F, tLag);

    cHatPred[start[j]:end[j]] = xPred[2, start[j]:end[j]]' / VPred[j];
    ceHatPred[start[j]:end[j]] = xPred[3, start[j]:end[j]]' / VPred[j];
    for(i in start[j]:end[j])
      respHatPred[i] = E0Pred[j] + EmaxPred[j] * 
        ceHatPred[i]^gamma / (EC50^gamma + ceHatPred[i]^gamma);
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
    respObsCond[i] = exp(normal_rng(log(fmax(machine_precision(),
      respHat[i])), sigmaPD));
    respObsPred[i] = exp(normal_rng(log(fmax(machine_precision(),
      respHatPred[i])), sigmaPD));
  }
}
