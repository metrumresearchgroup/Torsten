data{
  int<lower = 1> nSubjects;
  int<lower = 1> nt;
  int<lower = 1> nObs;
  int<lower = 1> iObs[nObs];
  real<lower = 0> amt[nt];
  int<lower = 1> cmt[nt];
  int<lower = 0> evid[nt];
  int<lower = 1> start[nSubjects];
  int<lower = 1> end[nSubjects];
  real<lower = 0> time[nt];
  vector<lower = 0>[nObs] cObs;
  vector[nObs] respObs;
  real<lower = 0> weight[nSubjects];
  real<lower = 0> rate[nt];
  real<lower = 0> ii[nt];
  int<lower = 0> addl[nt];
  int<lower = 0> ss[nt];
}

transformed data{
  vector[nObs] logCObs = log(cObs);
  int<lower = 1> nRandom = 5;
  int nCmt = 4;
  real biovar[nCmt] = rep_array(1.0, nCmt);
  real tlag[nCmt] = rep_array(0.0, nCmt);
}

parameters{
  real<lower = 0> CLHat;
  real<lower = 0> QHat;
  real<lower = 0> V1Hat;
  real<lower = 0> V2Hat;
  //  real<lower = 0> kaHat;
  real<lower = (CLHat / V1Hat + QHat / V1Hat + QHat / V2Hat +
		sqrt((CLHat / V1Hat + QHat / V1Hat + QHat / V2Hat)^2 -
		     4 * CLHat / V1Hat * QHat / V2Hat)) / 2> kaHat; // ka > lambda_1
  real<lower = 0> ke0Hat;
  real<lower = 0> EC50Hat;
  vector<lower = 0>[nRandom] omega;
  corr_matrix[nRandom] rho;
  real<lower = 0> omegaKe0;
  real<lower = 0> omegaEC50;
  real<lower = 0> sigma;
  real<lower = 0> sigmaResp;

  // reparameterization
  vector[nRandom] logtheta_raw[nSubjects];
  real logKe0_raw[nSubjects];
  real logEC50_raw[nSubjects];
}

transformed parameters{
  vector<lower = 0>[nRandom] thetaHat;
  cov_matrix[nRandom] Omega;
  real<lower = 0> CL[nSubjects];
  real<lower = 0> Q[nSubjects];
  real<lower = 0> V1[nSubjects];
  real<lower = 0> V2[nSubjects];
  real<lower = 0> ka[nSubjects];
  real<lower = 0> ke0[nSubjects];
  real<lower = 0> EC50[nSubjects];
  matrix[nCmt, nCmt] K;
  real k10;
  real k12;
  real k21;
  row_vector<lower = 0>[nt] cHat;
  row_vector<lower = 0>[nObs] cHatObs;
  row_vector<lower = 0>[nt] respHat;
  row_vector<lower = 0>[nObs] respHatObs;
  row_vector<lower = 0>[nt] ceHat;
  matrix[nCmt, nt] x;
  
  matrix[nRandom, nRandom] L;
  vector[nRandom] logtheta[nSubjects];
  real logKe0[nSubjects];
  real logEC50[nSubjects];

  thetaHat[1] = CLHat;
  thetaHat[2] = QHat;
  thetaHat[3] = V1Hat;
  thetaHat[4] = V2Hat;
  thetaHat[5] = kaHat;

  Omega = quad_form_diag(rho, omega); // diag_matrix(omega) * rho * diag_matrix(omega)
  L = cholesky_decompose(Omega);

  for(j in 1:nSubjects){
    logtheta[j] = log(thetaHat) + L * logtheta_raw[j];
    logKe0[j] = log(ke0Hat) + logKe0_raw[j] * omegaKe0;
    logEC50[j] = log(EC50Hat) + logEC50_raw[j] * omegaEC50;

    CL[j] = exp(logtheta[j, 1]) * (weight[j] / 70)^0.75;
    Q[j] = exp(logtheta[j, 2]) * (weight[j] / 70)^0.75;
    V1[j] = exp(logtheta[j, 3]) * weight[j] / 70;
    V2[j] = exp(logtheta[j, 4]) * weight[j] / 70;
    ka[j] = exp(logtheta[j, 5]);
    ke0[j] = exp(logKe0[j]);
    EC50[j] = exp(logEC50[j]);
    
    k10 = CL[j] / V1[j];
    k12 = Q[j] / V1[j];
    k21 = Q[j] / V2[j];

    K = rep_matrix(0, nCmt, nCmt);
    
    K[1, 1] = -ka[j];
    K[2, 1] = ka[j];
    K[2, 2] = -(k10 + k12);
    K[2, 3] = k21;
    K[3, 2] = k12;
    K[3, 3] = -k21;
    K[4, 2] = ke0[j];
    K[4, 4] = -ke0[j];
    
    x[, start[j]:end[j] ] = pmx_solve_linode(time[start[j]:end[j]], amt[start[j]:end[j]],
                                             rate[start[j]:end[j]], ii[start[j]:end[j]],
                                             evid[start[j]:end[j]], cmt[start[j]:end[j]],
                                             addl[start[j]:end[j]], ss[start[j]:end[j]], K, biovar, tlag);

    cHat[start[j]:end[j]] = 1000 * x[2, start[j]:end[j]] ./ V1[j];
    ceHat[start[j]:end[j]] = 1000 * x[4, start[j]:end[j]] ./ V1[j];
    respHat[start[j]:end[j]] = 100 * ceHat[start[j]:end[j]] ./ 
       (EC50[j] + ceHat[start[j]:end[j]]);
  }

  cHatObs = cHat[iObs];
  respHatObs = respHat[iObs];
}

model{
  // Prior
  CLHat ~ lognormal(log(10), 0.2);
  QHat ~ lognormal(log(15), 0.2);
  V1Hat ~ lognormal(log(30), 0.2);
  V2Hat ~ lognormal(log(100), 0.2);
  kaHat ~ lognormal(log(5), 0.25);
  ke0Hat ~ lognormal(log(10), 0.25);
  EC50Hat ~ lognormal(log(1.0), 0.2);
  omega ~ normal(0, 0.2);
  rho ~ lkj_corr(1); 
  omegaKe0 ~ normal(0, 0.2);
  omegaEC50 ~ normal(0, 0.2);
  sigma ~ cauchy(0, 0.2);
  sigmaResp ~ cauchy(0, 0.2);

  // Inter-individual variability
  for (i in 1:nSubjects) {
    logtheta_raw[i] ~ std_normal();      
  }
  logKe0_raw ~ std_normal();
  logEC50_raw ~ std_normal();

  // Likelihood
  logCObs ~ normal(log(cHatObs), sigma); 
  respObs ~ normal(respHatObs, sigmaResp); 
}

generated quantities{
  vector[nRandom] logthetaPred[nSubjects];
  real logKe0Pred[nSubjects];
  real logEC50Pred[nSubjects];
  row_vector<lower = 0>[nt] cHatPred;
  vector<lower = 0>[nt] cObsCond;
  vector<lower = 0>[nt] cObsPred;
  row_vector<lower = 0>[nt] respHatPred;
  row_vector<lower = 0>[nt] ceHatPred;
  vector[nt] respObsCond;
  vector[nt] respObsPred;
  real<lower = 0> CLPred[nSubjects];
  real<lower = 0> QPred[nSubjects];
  real<lower = 0> V1Pred[nSubjects];
  real<lower = 0> V2Pred[nSubjects];
  real<lower = 0> kaPred[nSubjects];
  real<lower = 0> ke0Pred[nSubjects];
  real<lower = 0> EC50Pred[nSubjects];
  matrix[nCmt, nt] xPred;
  matrix[nCmt, nCmt] KPred;
  real k10Pred;
  real k12Pred;
  real k21Pred;

  for(j in 1:nSubjects){
    logthetaPred[j] = multi_normal_rng(log(thetaHat), Omega);
    logKe0Pred[j] = normal_rng(log(ke0Hat), omegaKe0);
    logEC50Pred[j] = normal_rng(log(EC50Hat), omegaEC50);
    CLPred[j] = exp(logthetaPred[j, 1]) * (weight[j] / 70)^0.75;
    QPred[j] = exp(logthetaPred[j, 2]) * (weight[j] / 70)^0.75;
    V1Pred[j] = exp(logthetaPred[j, 3]) * weight[j] / 70;
    V2Pred[j] = exp(logthetaPred[j, 4]) * weight[j] / 70;
    kaPred[j] = exp(logthetaPred[j, 5]);
    ke0Pred[j] = exp(logKe0Pred[j]);
    EC50Pred[j] = exp(logEC50Pred[j]);
    
    k10Pred = CLPred[j] / V1Pred[j];
    k12Pred = QPred[j] / V1Pred[j];
    k21Pred = QPred[j] / V2Pred[j];

    KPred = rep_matrix(0, nCmt, nCmt);
    
    KPred[1, 1] = -kaPred[j];
    KPred[2, 1] = kaPred[j];
    KPred[2, 2] = -(k10Pred + k12Pred);
    KPred[2, 3] = k21Pred;
    KPred[3, 2] = k12Pred;
    KPred[3, 3] = -k21Pred;
    KPred[4, 2] = ke0Pred[j];
    KPred[4, 4] = -ke0Pred[j];
    
    xPred[, start[j]:end[j] ] = pmx_solve_linode(time[start[j]:end[j]], amt[start[j]:end[j]],
                                                 rate[start[j]:end[j]], ii[start[j]:end[j]],
                                                 evid[start[j]:end[j]], cmt[start[j]:end[j]],
                                                 addl[start[j]:end[j]], ss[start[j]:end[j]], KPred, biovar, tlag);

    cHatPred[start[j]:end[j]] = 1000 * xPred[2, start[j]:end[j]] ./ V1Pred[j];
    ceHatPred[start[j]:end[j]] = 1000 * xPred[4, start[j]:end[j]] ./ V1Pred[j];
    respHatPred[start[j]:end[j]] = 100 * ceHatPred[start[j]:end[j]] ./
      (EC50Pred[j] + ceHatPred[start[j]:end[j]]);
  }
  
  for(i in 1:nt){
    if(time[i] == 0){
      cObsCond[i] = 0;
      cObsPred[i] = 0;
      respObsCond[i] = 0;
      respObsPred[i] = 0;
    }else{
      cObsCond[i] = exp(normal_rng(log(cHat[i]), sigma));
      cObsPred[i] = exp(normal_rng(log(cHatPred[i]), sigma));
      respObsCond[i] = normal_rng(respHat[i], sigmaResp);
      respObsPred[i] = normal_rng(respHatPred[i], sigmaResp);
    }
  }
}
