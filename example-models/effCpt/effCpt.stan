data{
  int<lower = 1> nSubjects;
  int<lower = 1> nt;
  int<lower = 1> nObs;
  array[nObs] int<lower = 1> iObs;
  array[nt] real<lower = 0> amt;
  array[nt] int<lower = 1> cmt;
  array[nt] int<lower = 0> evid;
  array[nSubjects] int<lower = 1> start;
  array[nSubjects] int<lower = 1> end;
  array[nt] real<lower = 0> time;
  vector<lower = 0>[nObs] cObs;
  vector[nObs] respObs;
  array[nSubjects] real<lower = 0> weight;
  array[nt] real<lower = 0> rate;
  array[nt] real<lower = 0> ii;
  array[nt] int<lower = 0> addl;
  array[nt] int<lower = 0> ss;
}

transformed data{
  vector[nObs] logCObs = log(cObs);
  int<lower = 1> nRandom = 5;
  int nCmt = 4;
  array[nCmt] real biovar = rep_array(1.0, nCmt);
  array[nCmt] real tlag = rep_array(0.0, nCmt);
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
  array[nSubjects] vector[nRandom] logtheta_raw;
  array[nSubjects] real logKe0_raw;
  array[nSubjects] real logEC50_raw;
}

transformed parameters{
  vector<lower = 0>[nRandom] thetaHat;
  cov_matrix[nRandom] Omega;
  array[nSubjects] real<lower = 0> CL;
  array[nSubjects] real<lower = 0> Q;
  array[nSubjects] real<lower = 0> V1;
  array[nSubjects] real<lower = 0> V2;
  array[nSubjects] real<lower = 0> ka;
  array[nSubjects] real<lower = 0> ke0;
  array[nSubjects] real<lower = 0> EC50;
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
  array[nSubjects] vector[nRandom] logtheta;
  array[nSubjects] real logKe0;
  array[nSubjects] real logEC50;

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
  array[nSubjects] vector[nRandom] logthetaPred;
  array[nSubjects] real logKe0Pred;
  array[nSubjects] real logEC50Pred;
  row_vector<lower = 0>[nt] cHatPred;
  vector<lower = 0>[nt] cObsCond;
  vector<lower = 0>[nt] cObsPred;
  row_vector<lower = 0>[nt] respHatPred;
  row_vector<lower = 0>[nt] ceHatPred;
  vector[nt] respObsCond;
  vector[nt] respObsPred;
  array[nSubjects] real<lower = 0> CLPred;
  array[nSubjects] real<lower = 0> QPred;
  array[nSubjects] real<lower = 0> V1Pred;
  array[nSubjects] real<lower = 0> V2Pred;
  array[nSubjects] real<lower = 0> kaPred;
  array[nSubjects] real<lower = 0> ke0Pred;
  array[nSubjects] real<lower = 0> EC50Pred;
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
