data{
  int<lower = 1> nSubjects;
  int<lower = 1> nt;
  int<lower = 1> nObs;
  array[nObs] int<lower = 1> iObs;
  array[nt] real<lower = 0> amt;
  array[nt] real<lower = 0> rate;
  array[nt] real<lower = 0> ii;
  array[nt] int<lower = 1> cmt;
  array[nt] int<lower = 0> evid;
  array[nt] int<lower = 0> addl;
  array[nt] int<lower = 0> ss;
  array[nSubjects] int<lower = 1> start;
  array[nSubjects] int<lower = 1> end;
  array[nSubjects] real<lower = 0> weight;
  array[nt] real<lower = 0> time;
  vector<lower = 0>[nObs] cObs;
}

transformed data{
  vector[nObs] logCObs;
  int<lower = 1> nRandom;

  logCObs = log(cObs);

  nRandom = 5;
}

parameters{
  real<lower = 0> CLHat;
  real<lower = 0> QHat;
  real<lower = 0> V1Hat;
  real<lower = 0> V2Hat;
  real<lower = 0> kaHat;
  corr_matrix[nRandom] rho;
  vector<lower = 0>[nRandom] omega;
  real<lower = 0> sigma;
  array[nSubjects] vector[nRandom] logtheta;
}

transformed parameters{
  vector<lower = 0>[nRandom] thetaHat;
  cov_matrix[nRandom] Omega;
  array[nSubjects] real<lower = 0> CL;
  array[nSubjects] real<lower = 0> Q ;
  array[nSubjects] real<lower = 0> V1;
  array[nSubjects] real<lower = 0> V2;
  array[nSubjects] real<lower = 0> ka;
  vector<lower = 0>[nt] cHat;
  vector<lower = 0>[nObs] cHatObs;
  matrix<lower = 0>[3, nt] x;
  array[nSubjects,11] real parms;

  thetaHat[1] = CLHat;
  thetaHat[2] = QHat;
  thetaHat[3] = V1Hat;
  thetaHat[4] = V2Hat;
  thetaHat[5] = kaHat;

  Omega = quad_form_diag(rho, omega); //// diag_matrix(omega) * rho * diag_matrix(omega)

  for(j in 1:nSubjects){
    CL[j] = exp(logtheta[j, 1]) * (weight[j] / 70)^0.75;
    Q[j] = exp(logtheta[j, 2]) * (weight[j] / 70)^0.75;
    V1[j] = exp(logtheta[j, 3]) * weight[j] / 70;
    V2[j] = exp(logtheta[j, 4]) * weight[j] / 70;
    ka[j] = exp(logtheta[j, 5]);

    parms[j, 1] = CL[j];
    parms[j, 2] = Q[j];
    parms[j, 3] = V1[j];
    parms[j, 4] = V2[j];
    parms[j, 5] = ka[j];
    parms[j, 6] = 1; // F1
    parms[j, 7] = 1; // F2
    parms[j, 8] = 1; // F3
    parms[j, 9] = 0; // tlag1
    parms[j, 10] = 0; // tlag2
    parms[j, 11] = 0; // tlag3

    x[, start[j]:end[j]] = pmx_solve_twocpt(time[start[j]:end[j]],
                                            amt[start[j]:end[j]],
                                            rate[start[j]:end[j]],
                                            ii[start[j]:end[j]],
                                            evid[start[j]:end[j]],
                                            cmt[start[j]:end[j]],
                                            addl[start[j]:end[j]],
                                            ss[start[j]:end[j]],
                                            parms[j]);

    for(i in start[j]:end[j])
      cHat[i] = x[2, i] / V1[j];
  }

  cHatObs = cHat[iObs]; //// predictions for observed data records

}

model{
    CLHat ~ normal(0, 20);
    QHat ~ normal(0, 20);
    V1Hat ~ normal(0, 100);
    V2Hat ~ normal(0, 1000);
    kaHat ~ normal(0, 5);
    omega ~ cauchy(0, 2);
    rho ~ lkj_corr(1);
    sigma ~ cauchy(0, 5);

    //// Inter-individual variability
    logtheta ~ multi_normal(log(thetaHat), Omega);

    logCObs ~ normal(log(cHatObs), sigma); //// observed data likelihood
}

generated quantities{
  array[nSubjects] vector[nRandom] logthetaPred;
  vector<lower = 0>[nt] cHatPred;
  array[nt] real cObsCond;
  array[nt] real cObsPred;
  array[nSubjects] real<lower = 0> CLPred;
  array[nSubjects] real<lower = 0> QPred ;
  array[nSubjects] real<lower = 0> V1Pred;
  array[nSubjects] real<lower = 0> V2Pred;
  array[nSubjects] real<lower = 0> kaPred;
  matrix[3, nt] xPred;
  array[nSubjects,11] real parmsPred;

  for(j in 1:nSubjects){
    logthetaPred[j] = multi_normal_rng(log(thetaHat), Omega);
    CLPred[j] = exp(logthetaPred[j, 1]) * (weight[j] / 70)^0.75;
    QPred[j] = exp(logthetaPred[j, 2]) * (weight[j] / 70)^0.75;
    V1Pred[j] = exp(logthetaPred[j, 3]) * weight[j] / 70;
    V2Pred[j] = exp(logthetaPred[j, 4]) * weight[j] / 70;
    kaPred[j] = exp(logthetaPred[j, 5]);

    parmsPred[j, 1] = CLPred[j];
    parmsPred[j, 2] = QPred[j];
    parmsPred[j, 3] = V1Pred[j];
    parmsPred[j, 4] = V2Pred[j];
    parmsPred[j, 5] = kaPred[j];
    parmsPred[j, 6] = 1; // F1
    parmsPred[j, 7] = 1; // F2
    parmsPred[j, 8] = 1; // F3
    parmsPred[j, 9] = 0; // tlag1
    parmsPred[j, 10] = 0; // tlag2
    parmsPred[j, 11] = 0; // tlag3

    xPred[, start[j]:end[j]] = pmx_solve_twocpt(time[start[j]:end[j]],
                                                amt[start[j]:end[j]],
                                                rate[start[j]:end[j]],
                                                ii[start[j]:end[j]],
                                                evid[start[j]:end[j]],
                                                cmt[start[j]:end[j]],
                                                addl[start[j]:end[j]],
                                                ss[start[j]:end[j]],
                                                parmsPred[j]);

    for(i in start[j]:end[j])
      cHatPred[i] = xPred[2, i] / V1Pred[j];
  }

  for(i in 1:nt){
    if(time[i] == 0){
      cObsCond[i] = -99;
      cObsPred[i] = -99;
    }else{
      cObsCond[i] = exp(normal_rng(log(cHat[i]), sigma));
      cObsPred[i] = exp(normal_rng(log(cHatPred[i]), sigma));
    }
  }

}
