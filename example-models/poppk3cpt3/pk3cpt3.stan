functions{
  real normal_inv(real p, real mu, real sigma){
    return inv_Phi(p) * sigma + mu;
  }

  real normal_inv_trunc(real p, real mu, real sigma, real lo, real up){
    return normal_inv(p * normal_cdf(up | mu, sigma) + 
                (1 - p) * normal_cdf(lo | mu, sigma), mu, sigma);
  }
  
  real normal_trunc_rng(real mu, real sigma, real lo, real up){
    return normal_inv_trunc(uniform_rng(0, 1), mu, sigma, lo, up);
  }
}

data{
  // General data items
  int<lower = 1> nID;
  int<lower = 1> nt;
  int<lower = 1> nObs;
  array[nObs] int<lower = 1> iObs;
  array[nt] real<lower = 0> amt;
  array[nt] int<lower = 1> cmt;
  array[nt] int<lower = 0> evid;
  array[nt] real<lower = 0> rate;
  array[nt] real<lower = 0> ii;
  array[nt] int<lower = 0> addl;
  array[nt] int<lower = 0> ss;
  array[nID] real<lower = 0> weight;
  array[nID] int<lower = 1> start;
  array[nID] int<lower = 1> end;
  array[nt] real<lower = 0> time;
  vector<lower = 0>[nObs] cObs;
}

transformed data{

  // Integers required to specify dimensions
  int<lower = 1> nRandom = 7; // Number of random effects
  int<lower = 1> nCmt = 4; // Number of model compartments
  int<lower = 1> nParms = 7; // Number of parameters passed to Torsten function

  // Fixed value parameters, e.g.,
  array[nCmt] real F = rep_array(1.0, nCmt);
  array[nCmt] real tLag = rep_array(0.0, nCmt);
}

parameters{
  // Population-level model parameters
  // These are the parameters for which you specify prior distributions
  // and initial estimates, e.g., 
  real<lower = 0> CLHat;
  real<lower = 10> V2Hat;
  real<lower = 0> Q3Hat;
  real<lower = 0> V3Hat;
  real<lower = 0> Q4Hat;
  real<lower = 0> V4Hat;
  real<lower = 0> kaHat;
  real<lower = 0> gammaCL;
  real<lower = 0> gammaV2;
  
  //  corr_matrix[nRandom] rho;
  cholesky_factor_corr[2] L; // correlation only between CL & V2
  real<lower = 0> omegaCL;
  real<lower = 0> omegaV2;
  array[nRandom - 2] real<lower = 0> omegasq;
  real<lower = 0> sigma;
  
  // Individual-level model parameters directly sampled from the IIV
  // distribution
  //  vector[nRandom] logtheta[nID];
  matrix[nRandom, nID] eta;
}

transformed parameters{
  // Vector of PK parameter typical values -- only those with IIV
  vector<lower = 0>[nRandom] thetaHat = [CLHat, V2Hat, Q3Hat, V3Hat, Q4Hat, V4Hat, kaHat]';

  // Matrix of individual-level model parameters
  matrix<lower = 0>[nID, nRandom] theta;

  // Individual-level model parameters with recognizable names, e.g.,
  array[nID] real<lower = 0> CL;
  array[nID] real<lower = 0> V2;
  array[nID] real<lower = 0> Q3;
  array[nID] real<lower = 0> V3;
  array[nID] real<lower = 0> Q4;
  array[nID] real<lower = 0> V4;
  array[nID] real<lower = 0> ka;
  
  vector<lower = 0>[nRandom] omega;

  // Covariance matrix
  //  cov_matrix[nRandom] Omega;

  // Predicted concentrations (without residual variation)
  vector<lower = 0>[nt] cHat; // All events

  // Amounts in each compartment at each event
  matrix[nCmt, nt] x;

  // Matrix used to pass parameters to the Torsten function
  matrix[nCmt, nCmt] K;

  omega[1] = omegaCL; // sd(log(CL))
  omega[2] = omegaV2; // sd(log(V2))
  omega[3] = sqrt(omegasq[1]); // sd((log(Q3)))
  omega[4] = sqrt(omegasq[2]); // sd((log(V3)))
  omega[5] = sqrt(omegasq[3]); // sd((log(Q4)))
  omega[6] = sqrt(omegasq[4]); // sd((log(V4)))
  omega[7] = sqrt(omegasq[5]); // sd((log(ka)))

  //  Omega = quad_form_diag(rho, omega); // diag_matrix(omega) * rho * diag_matrix(omega)
  theta[, 1:2] = (rep_matrix(thetaHat[1:2], nID) .* 
          exp(diag_pre_multiply(omega[1:2], L * eta[1:2,])))';
  for(i in 3:nRandom)
    theta[, i] = thetaHat[i] * exp(omega[i] * eta[i,])';

  for(j in 1:nID){
    
    // Calculation of individual parameter values given theta and covariates, e.g.
    CL[j] = theta[j, 1] * (weight[j] / 70)^gammaCL;
    V2[j] = theta[j, 2] * (weight[j] / 70)^gammaV2;
    Q3[j] = theta[j, 3] * (weight[j] / 70)^0.75;
    V3[j] = theta[j, 4] * (weight[j] / 70);
    Q4[j] = theta[j, 5] * (weight[j] / 70)^0.75;
    V4[j] = theta[j, 6] * (weight[j] / 70);
    ka[j] = theta[j, 7];

    // Pack individual PK parameters into K matrix, e.g.
    
    K = rep_matrix(0, nCmt, nCmt);
    
    K[1, 1] = -ka[j];
    K[2, 1] = ka[j];
    K[2, 2] = -(CL[j] + Q3[j] + Q4[j]) / V2[j];
    K[2, 3] = Q3[j] / V3[j];
    K[2, 4] = Q4[j] / V4[j];
    K[3, 2] = Q3[j] / V2[j];
    K[3, 3] = -Q3[j] / V3[j];
    K[4, 2] = Q4[j] / V2[j];
    K[4, 4] = -Q4[j] / V4[j];

    x[, start[j]:end[j]] = pmx_solve_linode(time[start[j]:end[j]], 
                    amt[start[j]:end[j]],
                    rate[start[j]:end[j]],
                    ii[start[j]:end[j]],
                    evid[start[j]:end[j]],
                    cmt[start[j]:end[j]],
                    addl[start[j]:end[j]],
                    ss[start[j]:end[j]],
                    K, F, tLag);

    // Calculate target concentration for specified compartment.
    // Change compartment number and distribution volume as appropriate.

    cHat[start[j]:end[j]] = x[2, start[j]:end[j]]' ./ V2[j];
  }

}

model{
  // Priors
  CLHat ~ lognormal(log(138), 1); // normal(0, 100);
  V2Hat ~ lognormal(log(228), 1); // normal(0, 200);
  Q3Hat ~ lognormal(4.142, 0.0819);
  V3Hat ~ lognormal(6.590, 0.0864);
  Q4Hat ~ lognormal(5.679, 0.0967);
  V4Hat ~ lognormal(5.324, 0.0531);
  kaHat ~ lognormal(2.856, 0.0737);
  gammaCL ~ lognormal(log(0.75), 0.05);
  gammaV2 ~ lognormal(log(1), 0.05);

  omegaCL ~ lognormal(log(0.339), 0.5); // sd(log(CL))
  omegaV2 ~ lognormal(log(0.460), 0.5); // sd(log(V2))
  omegasq[1] ~ normal(0.347, 0.102); // var(log(Q3))
  omegasq[2] ~ normal(0.185, 0.0497); // var(log(V3))
  omegasq[3] ~ normal(0.146, 0.0525); // var(log(Q4))
  omegasq[4] ~ normal(0.0363, 0.0242); // var(log(V4))
  omegasq[5] ~ normal(0.246, 0.0531); // var(log(ka))

  //  rho ~ lkj_corr(2); 
  L ~ lkj_corr_cholesky(2);
  sigma ~ cauchy(0, 0.25);

  // Inter-individual variability
  //  logtheta ~ multi_normal(log(thetaHat), Omega);
  to_vector(eta) ~ normal(0, 1);

  cObs ~ normal(cHat[iObs], cHat[iObs] * sigma); // observed data likelihood
  // truncate below at 0
  target += -normal_lccdf(0.0 | cHat[iObs], cHat[iObs] * sigma);

}

generated quantities{
  matrix[nRandom, nID] etaPred;
  matrix<lower = 0>[nID, nRandom] thetaPred;
  corr_matrix[2] rho;
  vector<lower = 0>[nt] cHatPred;
  vector[nt] cObsCond;
  vector[nt] cObsPred;

  // Individual-level model parameters with recognizable names, e.g.,
  array[nID] real<lower = 0> CLPred;
  array[nID] real<lower = 0> V2Pred;
  array[nID] real<lower = 0> Q3Pred;
  array[nID] real<lower = 0> V3Pred;
  array[nID] real<lower = 0> Q4Pred;
  array[nID] real<lower = 0> V4Pred;
  array[nID] real<lower = 0> kaPred;

  matrix[nCmt, nt] xPred;
  matrix[nCmt, nCmt] KPred;

  rho = L * L';
  for(j in 1:nID) 
    for(i in 1:nRandom)
      etaPred[i, j] = normal_rng(0, 1);

  thetaPred[, 1:2] = (rep_matrix(thetaHat[1:2], nID) .* 
          exp(diag_pre_multiply(omega[1:2], L * etaPred[1:2,])))';
  for(i in 3:nRandom)
    thetaPred[, i] = thetaHat[i] * exp(omega[i] * etaPred[i,])';

  for(j in 1:nID){

    // Calculation of individual parameter values given theta and covariates, e.g.
    CLPred[j] = thetaPred[j, 1] * (weight[j] / 70)^gammaCL;
    V2Pred[j] = thetaPred[j, 2] * (weight[j] / 70)^gammaV2;
    Q3Pred[j] = thetaPred[j, 3] * (weight[j] / 70)^0.75;
    V3Pred[j] = thetaPred[j, 4] * (weight[j] / 70);
    Q4Pred[j] = thetaPred[j, 5] * (weight[j] / 70)^0.75;
    V4Pred[j] = thetaPred[j, 6] * (weight[j] / 70);
    kaPred[j] = thetaPred[j, 7];

    // Pack individual PK parameters into K matrix, e.g.
    
    KPred = rep_matrix(0, nCmt, nCmt);
    
    KPred[1, 1] = -kaPred[j];
    KPred[2, 1] = kaPred[j];
    KPred[2, 2] = -(CLPred[j] + Q3Pred[j] + Q4Pred[j]) / V2Pred[j];
    KPred[2, 3] = Q3Pred[j] / V3Pred[j];
    KPred[2, 4] = Q4Pred[j] / V4Pred[j];
    KPred[3, 2] = Q3Pred[j] / V2Pred[j];
    KPred[3, 3] = -Q3Pred[j] / V3Pred[j];
    KPred[4, 2] = Q4Pred[j] / V2Pred[j];
    KPred[4, 4] = -Q4Pred[j] / V4Pred[j];

    xPred[, start[j]:end[j]] = pmx_solve_linode(time[start[j]:end[j]], 
                    amt[start[j]:end[j]],
                    rate[start[j]:end[j]],
                    ii[start[j]:end[j]],
                    evid[start[j]:end[j]],
                    cmt[start[j]:end[j]],
                    addl[start[j]:end[j]],
                    ss[start[j]:end[j]],
                    KPred, F, tLag);

    // Calculate target concentration for specified compartment.
    // Change compartment number and distribution volume as appropriate.

    cHatPred[start[j]:end[j]] = xPred[2, start[j]:end[j]]' ./ V2Pred[j];
  }

  for(i in 1:nt){
    if(time[i] == 0){
      cObsCond[i] = 0;
      cObsPred[i] = 0;
    }else{
      cObsCond[i] = normal_trunc_rng(cHat[i], cHat[i] * sigma, 
         0.0, positive_infinity());
      cObsPred[i] = normal_trunc_rng(cHatPred[i], cHatPred[i] * sigma, 
         0.0, positive_infinity());
    }
  }
}
