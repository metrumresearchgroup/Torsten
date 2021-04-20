data{
  int<lower = 1> nId;
  int<lower = 1> nt;
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
  real<lower = 0> CLHat;
  real<lower = 0> VHat;
  real<lower = 0> kaHat;
  real<lower = 0> ke0Hat;
  real<lower = 0> EmaxHat;
  real<lower = 0> EC50;
  real<lower = 0> E0Hat;
  real<lower = 0> gamma;
  int<lower = 0> nRandom;
  corr_matrix[nRandom] rho;
  vector<lower = 0>[nRandom] omega;
  real<lower = 0> sigma;
  real<lower = 0> sigmaPD;
}

transformed data{
  int<lower = 1> nCmt = 3;
  real F[nCmt] = rep_array(1.0, nCmt);
  real tLag[nCmt] = rep_array(0.0, nCmt);
  vector<lower = 0>[nRandom] thetaHat =
    to_vector({CLHat, VHat, kaHat, ke0Hat, EmaxHat, E0Hat});
  cov_matrix[nRandom] Omega;

  Omega = quad_form_diag(rho, omega); // diag_matrix(omega) * rho * diag_matrix(omega)
}

parameters{
}

transformed parameters{
}

model{
}

generated quantities{
  vector[nRandom] logtheta[nId];
  real<lower = 0> CL[nId];
  real<lower = 0> V[nId];
  real<lower = 0> ka[nId];
  real<lower = 0> ke0[nId];
  real<lower = 0> Emax[nId];
  real<lower = 0> E0[nId];
  vector[nt] cHat;
  vector[nt] ceHat;
  vector[nt] respHat;
  real cObs[nt];
  real respObs[nt];
  matrix[nCmt, nCmt] K;
  matrix[nCmt, nt] x;

  for(j in 1:nId){
    logtheta[j] = multi_normal_rng(log(thetaHat), Omega);
    CL[j] = exp(logtheta[j, 1]) * (weight[j] / 70)^0.75;
    V[j] = exp(logtheta[j, 2]) * weight[j] / 70;
    ka[j] = exp(logtheta[j, 3]);
    ke0[j] = exp(logtheta[j, 4]);
    Emax[j] = exp(logtheta[j, 5]);
    E0[j] = exp(logtheta[j, 6]);

    K = rep_matrix(0, nCmt, nCmt);
    
    K[1, 1] = -ka[j];
    K[2, 1] = ka[j];
    K[2, 2] = -CL[j] / V[j];
    K[3, 2] = ke0[j];
    K[3, 3] = -ke0[j];

    x[, start[j]:end[j]] = pmx_solve_linode(time[start[j]:end[j]], 
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

  for(i in 1:nt){
    if(time[i] == 0){
      cObs[i] = 0;
    }else{
      cObs[i] = exp(normal_rng(log(fmax(machine_precision(), cHat[i])), sigma));
    }
    respObs[i] = exp(normal_rng(log(fmax(machine_precision(), respHat[i])), sigmaPD));
  }
}
