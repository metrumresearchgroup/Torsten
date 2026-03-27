functions{

    real[] oneCptPNODE(real t,
			real[] x,
			real[] parms,
			real[] rdummy,
			int[] idummy){
    real dxdt[3];
    real CL = parms[1];
    real V = parms[2];
    real ke0 = parms[3];
    real alpha = parms[4];
    real beta = parms[5];
    real Edrug = alpha * x[2];
    real hazard;

    dxdt[1] = -(CL / V) * x[1];
    dxdt[2] = ke0 * (x[1] / V - x[2]);
    hazard = beta * Edrug^beta * t^(beta - 1);
    dxdt[3] = hazard;
    
    return dxdt;
  }

}

data{
  int<lower = 1> nId;
  int<lower = 1> nt;
  int<lower = 1> nPKObs;
  int<lower = 1> iPKObs[nPKObs];
  real<lower = 0> amt[nt];
  real<lower = 0> rate[nt];
  real<lower = 0> ii[nt];
  int<lower = 0> addl[nt];
  int<lower = 1> cmt[nt];
  int<lower = 0> evid[nt];
  int<lower = 1> start[nId];
  int<lower = 1> end[nId];
  real<lower = 0> time[nt];
  real<lower = 0> CLHat;
  real<lower = 0> VHat;
  real<lower = 0> ke0;
  real<lower = 0> alpha;
  real<lower = 0> beta;
  int<lower = 1> nRandom;
  vector<lower = 0>[nRandom] omega;
  corr_matrix[nRandom] rho;
  real<lower = 0> sigma;
}

transformed data{
  int<lower = 1> nCmt = 3;
  int<lower = 0> ss[nt] = rep_array(0, nt);
  real F[nCmt] = rep_array(1.0, nCmt);
  real tLag[nCmt] = rep_array(0.0, nCmt);
  cov_matrix[nRandom] Omega;
  
  Omega = quad_form_diag(rho, omega);
}

parameters{
}

transformed parameters{
}

model{
}

generated quantities{
  vector<lower = 0>[nRandom] thetaHat = to_vector({CLHat, VHat});
  real<lower = 0> CL[nId];
  real<lower = 0> V[nId];
  vector<lower = 0>[nRandom] theta[nId];

  vector[nt] cHat;
  vector[nt] ce;
  vector[nPKObs] cObs;
  vector<lower = 0, upper = 1>[nt] cdf;
  matrix[nt, 3] x;
  real<lower = 0> parms[5];

  for(j in 1:nId){
    theta[j,] = exp(multi_normal_rng(log(thetaHat), Omega));
    CL[j] = theta[j, 1];
    V[j] = theta[j, 2];

    parms = {CL[j], V[j], ke0, alpha, beta};

    x[start[j]:end[j],] = generalOdeModel_rk45(oneCptPNODE, 3,
			     time[start[j]:end[j]], 
			     amt[start[j]:end[j]],
			     rate[start[j]:end[j]],
			     ii[start[j]:end[j]],
			     evid[start[j]:end[j]],
			     cmt[start[j]:end[j]],
			     addl[start[j]:end[j]],
			     ss[start[j]:end[j]],
			     parms, F, tLag,
			     1e-6, 1e-6, 1e8);

    cHat[start[j]:end[j]] = x[start[j]:end[j], 1] / V[j];
    ce[start[j]:end[j]] = x[start[j]:end[j], 2];
  }

  cdf = 1 - exp(-x[, 3]);
  for(i in 1:nPKObs)
    cObs[i] = exp(normal_rng(log(cHat[iPKObs[i]]), sigma));
}
