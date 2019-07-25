functions {
  
  // define ODE system for two compartmnt model
  real[] twoCptModelODE(real t,
			real[] x,
			real[] parms,
			real[] rate,  // in this example, rate is treated as data
			int[] dummy){
			  
    // Parameters
    real CL = parms[1];
    real Q = parms[2];
    real V1 = parms[3];
    real V2 = parms[4];
    real ka = parms[5];
    
    // Re-parametrization
    real k10 = CL / V1;
    real k12 = Q / V1;
    real k21 = Q / V2;
    
    // Return object (derivative)
    real y[3];  // 1 element per compartment of
                // the model

    // PK component of the ODE system
    y[1] = -ka*x[1];
    y[2] = ka*x[1] - (k10 + k12)*x[2] + k21*x[3];
    y[3] = k12*x[2] - k21*x[3];

    return y;
  }
  
}

data{
  int<lower = 1> nt;  // number of events
  int<lower = 1> nObs;  // number of observation
  int<lower = 1> iObs[nObs];  // index of observation
  
  // NONMEM data
  int<lower = 1> cmt[nt];
  int evid[nt];
  int addl[nt];
  int ss[nt];
  real amt[nt];
  real time[nt];
  real rate[nt];
  real ii[nt];
}

transformed data{
  int nTheta = 5;  // number of ODE parameters in Two Compartment Model
  int nCmt = 3;  // number of compartments in model

  // Since we're not trying to evaluate the bio-variability (F) and 
  // the lag times, we declare them as data.
  real theta_data[nt, nTheta];
  real biovar_data[nt, nCmt];
  real tlag_data[nt, nCmt];

  real<lower = 0> CL_data;
  real<lower = 0> Q_data;
  real<lower = 0> V2_data;
  real<lower = 0> V3_data;
  real<lower = 0> ka_data;

  real rtol;
  real atol;
  int max_step;

  rtol = 1e-8;
  atol = 1e-8;
  max_step = 100000;

  CL_data    = 5;
  ka_data    = 1.2;
  Q_data     = 8;
  V2_data    = 20;
  V3_data    = 70;

  for (i in 1:nt) {
    theta_data[i] = { CL_data, Q_data, V2_data, V3_data, ka_data };    
    biovar_data[i] = { 1, 1, 1 };
    tlag_data[i] = { 0, 0, 0 };
  }

}

generated quantities {
  real theta[nt, nTheta];  // ODE parameters
  real biovar[nt, nCmt];
  real tlag[nt, nCmt];

  matrix<lower = 0>[nt, nCmt] x; 
  vector<lower = 0>[nt] cHat;

  real<lower = 0> CL;
  real<lower = 0> Q;
  real<lower = 0> V2;
  real<lower = 0> V3;
  real<lower = 0> ka;
  real<lower = 0> sigma;

  vector<lower = 0>[nObs] cHatDat1;

  CL    = 5;
  ka    = 1.2;
  Q     = 8;
  sigma = 0.01;
  V2    = 20;                   /* central cpt */
  V3    = 70;

  for (i in 1:nt) {
    theta[i] = { CL, Q, V2, V3, ka };
    biovar[i] = { 1, 1, 1 };
    tlag[i] = { 0, 0, 0 };
  }

  /* data args */
  x = generalOdeModel_rk45(twoCptModelODE, nCmt, time, amt, rate, ii, evid, cmt, addl, ss, theta_data    , biovar_data    , tlag_data);  cHat = col(x, 2) ./ V2;  cHatDat1 = cHat[iObs];
}

