functions{

  real[] twoCptIndirectModelODE(real t,
				real[] x,
				real[] x_pk,
				real[] parms,
				real[] rdummy,
				int[] idummy){
    real V1      = parms[3];
    real kout    = parms[6];
    real effect0 = parms[7];
    real EC50    = parms[8];

/*   parms = {CL, Q, V1, V2, ka, kout, effect0, ec50}; */
    real conc = x_pk[2] / V1;
    real Edrug = conc / (EC50 + conc);
    real kin0 = kout * effect0;
    real kin = kin0 * (1 - Edrug);

    real effect = x[1] + effect0;

    real dxdt[1];
  
    dxdt[1] = kin - kout * effect;

    return dxdt;
  }
}

data{
  int<lower = 1> nt;
  int<lower = 1> nObs;  // number of observation
  int<lower = 1> iObs[nObs];  // index of observation
  real<lower = 0> amt[nt]; // mcg -- try 80000 mcg
  int<lower = 1> cmt[nt];
  int<lower = 0> evid[nt];
  real<lower = 0> time[nt];
  real<lower = 0> rate[nt];
  real<lower = 0> ii[nt];
  int<lower = 0> addl[nt];
  int<lower = 0> ss[nt];
}

transformed data{
  int nTheta = 8;
  real ka = 2.0; // 1/h
  real CL = 10; // L/h
  real V1 = 35; // L
  real V2 = 105; // L
  real Q = 15; // L/h
  real kout = 0.05; // 1/h
  real effect0 = 10; // units
  real EC50 = 400; // ng/mL
  int<lower = 1> nCmt = 1;

  real theta_data[nt, nTheta];
  real biovar_data[nt, nCmt];
  real tlag_data[nt, nCmt];

  real rtol;
  real atol;
  int max_step;

  for (i in 1:nt) {
    theta_data[i] = { CL, Q, V1, V2, ka, kout, effect0, EC50};
    biovar_data[i] = rep_array(1.0, nCmt);
    tlag_data[i] = rep_array(0.0, nCmt);
  }

  rtol = 1e-8;
  atol = 1e-8;
  max_step = 100000;
}

parameters{
}

transformed parameters{
}

model{
}

/* generated quantities{ */
/*   vector[nt] cHat; */
/*   vector[nt] effectHat; */
/*   matrix[nt, nCmt] x; */
/*   real<lower = 0> parms[8]; */

/*   parms = {CL, Q, V1, V2, ka, kout, effect0, ec50}; */

/*   x =  mixOde2CptModel_rk45(twoCptIndirectModelODE, 1, */
/* 			   time, amt, rate, ii, evid, cmt, addl, ss, */
/* 			   parms, F, tLag,1e-6, 1e-6, 1e8); */

/*   cHat = x[, 2] / V1; */
/*   effectHat = x[, 4] + effect0; */

/* } */

generated quantities {
  real theta[nt, nTheta];  // ODE parameters
  real biovar[nt, nCmt];
  real tlag[nt, nCmt];

  matrix[nt, nCmt + 3] x; 
  vector[nt] cHat;

  vector[nObs] cHatDat1;
  vector[nObs] cHatDat2;
  vector[nObs] cHatDat3;
  vector[nObs] cHatDat4;
  vector[nObs] cHatDat5;
  vector[nObs] cHatDat6;
  vector[nObs] cHatDat7;
  vector[nObs] cHatDat8;

  vector[nObs] cHatObs1;
  vector[nObs] cHatObs2;
  vector[nObs] cHatObs3;
  vector[nObs] cHatObs4;
  vector[nObs] cHatObs5;
  vector[nObs] cHatObs6;
  vector[nObs] cHatObs7;
  vector[nObs] cHatObs8;

  /* theiovar: param, tlag: param */
  vector[nObs] cHatDPP1;
  vector[nObs] cHatDPP2;
  vector[nObs] cHatDPP3;
  vector[nObs] cHatDPP4;
  vector[nObs] cHatDPP5;
  vector[nObs] cHatDPP6;
  vector[nObs] cHatDPP7;
  vector[nObs] cHatDPP8;

  /* thebiovar: data, tlag: param */
  vector[nObs] cHatPDP1;
  vector[nObs] cHatPDP2;
  vector[nObs] cHatPDP3;
  vector[nObs] cHatPDP4;
  vector[nObs] cHatPDP5;
  vector[nObs] cHatPDP6;
  vector[nObs] cHatPDP7;
  vector[nObs] cHatPDP8;

  /* thebiovar: param, tlag: data */
  vector[nObs] cHatPPD1;
  vector[nObs] cHatPPD2;
  vector[nObs] cHatPPD3;
  vector[nObs] cHatPPD4;
  vector[nObs] cHatPPD5;
  vector[nObs] cHatPPD6;
  vector[nObs] cHatPPD7;
  vector[nObs] cHatPPD8;

  /* theiovar: data, tlag: parm */
  vector[nObs] cHatDDP1;
  vector[nObs] cHatDDP2;
  vector[nObs] cHatDDP3;
  vector[nObs] cHatDDP4;
  vector[nObs] cHatDDP5;
  vector[nObs] cHatDDP6;
  vector[nObs] cHatDDP7;
  vector[nObs] cHatDDP8;

  /* theiovar: param, tlag: data */
  vector[nObs] cHatDPD1;
  vector[nObs] cHatDPD2;
  vector[nObs] cHatDPD3;
  vector[nObs] cHatDPD4;
  vector[nObs] cHatDPD5;
  vector[nObs] cHatDPD6;
  vector[nObs] cHatDPD7;
  vector[nObs] cHatDPD8;

  /* thebiovar: data, tlag: data */
  vector[nObs] cHatPDD1;
  vector[nObs] cHatPDD2;
  vector[nObs] cHatPDD3;
  vector[nObs] cHatPDD4;
  vector[nObs] cHatPDD5;
  vector[nObs] cHatPDD6;
  vector[nObs] cHatPDD7;
  vector[nObs] cHatPDD8;

  for (i in 1:nt) {
    theta[i] = { CL, Q, V1, V2, ka, kout, effect0, EC50};
    biovar[i] = rep_array(1.0, nCmt);
    tlag[i] = rep_array(0.0, nCmt);
  }
  
  /* /\* data args *\/ */
  /* x = mixOde2CptModel_rk45(twoCptIndirectModelODE, nCmt, time, amt, rate, ii, evid, cmt, addl, ss, theta_data    , biovar_data    , tlag_data,    rtol, atol, max_step);  cHat = col(x, 4);  cHatDat1 = cHat[iObs]; */

  /* data args */
  x = mixOde2CptModel_rk45(twoCptIndirectModelODE, nCmt, time, amt, rate, ii, evid, cmt, addl, ss, theta_data    , biovar_data    , tlag_data,    rtol, atol, max_step);  cHat = col(x, 4);  cHatDat1 = cHat[iObs];
  x = mixOde2CptModel_rk45(twoCptIndirectModelODE, nCmt, time, amt, rate, ii, evid, cmt, addl, ss, theta_data[1] , biovar_data    , tlag_data,    rtol, atol, max_step);  cHat = col(x, 4);  cHatDat2 = cHat[iObs];
  x = mixOde2CptModel_rk45(twoCptIndirectModelODE, nCmt, time, amt, rate, ii, evid, cmt, addl, ss, theta_data[1] , biovar_data[1] , tlag_data,    rtol, atol, max_step);  cHat = col(x, 4);  cHatDat3 = cHat[iObs];
  x = mixOde2CptModel_rk45(twoCptIndirectModelODE, nCmt, time, amt, rate, ii, evid, cmt, addl, ss, theta_data[1] , biovar_data[1] , tlag_data[1], rtol, atol, max_step);  cHat = col(x, 4);  cHatDat4 = cHat[iObs];
  x = mixOde2CptModel_rk45(twoCptIndirectModelODE, nCmt, time, amt, rate, ii, evid, cmt, addl, ss, theta_data[1] , biovar_data    , tlag_data[1], rtol, atol, max_step);  cHat = col(x, 4);  cHatDat5 = cHat[iObs];
  x = mixOde2CptModel_rk45(twoCptIndirectModelODE, nCmt, time, amt, rate, ii, evid, cmt, addl, ss, theta_data    , biovar_data[1] , tlag_data,    rtol, atol, max_step);  cHat = col(x, 4);  cHatDat6 = cHat[iObs];
  x = mixOde2CptModel_rk45(twoCptIndirectModelODE, nCmt, time, amt, rate, ii, evid, cmt, addl, ss, theta_data    , biovar_data[1] , tlag_data[1], rtol, atol, max_step);  cHat = col(x, 4);  cHatDat7 = cHat[iObs];
  x = mixOde2CptModel_rk45(twoCptIndirectModelODE, nCmt, time, amt, rate, ii, evid, cmt, addl, ss, theta_data    , biovar_data    , tlag_data,    rtol, atol, max_step);  cHat = col(x, 4);  cHatDat8 = cHat[iObs];

  /* param args */
  x = mixOde2CptModel_rk45(twoCptIndirectModelODE, nCmt, time, amt, rate, ii, evid, cmt, addl, ss, theta         , biovar         , tlag,         rtol, atol, max_step);  cHat = col(x, 4);  cHatObs1 = cHat[iObs];
  x = mixOde2CptModel_rk45(twoCptIndirectModelODE, nCmt, time, amt, rate, ii, evid, cmt, addl, ss, theta[1]      , biovar         , tlag,         rtol, atol, max_step);  cHat = col(x, 4);  cHatObs2 = cHat[iObs];  
  x = mixOde2CptModel_rk45(twoCptIndirectModelODE, nCmt, time, amt, rate, ii, evid, cmt, addl, ss, theta[1]      , biovar[1]      , tlag,         rtol, atol, max_step);  cHat = col(x, 4);  cHatObs3 = cHat[iObs];
  x = mixOde2CptModel_rk45(twoCptIndirectModelODE, nCmt, time, amt, rate, ii, evid, cmt, addl, ss, theta[1]      , biovar[1]      , tlag[1],      rtol, atol, max_step);  cHat = col(x, 4);  cHatObs4 = cHat[iObs];
  x = mixOde2CptModel_rk45(twoCptIndirectModelODE, nCmt, time, amt, rate, ii, evid, cmt, addl, ss, theta[1]      , biovar         , tlag[1],      rtol, atol, max_step);  cHat = col(x, 4);  cHatObs5 = cHat[iObs];
  x = mixOde2CptModel_rk45(twoCptIndirectModelODE, nCmt, time, amt, rate, ii, evid, cmt, addl, ss, theta         , biovar[1]      , tlag,         rtol, atol, max_step);  cHat = col(x, 4);  cHatObs6 = cHat[iObs];
  x = mixOde2CptModel_rk45(twoCptIndirectModelODE, nCmt, time, amt, rate, ii, evid, cmt, addl, ss, theta         , biovar[1]      , tlag[1],      rtol, atol, max_step);  cHat = col(x, 4);  cHatObs7 = cHat[iObs];
  x = mixOde2CptModel_rk45(twoCptIndirectModelODE, nCmt, time, amt, rate, ii, evid, cmt, addl, ss, theta         , biovar         , tlag[1],      rtol, atol, max_step);  cHat = col(x, 4);  cHatObs8 = cHat[iObs];

  /* theta: data                                                            , biovar: param  , tlag: param */
  x = mixOde2CptModel_rk45(twoCptIndirectModelODE, nCmt, time, amt, rate, ii, evid, cmt, addl, ss, theta_data    , biovar         , tlag,         rtol, atol, max_step);  cHat = col(x, 4);  cHatDPP1 = cHat[iObs];
  x = mixOde2CptModel_rk45(twoCptIndirectModelODE, nCmt, time, amt, rate, ii, evid, cmt, addl, ss, theta_data[1] , biovar         , tlag,         rtol, atol, max_step);  cHat = col(x, 4);  cHatDPP2 = cHat[iObs];  
  x = mixOde2CptModel_rk45(twoCptIndirectModelODE, nCmt, time, amt, rate, ii, evid, cmt, addl, ss, theta_data[1] , biovar[1]      , tlag,         rtol, atol, max_step);  cHat = col(x, 4);  cHatDPP3 = cHat[iObs];
  x = mixOde2CptModel_rk45(twoCptIndirectModelODE, nCmt, time, amt, rate, ii, evid, cmt, addl, ss, theta_data[1] , biovar[1]      , tlag[1],      rtol, atol, max_step);  cHat = col(x, 4);  cHatDPP4 = cHat[iObs];
  x = mixOde2CptModel_rk45(twoCptIndirectModelODE, nCmt, time, amt, rate, ii, evid, cmt, addl, ss, theta_data[1] , biovar         , tlag[1],      rtol, atol, max_step);  cHat = col(x, 4);  cHatDPP5 = cHat[iObs];
  x = mixOde2CptModel_rk45(twoCptIndirectModelODE, nCmt, time, amt, rate, ii, evid, cmt, addl, ss, theta_data    , biovar[1]      , tlag,         rtol, atol, max_step);  cHat = col(x, 4);  cHatDPP6 = cHat[iObs];
  x = mixOde2CptModel_rk45(twoCptIndirectModelODE, nCmt, time, amt, rate, ii, evid, cmt, addl, ss, theta_data    , biovar[1]      , tlag[1],      rtol, atol, max_step);  cHat = col(x, 4);  cHatDPP7 = cHat[iObs];
  x = mixOde2CptModel_rk45(twoCptIndirectModelODE, nCmt, time, amt, rate, ii, evid, cmt, addl, ss, theta_data    , biovar         , tlag[1],      rtol, atol, max_step);  cHat = col(x, 4);  cHatDPP8 = cHat[iObs];

  /* theta: parm                                                            , biovar: data   , tlag: param */
  x = mixOde2CptModel_rk45(twoCptIndirectModelODE, nCmt, time, amt, rate, ii, evid, cmt, addl, ss, theta         , biovar_data    , tlag,         rtol, atol, max_step);  cHat = col(x, 4);  cHatPDP1 = cHat[iObs];
  x = mixOde2CptModel_rk45(twoCptIndirectModelODE, nCmt, time, amt, rate, ii, evid, cmt, addl, ss, theta[1]      , biovar_data    , tlag,         rtol, atol, max_step);  cHat = col(x, 4);  cHatPDP2 = cHat[iObs];  
  x = mixOde2CptModel_rk45(twoCptIndirectModelODE, nCmt, time, amt, rate, ii, evid, cmt, addl, ss, theta[1]      , biovar_data[1] , tlag,         rtol, atol, max_step);  cHat = col(x, 4);  cHatPDP3 = cHat[iObs];
  x = mixOde2CptModel_rk45(twoCptIndirectModelODE, nCmt, time, amt, rate, ii, evid, cmt, addl, ss, theta[1]      , biovar_data[1] , tlag[1],      rtol, atol, max_step);  cHat = col(x, 4);  cHatPDP4 = cHat[iObs];
  x = mixOde2CptModel_rk45(twoCptIndirectModelODE, nCmt, time, amt, rate, ii, evid, cmt, addl, ss, theta[1]      , biovar_data    , tlag[1],      rtol, atol, max_step);  cHat = col(x, 4);  cHatPDP5 = cHat[iObs];
  x = mixOde2CptModel_rk45(twoCptIndirectModelODE, nCmt, time, amt, rate, ii, evid, cmt, addl, ss, theta         , biovar_data[1] , tlag,         rtol, atol, max_step);  cHat = col(x, 4);  cHatPDP6 = cHat[iObs];
  x = mixOde2CptModel_rk45(twoCptIndirectModelODE, nCmt, time, amt, rate, ii, evid, cmt, addl, ss, theta         , biovar_data[1] , tlag[1],      rtol, atol, max_step);  cHat = col(x, 4);  cHatPDP7 = cHat[iObs];
  x = mixOde2CptModel_rk45(twoCptIndirectModelODE, nCmt, time, amt, rate, ii, evid, cmt, addl, ss, theta         , biovar_data    , tlag[1],      rtol, atol, max_step);  cHat = col(x, 4);  cHatPDP8 = cHat[iObs];

  /* theta: parm                                                            , biovar: parm   , tlag: data */
  x = mixOde2CptModel_rk45(twoCptIndirectModelODE, nCmt, time, amt, rate, ii, evid, cmt, addl, ss, theta         , biovar         , tlag_data,    rtol, atol, max_step);  cHat = col(x, 4);  cHatPPD1 = cHat[iObs];
  x = mixOde2CptModel_rk45(twoCptIndirectModelODE, nCmt, time, amt, rate, ii, evid, cmt, addl, ss, theta[1]      , biovar         , tlag_data,    rtol, atol, max_step);  cHat = col(x, 4);  cHatPPD2 = cHat[iObs];  
  x = mixOde2CptModel_rk45(twoCptIndirectModelODE, nCmt, time, amt, rate, ii, evid, cmt, addl, ss, theta[1]      , biovar[1]      , tlag_data,    rtol, atol, max_step);  cHat = col(x, 4);  cHatPPD3 = cHat[iObs];
  x = mixOde2CptModel_rk45(twoCptIndirectModelODE, nCmt, time, amt, rate, ii, evid, cmt, addl, ss, theta[1]      , biovar[1]      , tlag_data[1], rtol, atol, max_step);  cHat = col(x, 4);  cHatPPD4 = cHat[iObs];
  x = mixOde2CptModel_rk45(twoCptIndirectModelODE, nCmt, time, amt, rate, ii, evid, cmt, addl, ss, theta[1]      , biovar         , tlag_data[1], rtol, atol, max_step);  cHat = col(x, 4);  cHatPPD5 = cHat[iObs];
  x = mixOde2CptModel_rk45(twoCptIndirectModelODE, nCmt, time, amt, rate, ii, evid, cmt, addl, ss, theta         , biovar[1]      , tlag_data,    rtol, atol, max_step);  cHat = col(x, 4);  cHatPPD6 = cHat[iObs];
  x = mixOde2CptModel_rk45(twoCptIndirectModelODE, nCmt, time, amt, rate, ii, evid, cmt, addl, ss, theta         , biovar[1]      , tlag_data[1], rtol, atol, max_step);  cHat = col(x, 4);  cHatPPD7 = cHat[iObs];
  x = mixOde2CptModel_rk45(twoCptIndirectModelODE, nCmt, time, amt, rate, ii, evid, cmt, addl, ss, theta         , biovar         , tlag_data[1], rtol, atol, max_step);  cHat = col(x, 4);  cHatPPD8 = cHat[iObs];

  /* theta: data                                                            , biovar: data   , tlag: param */
  x = mixOde2CptModel_rk45(twoCptIndirectModelODE, nCmt, time, amt, rate, ii, evid, cmt, addl, ss, theta_data    , biovar_data    , tlag,         rtol, atol, max_step);  cHat = col(x, 4);  cHatDDP1 = cHat[iObs];
  x = mixOde2CptModel_rk45(twoCptIndirectModelODE, nCmt, time, amt, rate, ii, evid, cmt, addl, ss, theta_data[1] , biovar_data    , tlag,         rtol, atol, max_step);  cHat = col(x, 4);  cHatDDP2 = cHat[iObs];  
  x = mixOde2CptModel_rk45(twoCptIndirectModelODE, nCmt, time, amt, rate, ii, evid, cmt, addl, ss, theta_data[1] , biovar_data[1] , tlag,         rtol, atol, max_step);  cHat = col(x, 4);  cHatDDP3 = cHat[iObs];
  x = mixOde2CptModel_rk45(twoCptIndirectModelODE, nCmt, time, amt, rate, ii, evid, cmt, addl, ss, theta_data[1] , biovar_data[1] , tlag[1],      rtol, atol, max_step);  cHat = col(x, 4);  cHatDDP4 = cHat[iObs];
  x = mixOde2CptModel_rk45(twoCptIndirectModelODE, nCmt, time, amt, rate, ii, evid, cmt, addl, ss, theta_data[1] , biovar_data    , tlag[1],      rtol, atol, max_step);  cHat = col(x, 4);  cHatDDP5 = cHat[iObs];
  x = mixOde2CptModel_rk45(twoCptIndirectModelODE, nCmt, time, amt, rate, ii, evid, cmt, addl, ss, theta_data    , biovar_data[1] , tlag,         rtol, atol, max_step);  cHat = col(x, 4);  cHatDDP6 = cHat[iObs];
  x = mixOde2CptModel_rk45(twoCptIndirectModelODE, nCmt, time, amt, rate, ii, evid, cmt, addl, ss, theta_data    , biovar_data[1] , tlag[1],      rtol, atol, max_step);  cHat = col(x, 4);  cHatDDP7 = cHat[iObs];
  x = mixOde2CptModel_rk45(twoCptIndirectModelODE, nCmt, time, amt, rate, ii, evid, cmt, addl, ss, theta_data    , biovar_data    , tlag[1],      rtol, atol, max_step);  cHat = col(x, 4);  cHatDDP8 = cHat[iObs];

  /* theta: data                                                            , biovar: param  , tlag: data */
  x = mixOde2CptModel_rk45(twoCptIndirectModelODE, nCmt, time, amt, rate, ii, evid, cmt, addl, ss, theta_data    , biovar         , tlag_data,    rtol, atol, max_step);  cHat = col(x, 4);  cHatDPD1 = cHat[iObs];
  x = mixOde2CptModel_rk45(twoCptIndirectModelODE, nCmt, time, amt, rate, ii, evid, cmt, addl, ss, theta_data[1] , biovar         , tlag_data,    rtol, atol, max_step);  cHat = col(x, 4);  cHatDPD2 = cHat[iObs];  
  x = mixOde2CptModel_rk45(twoCptIndirectModelODE, nCmt, time, amt, rate, ii, evid, cmt, addl, ss, theta_data[1] , biovar[1]      , tlag_data,    rtol, atol, max_step);  cHat = col(x, 4);  cHatDPD3 = cHat[iObs];
  x = mixOde2CptModel_rk45(twoCptIndirectModelODE, nCmt, time, amt, rate, ii, evid, cmt, addl, ss, theta_data[1] , biovar[1]      , tlag_data[1], rtol, atol, max_step);  cHat = col(x, 4);  cHatDPD4 = cHat[iObs];
  x = mixOde2CptModel_rk45(twoCptIndirectModelODE, nCmt, time, amt, rate, ii, evid, cmt, addl, ss, theta_data[1] , biovar         , tlag_data[1], rtol, atol, max_step);  cHat = col(x, 4);  cHatDPD5 = cHat[iObs];
  x = mixOde2CptModel_rk45(twoCptIndirectModelODE, nCmt, time, amt, rate, ii, evid, cmt, addl, ss, theta_data    , biovar[1]      , tlag_data,    rtol, atol, max_step);  cHat = col(x, 4);  cHatDPD6 = cHat[iObs];
  x = mixOde2CptModel_rk45(twoCptIndirectModelODE, nCmt, time, amt, rate, ii, evid, cmt, addl, ss, theta_data    , biovar[1]      , tlag_data[1], rtol, atol, max_step);  cHat = col(x, 4);  cHatDPD7 = cHat[iObs];
  x = mixOde2CptModel_rk45(twoCptIndirectModelODE, nCmt, time, amt, rate, ii, evid, cmt, addl, ss, theta_data    , biovar         , tlag_data[1], rtol, atol, max_step);  cHat = col(x, 4);  cHatDPD8 = cHat[iObs];

  /* theta: param                                                           , biovar: data   , tlag: data */
  x = mixOde2CptModel_rk45(twoCptIndirectModelODE, nCmt, time, amt, rate, ii, evid, cmt, addl, ss, theta         , biovar_data    , tlag_data,    rtol, atol, max_step);  cHat = col(x, 4);  cHatPDD1 = cHat[iObs];
  x = mixOde2CptModel_rk45(twoCptIndirectModelODE, nCmt, time, amt, rate, ii, evid, cmt, addl, ss, theta[1]      , biovar_data    , tlag_data,    rtol, atol, max_step);  cHat = col(x, 4);  cHatPDD2 = cHat[iObs];  
  x = mixOde2CptModel_rk45(twoCptIndirectModelODE, nCmt, time, amt, rate, ii, evid, cmt, addl, ss, theta[1]      , biovar_data[1] , tlag_data,    rtol, atol, max_step);  cHat = col(x, 4);  cHatPDD3 = cHat[iObs];
  x = mixOde2CptModel_rk45(twoCptIndirectModelODE, nCmt, time, amt, rate, ii, evid, cmt, addl, ss, theta[1]      , biovar_data[1] , tlag_data[1], rtol, atol, max_step);  cHat = col(x, 4);  cHatPDD4 = cHat[iObs];
  x = mixOde2CptModel_rk45(twoCptIndirectModelODE, nCmt, time, amt, rate, ii, evid, cmt, addl, ss, theta[1]      , biovar_data    , tlag_data[1], rtol, atol, max_step);  cHat = col(x, 4);  cHatPDD5 = cHat[iObs];
  x = mixOde2CptModel_rk45(twoCptIndirectModelODE, nCmt, time, amt, rate, ii, evid, cmt, addl, ss, theta         , biovar_data[1] , tlag_data,    rtol, atol, max_step);  cHat = col(x, 4);  cHatPDD6 = cHat[iObs];
  x = mixOde2CptModel_rk45(twoCptIndirectModelODE, nCmt, time, amt, rate, ii, evid, cmt, addl, ss, theta         , biovar_data[1] , tlag_data[1], rtol, atol, max_step);  cHat = col(x, 4);  cHatPDD7 = cHat[iObs];
  x = mixOde2CptModel_rk45(twoCptIndirectModelODE, nCmt, time, amt, rate, ii, evid, cmt, addl, ss, theta         , biovar_data    , tlag_data[1], rtol, atol, max_step);  cHat = col(x, 4);  cHatPDD8 = cHat[iObs];
}

