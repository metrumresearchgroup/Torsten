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

generated quantities {
  real theta[nt, nTheta];  // ODE parameters
  real biovar[nt, nCmt];
  real tlag[nt, nCmt];

  matrix[nt, nCmt + 3] x; 
  vector[nt] cHat;

  for (i in 1:nt) {
    theta[i] = { CL, Q, V1, V2, ka, kout, effect0, EC50};
    biovar[i] = rep_array(1.0, nCmt);
    tlag[i] = rep_array(0.0, nCmt);
  }
  
  /* /\* data args *\/ */
  /* x = mixOde2CptModel_rk45(twoCptIndirectModelODE, nCmt, time, amt, rate, ii, evid, cmt, addl, ss, theta_data    , biovar_data    , tlag_data,    rtol, atol, max_step);  cHat = col(x, 4);  cHatDat1 = cHat[iObs]; */

  /* data args */
  x = mixOde2CptModel_rk45(twoCptIndirectModelODE, nCmt, time, amt, rate, ii, evid, cmt, addl, ss, theta_data    , biovar_data    , tlag_data,    rtol, atol, max_step);
  x = mixOde2CptModel_rk45(twoCptIndirectModelODE, nCmt, time, amt, rate, ii, evid, cmt, addl, ss, theta_data[1] , biovar_data    , tlag_data,    rtol, atol, max_step);
  x = mixOde2CptModel_rk45(twoCptIndirectModelODE, nCmt, time, amt, rate, ii, evid, cmt, addl, ss, theta_data[1] , biovar_data[1] , tlag_data,    rtol, atol, max_step);
  x = mixOde2CptModel_rk45(twoCptIndirectModelODE, nCmt, time, amt, rate, ii, evid, cmt, addl, ss, theta_data[1] , biovar_data[1] , tlag_data[1], rtol, atol, max_step);
  x = mixOde2CptModel_rk45(twoCptIndirectModelODE, nCmt, time, amt, rate, ii, evid, cmt, addl, ss, theta_data[1] , biovar_data    , tlag_data[1], rtol, atol, max_step);
  x = mixOde2CptModel_rk45(twoCptIndirectModelODE, nCmt, time, amt, rate, ii, evid, cmt, addl, ss, theta_data    , biovar_data[1] , tlag_data,    rtol, atol, max_step);
  x = mixOde2CptModel_rk45(twoCptIndirectModelODE, nCmt, time, amt, rate, ii, evid, cmt, addl, ss, theta_data    , biovar_data[1] , tlag_data[1], rtol, atol, max_step);
  x = mixOde2CptModel_rk45(twoCptIndirectModelODE, nCmt, time, amt, rate, ii, evid, cmt, addl, ss, theta_data    , biovar_data    , tlag_data,    rtol, atol, max_step);

  /* param args */
  x = mixOde2CptModel_rk45(twoCptIndirectModelODE, nCmt, time, amt, rate, ii, evid, cmt, addl, ss, theta         , biovar         , tlag,         rtol, atol, max_step);
  x = mixOde2CptModel_rk45(twoCptIndirectModelODE, nCmt, time, amt, rate, ii, evid, cmt, addl, ss, theta[1]      , biovar         , tlag,         rtol, atol, max_step);  
  x = mixOde2CptModel_rk45(twoCptIndirectModelODE, nCmt, time, amt, rate, ii, evid, cmt, addl, ss, theta[1]      , biovar[1]      , tlag,         rtol, atol, max_step);
  x = mixOde2CptModel_rk45(twoCptIndirectModelODE, nCmt, time, amt, rate, ii, evid, cmt, addl, ss, theta[1]      , biovar[1]      , tlag[1],      rtol, atol, max_step);
  x = mixOde2CptModel_rk45(twoCptIndirectModelODE, nCmt, time, amt, rate, ii, evid, cmt, addl, ss, theta[1]      , biovar         , tlag[1],      rtol, atol, max_step);
  x = mixOde2CptModel_rk45(twoCptIndirectModelODE, nCmt, time, amt, rate, ii, evid, cmt, addl, ss, theta         , biovar[1]      , tlag,         rtol, atol, max_step);
  x = mixOde2CptModel_rk45(twoCptIndirectModelODE, nCmt, time, amt, rate, ii, evid, cmt, addl, ss, theta         , biovar[1]      , tlag[1],      rtol, atol, max_step);
  x = mixOde2CptModel_rk45(twoCptIndirectModelODE, nCmt, time, amt, rate, ii, evid, cmt, addl, ss, theta         , biovar         , tlag[1],      rtol, atol, max_step);

  /* theta: data                                                            , biovar: param  , tlag: param */
  x = mixOde2CptModel_rk45(twoCptIndirectModelODE, nCmt, time, amt, rate, ii, evid, cmt, addl, ss, theta_data    , biovar         , tlag,         rtol, atol, max_step);
  x = mixOde2CptModel_rk45(twoCptIndirectModelODE, nCmt, time, amt, rate, ii, evid, cmt, addl, ss, theta_data[1] , biovar         , tlag,         rtol, atol, max_step);  
  x = mixOde2CptModel_rk45(twoCptIndirectModelODE, nCmt, time, amt, rate, ii, evid, cmt, addl, ss, theta_data[1] , biovar[1]      , tlag,         rtol, atol, max_step);
  x = mixOde2CptModel_rk45(twoCptIndirectModelODE, nCmt, time, amt, rate, ii, evid, cmt, addl, ss, theta_data[1] , biovar[1]      , tlag[1],      rtol, atol, max_step);
  x = mixOde2CptModel_rk45(twoCptIndirectModelODE, nCmt, time, amt, rate, ii, evid, cmt, addl, ss, theta_data[1] , biovar         , tlag[1],      rtol, atol, max_step);
  x = mixOde2CptModel_rk45(twoCptIndirectModelODE, nCmt, time, amt, rate, ii, evid, cmt, addl, ss, theta_data    , biovar[1]      , tlag,         rtol, atol, max_step);
  x = mixOde2CptModel_rk45(twoCptIndirectModelODE, nCmt, time, amt, rate, ii, evid, cmt, addl, ss, theta_data    , biovar[1]      , tlag[1],      rtol, atol, max_step);
  x = mixOde2CptModel_rk45(twoCptIndirectModelODE, nCmt, time, amt, rate, ii, evid, cmt, addl, ss, theta_data    , biovar         , tlag[1],      rtol, atol, max_step);

  /* theta: parm                                                            , biovar: data   , tlag: param */
  x = mixOde2CptModel_rk45(twoCptIndirectModelODE, nCmt, time, amt, rate, ii, evid, cmt, addl, ss, theta         , biovar_data    , tlag,         rtol, atol, max_step);
  x = mixOde2CptModel_rk45(twoCptIndirectModelODE, nCmt, time, amt, rate, ii, evid, cmt, addl, ss, theta[1]      , biovar_data    , tlag,         rtol, atol, max_step);  
  x = mixOde2CptModel_rk45(twoCptIndirectModelODE, nCmt, time, amt, rate, ii, evid, cmt, addl, ss, theta[1]      , biovar_data[1] , tlag,         rtol, atol, max_step);
  x = mixOde2CptModel_rk45(twoCptIndirectModelODE, nCmt, time, amt, rate, ii, evid, cmt, addl, ss, theta[1]      , biovar_data[1] , tlag[1],      rtol, atol, max_step);
  x = mixOde2CptModel_rk45(twoCptIndirectModelODE, nCmt, time, amt, rate, ii, evid, cmt, addl, ss, theta[1]      , biovar_data    , tlag[1],      rtol, atol, max_step);
  x = mixOde2CptModel_rk45(twoCptIndirectModelODE, nCmt, time, amt, rate, ii, evid, cmt, addl, ss, theta         , biovar_data[1] , tlag,         rtol, atol, max_step);
  x = mixOde2CptModel_rk45(twoCptIndirectModelODE, nCmt, time, amt, rate, ii, evid, cmt, addl, ss, theta         , biovar_data[1] , tlag[1],      rtol, atol, max_step);
  x = mixOde2CptModel_rk45(twoCptIndirectModelODE, nCmt, time, amt, rate, ii, evid, cmt, addl, ss, theta         , biovar_data    , tlag[1],      rtol, atol, max_step);

  /* theta: parm                                                            , biovar: parm   , tlag: data */
  x = mixOde2CptModel_rk45(twoCptIndirectModelODE, nCmt, time, amt, rate, ii, evid, cmt, addl, ss, theta         , biovar         , tlag_data,    rtol, atol, max_step);
  x = mixOde2CptModel_rk45(twoCptIndirectModelODE, nCmt, time, amt, rate, ii, evid, cmt, addl, ss, theta[1]      , biovar         , tlag_data,    rtol, atol, max_step);  
  x = mixOde2CptModel_rk45(twoCptIndirectModelODE, nCmt, time, amt, rate, ii, evid, cmt, addl, ss, theta[1]      , biovar[1]      , tlag_data,    rtol, atol, max_step);
  x = mixOde2CptModel_rk45(twoCptIndirectModelODE, nCmt, time, amt, rate, ii, evid, cmt, addl, ss, theta[1]      , biovar[1]      , tlag_data[1], rtol, atol, max_step);
  x = mixOde2CptModel_rk45(twoCptIndirectModelODE, nCmt, time, amt, rate, ii, evid, cmt, addl, ss, theta[1]      , biovar         , tlag_data[1], rtol, atol, max_step);
  x = mixOde2CptModel_rk45(twoCptIndirectModelODE, nCmt, time, amt, rate, ii, evid, cmt, addl, ss, theta         , biovar[1]      , tlag_data,    rtol, atol, max_step);
  x = mixOde2CptModel_rk45(twoCptIndirectModelODE, nCmt, time, amt, rate, ii, evid, cmt, addl, ss, theta         , biovar[1]      , tlag_data[1], rtol, atol, max_step);
  x = mixOde2CptModel_rk45(twoCptIndirectModelODE, nCmt, time, amt, rate, ii, evid, cmt, addl, ss, theta         , biovar         , tlag_data[1], rtol, atol, max_step);

  /* theta: data                                                            , biovar: data   , tlag: param */
  x = mixOde2CptModel_rk45(twoCptIndirectModelODE, nCmt, time, amt, rate, ii, evid, cmt, addl, ss, theta_data    , biovar_data    , tlag,         rtol, atol, max_step);
  x = mixOde2CptModel_rk45(twoCptIndirectModelODE, nCmt, time, amt, rate, ii, evid, cmt, addl, ss, theta_data[1] , biovar_data    , tlag,         rtol, atol, max_step);  
  x = mixOde2CptModel_rk45(twoCptIndirectModelODE, nCmt, time, amt, rate, ii, evid, cmt, addl, ss, theta_data[1] , biovar_data[1] , tlag,         rtol, atol, max_step);
  x = mixOde2CptModel_rk45(twoCptIndirectModelODE, nCmt, time, amt, rate, ii, evid, cmt, addl, ss, theta_data[1] , biovar_data[1] , tlag[1],      rtol, atol, max_step);
  x = mixOde2CptModel_rk45(twoCptIndirectModelODE, nCmt, time, amt, rate, ii, evid, cmt, addl, ss, theta_data[1] , biovar_data    , tlag[1],      rtol, atol, max_step);
  x = mixOde2CptModel_rk45(twoCptIndirectModelODE, nCmt, time, amt, rate, ii, evid, cmt, addl, ss, theta_data    , biovar_data[1] , tlag,         rtol, atol, max_step);
  x = mixOde2CptModel_rk45(twoCptIndirectModelODE, nCmt, time, amt, rate, ii, evid, cmt, addl, ss, theta_data    , biovar_data[1] , tlag[1],      rtol, atol, max_step);
  x = mixOde2CptModel_rk45(twoCptIndirectModelODE, nCmt, time, amt, rate, ii, evid, cmt, addl, ss, theta_data    , biovar_data    , tlag[1],      rtol, atol, max_step);

  /* theta: data                                                            , biovar: param  , tlag: data */
  x = mixOde2CptModel_rk45(twoCptIndirectModelODE, nCmt, time, amt, rate, ii, evid, cmt, addl, ss, theta_data    , biovar         , tlag_data,    rtol, atol, max_step);
  x = mixOde2CptModel_rk45(twoCptIndirectModelODE, nCmt, time, amt, rate, ii, evid, cmt, addl, ss, theta_data[1] , biovar         , tlag_data,    rtol, atol, max_step);  
  x = mixOde2CptModel_rk45(twoCptIndirectModelODE, nCmt, time, amt, rate, ii, evid, cmt, addl, ss, theta_data[1] , biovar[1]      , tlag_data,    rtol, atol, max_step);
  x = mixOde2CptModel_rk45(twoCptIndirectModelODE, nCmt, time, amt, rate, ii, evid, cmt, addl, ss, theta_data[1] , biovar[1]      , tlag_data[1], rtol, atol, max_step);
  x = mixOde2CptModel_rk45(twoCptIndirectModelODE, nCmt, time, amt, rate, ii, evid, cmt, addl, ss, theta_data[1] , biovar         , tlag_data[1], rtol, atol, max_step);
  x = mixOde2CptModel_rk45(twoCptIndirectModelODE, nCmt, time, amt, rate, ii, evid, cmt, addl, ss, theta_data    , biovar[1]      , tlag_data,    rtol, atol, max_step);
  x = mixOde2CptModel_rk45(twoCptIndirectModelODE, nCmt, time, amt, rate, ii, evid, cmt, addl, ss, theta_data    , biovar[1]      , tlag_data[1], rtol, atol, max_step);
  x = mixOde2CptModel_rk45(twoCptIndirectModelODE, nCmt, time, amt, rate, ii, evid, cmt, addl, ss, theta_data    , biovar         , tlag_data[1], rtol, atol, max_step);

  /* theta: param                                                           , biovar: data   , tlag: data */
  x = mixOde2CptModel_rk45(twoCptIndirectModelODE, nCmt, time, amt, rate, ii, evid, cmt, addl, ss, theta         , biovar_data    , tlag_data,    rtol, atol, max_step);
  x = mixOde2CptModel_rk45(twoCptIndirectModelODE, nCmt, time, amt, rate, ii, evid, cmt, addl, ss, theta[1]      , biovar_data    , tlag_data,    rtol, atol, max_step);  
  x = mixOde2CptModel_rk45(twoCptIndirectModelODE, nCmt, time, amt, rate, ii, evid, cmt, addl, ss, theta[1]      , biovar_data[1] , tlag_data,    rtol, atol, max_step);
  x = mixOde2CptModel_rk45(twoCptIndirectModelODE, nCmt, time, amt, rate, ii, evid, cmt, addl, ss, theta[1]      , biovar_data[1] , tlag_data[1], rtol, atol, max_step);
  x = mixOde2CptModel_rk45(twoCptIndirectModelODE, nCmt, time, amt, rate, ii, evid, cmt, addl, ss, theta[1]      , biovar_data    , tlag_data[1], rtol, atol, max_step);
  x = mixOde2CptModel_rk45(twoCptIndirectModelODE, nCmt, time, amt, rate, ii, evid, cmt, addl, ss, theta         , biovar_data[1] , tlag_data,    rtol, atol, max_step);
  x = mixOde2CptModel_rk45(twoCptIndirectModelODE, nCmt, time, amt, rate, ii, evid, cmt, addl, ss, theta         , biovar_data[1] , tlag_data[1], rtol, atol, max_step);
  x = mixOde2CptModel_rk45(twoCptIndirectModelODE, nCmt, time, amt, rate, ii, evid, cmt, addl, ss, theta         , biovar_data    , tlag_data[1], rtol, atol, max_step);

  /* bdf */
  /* data args */
  x = mixOde2CptModel_bdf(twoCptIndirectModelODE, nCmt, time, amt, rate, ii, evid, cmt, addl, ss, theta_data    , biovar_data    , tlag_data,    rtol, atol, max_step);
  x = mixOde2CptModel_bdf(twoCptIndirectModelODE, nCmt, time, amt, rate, ii, evid, cmt, addl, ss, theta_data[1] , biovar_data    , tlag_data,    rtol, atol, max_step);
  x = mixOde2CptModel_bdf(twoCptIndirectModelODE, nCmt, time, amt, rate, ii, evid, cmt, addl, ss, theta_data[1] , biovar_data[1] , tlag_data,    rtol, atol, max_step);
  x = mixOde2CptModel_bdf(twoCptIndirectModelODE, nCmt, time, amt, rate, ii, evid, cmt, addl, ss, theta_data[1] , biovar_data[1] , tlag_data[1], rtol, atol, max_step);
  x = mixOde2CptModel_bdf(twoCptIndirectModelODE, nCmt, time, amt, rate, ii, evid, cmt, addl, ss, theta_data[1] , biovar_data    , tlag_data[1], rtol, atol, max_step);
  x = mixOde2CptModel_bdf(twoCptIndirectModelODE, nCmt, time, amt, rate, ii, evid, cmt, addl, ss, theta_data    , biovar_data[1] , tlag_data,    rtol, atol, max_step);
  x = mixOde2CptModel_bdf(twoCptIndirectModelODE, nCmt, time, amt, rate, ii, evid, cmt, addl, ss, theta_data    , biovar_data[1] , tlag_data[1], rtol, atol, max_step);
  x = mixOde2CptModel_bdf(twoCptIndirectModelODE, nCmt, time, amt, rate, ii, evid, cmt, addl, ss, theta_data    , biovar_data    , tlag_data,    rtol, atol, max_step);

  /* param args */
  x = mixOde2CptModel_bdf(twoCptIndirectModelODE, nCmt, time, amt, rate, ii, evid, cmt, addl, ss, theta         , biovar         , tlag,         rtol, atol, max_step);
  x = mixOde2CptModel_bdf(twoCptIndirectModelODE, nCmt, time, amt, rate, ii, evid, cmt, addl, ss, theta[1]      , biovar         , tlag,         rtol, atol, max_step);  
  x = mixOde2CptModel_bdf(twoCptIndirectModelODE, nCmt, time, amt, rate, ii, evid, cmt, addl, ss, theta[1]      , biovar[1]      , tlag,         rtol, atol, max_step);
  x = mixOde2CptModel_bdf(twoCptIndirectModelODE, nCmt, time, amt, rate, ii, evid, cmt, addl, ss, theta[1]      , biovar[1]      , tlag[1],      rtol, atol, max_step);
  x = mixOde2CptModel_bdf(twoCptIndirectModelODE, nCmt, time, amt, rate, ii, evid, cmt, addl, ss, theta[1]      , biovar         , tlag[1],      rtol, atol, max_step);
  x = mixOde2CptModel_bdf(twoCptIndirectModelODE, nCmt, time, amt, rate, ii, evid, cmt, addl, ss, theta         , biovar[1]      , tlag,         rtol, atol, max_step);
  x = mixOde2CptModel_bdf(twoCptIndirectModelODE, nCmt, time, amt, rate, ii, evid, cmt, addl, ss, theta         , biovar[1]      , tlag[1],      rtol, atol, max_step);
  x = mixOde2CptModel_bdf(twoCptIndirectModelODE, nCmt, time, amt, rate, ii, evid, cmt, addl, ss, theta         , biovar         , tlag[1],      rtol, atol, max_step);

  /* theta: data                                                            , biovar: param  , tlag: param */
  x = mixOde2CptModel_bdf(twoCptIndirectModelODE, nCmt, time, amt, rate, ii, evid, cmt, addl, ss, theta_data    , biovar         , tlag,         rtol, atol, max_step);
  x = mixOde2CptModel_bdf(twoCptIndirectModelODE, nCmt, time, amt, rate, ii, evid, cmt, addl, ss, theta_data[1] , biovar         , tlag,         rtol, atol, max_step);  
  x = mixOde2CptModel_bdf(twoCptIndirectModelODE, nCmt, time, amt, rate, ii, evid, cmt, addl, ss, theta_data[1] , biovar[1]      , tlag,         rtol, atol, max_step);
  x = mixOde2CptModel_bdf(twoCptIndirectModelODE, nCmt, time, amt, rate, ii, evid, cmt, addl, ss, theta_data[1] , biovar[1]      , tlag[1],      rtol, atol, max_step);
  x = mixOde2CptModel_bdf(twoCptIndirectModelODE, nCmt, time, amt, rate, ii, evid, cmt, addl, ss, theta_data[1] , biovar         , tlag[1],      rtol, atol, max_step);
  x = mixOde2CptModel_bdf(twoCptIndirectModelODE, nCmt, time, amt, rate, ii, evid, cmt, addl, ss, theta_data    , biovar[1]      , tlag,         rtol, atol, max_step);
  x = mixOde2CptModel_bdf(twoCptIndirectModelODE, nCmt, time, amt, rate, ii, evid, cmt, addl, ss, theta_data    , biovar[1]      , tlag[1],      rtol, atol, max_step);
  x = mixOde2CptModel_bdf(twoCptIndirectModelODE, nCmt, time, amt, rate, ii, evid, cmt, addl, ss, theta_data    , biovar         , tlag[1],      rtol, atol, max_step);

  /* theta: parm                                                            , biovar: data   , tlag: param */
  x = mixOde2CptModel_bdf(twoCptIndirectModelODE, nCmt, time, amt, rate, ii, evid, cmt, addl, ss, theta         , biovar_data    , tlag,         rtol, atol, max_step);
  x = mixOde2CptModel_bdf(twoCptIndirectModelODE, nCmt, time, amt, rate, ii, evid, cmt, addl, ss, theta[1]      , biovar_data    , tlag,         rtol, atol, max_step);  
  x = mixOde2CptModel_bdf(twoCptIndirectModelODE, nCmt, time, amt, rate, ii, evid, cmt, addl, ss, theta[1]      , biovar_data[1] , tlag,         rtol, atol, max_step);
  x = mixOde2CptModel_bdf(twoCptIndirectModelODE, nCmt, time, amt, rate, ii, evid, cmt, addl, ss, theta[1]      , biovar_data[1] , tlag[1],      rtol, atol, max_step);
  x = mixOde2CptModel_bdf(twoCptIndirectModelODE, nCmt, time, amt, rate, ii, evid, cmt, addl, ss, theta[1]      , biovar_data    , tlag[1],      rtol, atol, max_step);
  x = mixOde2CptModel_bdf(twoCptIndirectModelODE, nCmt, time, amt, rate, ii, evid, cmt, addl, ss, theta         , biovar_data[1] , tlag,         rtol, atol, max_step);
  x = mixOde2CptModel_bdf(twoCptIndirectModelODE, nCmt, time, amt, rate, ii, evid, cmt, addl, ss, theta         , biovar_data[1] , tlag[1],      rtol, atol, max_step);
  x = mixOde2CptModel_bdf(twoCptIndirectModelODE, nCmt, time, amt, rate, ii, evid, cmt, addl, ss, theta         , biovar_data    , tlag[1],      rtol, atol, max_step);

  /* theta: parm                                                            , biovar: parm   , tlag: data */
  x = mixOde2CptModel_bdf(twoCptIndirectModelODE, nCmt, time, amt, rate, ii, evid, cmt, addl, ss, theta         , biovar         , tlag_data,    rtol, atol, max_step);
  x = mixOde2CptModel_bdf(twoCptIndirectModelODE, nCmt, time, amt, rate, ii, evid, cmt, addl, ss, theta[1]      , biovar         , tlag_data,    rtol, atol, max_step);  
  x = mixOde2CptModel_bdf(twoCptIndirectModelODE, nCmt, time, amt, rate, ii, evid, cmt, addl, ss, theta[1]      , biovar[1]      , tlag_data,    rtol, atol, max_step);
  x = mixOde2CptModel_bdf(twoCptIndirectModelODE, nCmt, time, amt, rate, ii, evid, cmt, addl, ss, theta[1]      , biovar[1]      , tlag_data[1], rtol, atol, max_step);
  x = mixOde2CptModel_bdf(twoCptIndirectModelODE, nCmt, time, amt, rate, ii, evid, cmt, addl, ss, theta[1]      , biovar         , tlag_data[1], rtol, atol, max_step);
  x = mixOde2CptModel_bdf(twoCptIndirectModelODE, nCmt, time, amt, rate, ii, evid, cmt, addl, ss, theta         , biovar[1]      , tlag_data,    rtol, atol, max_step);
  x = mixOde2CptModel_bdf(twoCptIndirectModelODE, nCmt, time, amt, rate, ii, evid, cmt, addl, ss, theta         , biovar[1]      , tlag_data[1], rtol, atol, max_step);
  x = mixOde2CptModel_bdf(twoCptIndirectModelODE, nCmt, time, amt, rate, ii, evid, cmt, addl, ss, theta         , biovar         , tlag_data[1], rtol, atol, max_step);

  /* theta: data                                                            , biovar: data   , tlag: param */
  x = mixOde2CptModel_bdf(twoCptIndirectModelODE, nCmt, time, amt, rate, ii, evid, cmt, addl, ss, theta_data    , biovar_data    , tlag,         rtol, atol, max_step);
  x = mixOde2CptModel_bdf(twoCptIndirectModelODE, nCmt, time, amt, rate, ii, evid, cmt, addl, ss, theta_data[1] , biovar_data    , tlag,         rtol, atol, max_step);  
  x = mixOde2CptModel_bdf(twoCptIndirectModelODE, nCmt, time, amt, rate, ii, evid, cmt, addl, ss, theta_data[1] , biovar_data[1] , tlag,         rtol, atol, max_step);
  x = mixOde2CptModel_bdf(twoCptIndirectModelODE, nCmt, time, amt, rate, ii, evid, cmt, addl, ss, theta_data[1] , biovar_data[1] , tlag[1],      rtol, atol, max_step);
  x = mixOde2CptModel_bdf(twoCptIndirectModelODE, nCmt, time, amt, rate, ii, evid, cmt, addl, ss, theta_data[1] , biovar_data    , tlag[1],      rtol, atol, max_step);
  x = mixOde2CptModel_bdf(twoCptIndirectModelODE, nCmt, time, amt, rate, ii, evid, cmt, addl, ss, theta_data    , biovar_data[1] , tlag,         rtol, atol, max_step);
  x = mixOde2CptModel_bdf(twoCptIndirectModelODE, nCmt, time, amt, rate, ii, evid, cmt, addl, ss, theta_data    , biovar_data[1] , tlag[1],      rtol, atol, max_step);
  x = mixOde2CptModel_bdf(twoCptIndirectModelODE, nCmt, time, amt, rate, ii, evid, cmt, addl, ss, theta_data    , biovar_data    , tlag[1],      rtol, atol, max_step);

  /* theta: data                                                            , biovar: param  , tlag: data */
  x = mixOde2CptModel_bdf(twoCptIndirectModelODE, nCmt, time, amt, rate, ii, evid, cmt, addl, ss, theta_data    , biovar         , tlag_data,    rtol, atol, max_step);
  x = mixOde2CptModel_bdf(twoCptIndirectModelODE, nCmt, time, amt, rate, ii, evid, cmt, addl, ss, theta_data[1] , biovar         , tlag_data,    rtol, atol, max_step);  
  x = mixOde2CptModel_bdf(twoCptIndirectModelODE, nCmt, time, amt, rate, ii, evid, cmt, addl, ss, theta_data[1] , biovar[1]      , tlag_data,    rtol, atol, max_step);
  x = mixOde2CptModel_bdf(twoCptIndirectModelODE, nCmt, time, amt, rate, ii, evid, cmt, addl, ss, theta_data[1] , biovar[1]      , tlag_data[1], rtol, atol, max_step);
  x = mixOde2CptModel_bdf(twoCptIndirectModelODE, nCmt, time, amt, rate, ii, evid, cmt, addl, ss, theta_data[1] , biovar         , tlag_data[1], rtol, atol, max_step);
  x = mixOde2CptModel_bdf(twoCptIndirectModelODE, nCmt, time, amt, rate, ii, evid, cmt, addl, ss, theta_data    , biovar[1]      , tlag_data,    rtol, atol, max_step);
  x = mixOde2CptModel_bdf(twoCptIndirectModelODE, nCmt, time, amt, rate, ii, evid, cmt, addl, ss, theta_data    , biovar[1]      , tlag_data[1], rtol, atol, max_step);
  x = mixOde2CptModel_bdf(twoCptIndirectModelODE, nCmt, time, amt, rate, ii, evid, cmt, addl, ss, theta_data    , biovar         , tlag_data[1], rtol, atol, max_step);

  /* theta: param                                                           , biovar: data   , tlag: data */
  x = mixOde2CptModel_bdf(twoCptIndirectModelODE, nCmt, time, amt, rate, ii, evid, cmt, addl, ss, theta         , biovar_data    , tlag_data,    rtol, atol, max_step);
  x = mixOde2CptModel_bdf(twoCptIndirectModelODE, nCmt, time, amt, rate, ii, evid, cmt, addl, ss, theta[1]      , biovar_data    , tlag_data,    rtol, atol, max_step);  
  x = mixOde2CptModel_bdf(twoCptIndirectModelODE, nCmt, time, amt, rate, ii, evid, cmt, addl, ss, theta[1]      , biovar_data[1] , tlag_data,    rtol, atol, max_step);
  x = mixOde2CptModel_bdf(twoCptIndirectModelODE, nCmt, time, amt, rate, ii, evid, cmt, addl, ss, theta[1]      , biovar_data[1] , tlag_data[1], rtol, atol, max_step);
  x = mixOde2CptModel_bdf(twoCptIndirectModelODE, nCmt, time, amt, rate, ii, evid, cmt, addl, ss, theta[1]      , biovar_data    , tlag_data[1], rtol, atol, max_step);
  x = mixOde2CptModel_bdf(twoCptIndirectModelODE, nCmt, time, amt, rate, ii, evid, cmt, addl, ss, theta         , biovar_data[1] , tlag_data,    rtol, atol, max_step);
  x = mixOde2CptModel_bdf(twoCptIndirectModelODE, nCmt, time, amt, rate, ii, evid, cmt, addl, ss, theta         , biovar_data[1] , tlag_data[1], rtol, atol, max_step);
  x = mixOde2CptModel_bdf(twoCptIndirectModelODE, nCmt, time, amt, rate, ii, evid, cmt, addl, ss, theta         , biovar_data    , tlag_data[1], rtol, atol, max_step);

  /* mixode + 1cpt */
  /* data args */
  x = mixOde1CptModel_rk45(twoCptIndirectModelODE, nCmt, time, amt, rate, ii, evid, cmt, addl, ss, theta_data    , biovar_data    , tlag_data,    rtol, atol, max_step);
  x = mixOde1CptModel_rk45(twoCptIndirectModelODE, nCmt, time, amt, rate, ii, evid, cmt, addl, ss, theta_data[1] , biovar_data    , tlag_data,    rtol, atol, max_step);
  x = mixOde1CptModel_rk45(twoCptIndirectModelODE, nCmt, time, amt, rate, ii, evid, cmt, addl, ss, theta_data[1] , biovar_data[1] , tlag_data,    rtol, atol, max_step);
  x = mixOde1CptModel_rk45(twoCptIndirectModelODE, nCmt, time, amt, rate, ii, evid, cmt, addl, ss, theta_data[1] , biovar_data[1] , tlag_data[1], rtol, atol, max_step);
  x = mixOde1CptModel_rk45(twoCptIndirectModelODE, nCmt, time, amt, rate, ii, evid, cmt, addl, ss, theta_data[1] , biovar_data    , tlag_data[1], rtol, atol, max_step);
  x = mixOde1CptModel_rk45(twoCptIndirectModelODE, nCmt, time, amt, rate, ii, evid, cmt, addl, ss, theta_data    , biovar_data[1] , tlag_data,    rtol, atol, max_step);
  x = mixOde1CptModel_rk45(twoCptIndirectModelODE, nCmt, time, amt, rate, ii, evid, cmt, addl, ss, theta_data    , biovar_data[1] , tlag_data[1], rtol, atol, max_step);
  x = mixOde1CptModel_rk45(twoCptIndirectModelODE, nCmt, time, amt, rate, ii, evid, cmt, addl, ss, theta_data    , biovar_data    , tlag_data,    rtol, atol, max_step);

  /* param args */
  x = mixOde1CptModel_rk45(twoCptIndirectModelODE, nCmt, time, amt, rate, ii, evid, cmt, addl, ss, theta         , biovar         , tlag,         rtol, atol, max_step);
  x = mixOde1CptModel_rk45(twoCptIndirectModelODE, nCmt, time, amt, rate, ii, evid, cmt, addl, ss, theta[1]      , biovar         , tlag,         rtol, atol, max_step);  
  x = mixOde1CptModel_rk45(twoCptIndirectModelODE, nCmt, time, amt, rate, ii, evid, cmt, addl, ss, theta[1]      , biovar[1]      , tlag,         rtol, atol, max_step);
  x = mixOde1CptModel_rk45(twoCptIndirectModelODE, nCmt, time, amt, rate, ii, evid, cmt, addl, ss, theta[1]      , biovar[1]      , tlag[1],      rtol, atol, max_step);
  x = mixOde1CptModel_rk45(twoCptIndirectModelODE, nCmt, time, amt, rate, ii, evid, cmt, addl, ss, theta[1]      , biovar         , tlag[1],      rtol, atol, max_step);
  x = mixOde1CptModel_rk45(twoCptIndirectModelODE, nCmt, time, amt, rate, ii, evid, cmt, addl, ss, theta         , biovar[1]      , tlag,         rtol, atol, max_step);
  x = mixOde1CptModel_rk45(twoCptIndirectModelODE, nCmt, time, amt, rate, ii, evid, cmt, addl, ss, theta         , biovar[1]      , tlag[1],      rtol, atol, max_step);
  x = mixOde1CptModel_rk45(twoCptIndirectModelODE, nCmt, time, amt, rate, ii, evid, cmt, addl, ss, theta         , biovar         , tlag[1],      rtol, atol, max_step);

  /* theta: data                                                            , biovar: param  , tlag: param */
  x = mixOde1CptModel_rk45(twoCptIndirectModelODE, nCmt, time, amt, rate, ii, evid, cmt, addl, ss, theta_data    , biovar         , tlag,         rtol, atol, max_step);
  x = mixOde1CptModel_rk45(twoCptIndirectModelODE, nCmt, time, amt, rate, ii, evid, cmt, addl, ss, theta_data[1] , biovar         , tlag,         rtol, atol, max_step);  
  x = mixOde1CptModel_rk45(twoCptIndirectModelODE, nCmt, time, amt, rate, ii, evid, cmt, addl, ss, theta_data[1] , biovar[1]      , tlag,         rtol, atol, max_step);
  x = mixOde1CptModel_rk45(twoCptIndirectModelODE, nCmt, time, amt, rate, ii, evid, cmt, addl, ss, theta_data[1] , biovar[1]      , tlag[1],      rtol, atol, max_step);
  x = mixOde1CptModel_rk45(twoCptIndirectModelODE, nCmt, time, amt, rate, ii, evid, cmt, addl, ss, theta_data[1] , biovar         , tlag[1],      rtol, atol, max_step);
  x = mixOde1CptModel_rk45(twoCptIndirectModelODE, nCmt, time, amt, rate, ii, evid, cmt, addl, ss, theta_data    , biovar[1]      , tlag,         rtol, atol, max_step);
  x = mixOde1CptModel_rk45(twoCptIndirectModelODE, nCmt, time, amt, rate, ii, evid, cmt, addl, ss, theta_data    , biovar[1]      , tlag[1],      rtol, atol, max_step);
  x = mixOde1CptModel_rk45(twoCptIndirectModelODE, nCmt, time, amt, rate, ii, evid, cmt, addl, ss, theta_data    , biovar         , tlag[1],      rtol, atol, max_step);

  /* theta: parm                                                            , biovar: data   , tlag: param */
  x = mixOde1CptModel_rk45(twoCptIndirectModelODE, nCmt, time, amt, rate, ii, evid, cmt, addl, ss, theta         , biovar_data    , tlag,         rtol, atol, max_step);
  x = mixOde1CptModel_rk45(twoCptIndirectModelODE, nCmt, time, amt, rate, ii, evid, cmt, addl, ss, theta[1]      , biovar_data    , tlag,         rtol, atol, max_step);  
  x = mixOde1CptModel_rk45(twoCptIndirectModelODE, nCmt, time, amt, rate, ii, evid, cmt, addl, ss, theta[1]      , biovar_data[1] , tlag,         rtol, atol, max_step);
  x = mixOde1CptModel_rk45(twoCptIndirectModelODE, nCmt, time, amt, rate, ii, evid, cmt, addl, ss, theta[1]      , biovar_data[1] , tlag[1],      rtol, atol, max_step);
  x = mixOde1CptModel_rk45(twoCptIndirectModelODE, nCmt, time, amt, rate, ii, evid, cmt, addl, ss, theta[1]      , biovar_data    , tlag[1],      rtol, atol, max_step);
  x = mixOde1CptModel_rk45(twoCptIndirectModelODE, nCmt, time, amt, rate, ii, evid, cmt, addl, ss, theta         , biovar_data[1] , tlag,         rtol, atol, max_step);
  x = mixOde1CptModel_rk45(twoCptIndirectModelODE, nCmt, time, amt, rate, ii, evid, cmt, addl, ss, theta         , biovar_data[1] , tlag[1],      rtol, atol, max_step);
  x = mixOde1CptModel_rk45(twoCptIndirectModelODE, nCmt, time, amt, rate, ii, evid, cmt, addl, ss, theta         , biovar_data    , tlag[1],      rtol, atol, max_step);

  /* theta: parm                                                            , biovar: parm   , tlag: data */
  x = mixOde1CptModel_rk45(twoCptIndirectModelODE, nCmt, time, amt, rate, ii, evid, cmt, addl, ss, theta         , biovar         , tlag_data,    rtol, atol, max_step);
  x = mixOde1CptModel_rk45(twoCptIndirectModelODE, nCmt, time, amt, rate, ii, evid, cmt, addl, ss, theta[1]      , biovar         , tlag_data,    rtol, atol, max_step);  
  x = mixOde1CptModel_rk45(twoCptIndirectModelODE, nCmt, time, amt, rate, ii, evid, cmt, addl, ss, theta[1]      , biovar[1]      , tlag_data,    rtol, atol, max_step);
  x = mixOde1CptModel_rk45(twoCptIndirectModelODE, nCmt, time, amt, rate, ii, evid, cmt, addl, ss, theta[1]      , biovar[1]      , tlag_data[1], rtol, atol, max_step);
  x = mixOde1CptModel_rk45(twoCptIndirectModelODE, nCmt, time, amt, rate, ii, evid, cmt, addl, ss, theta[1]      , biovar         , tlag_data[1], rtol, atol, max_step);
  x = mixOde1CptModel_rk45(twoCptIndirectModelODE, nCmt, time, amt, rate, ii, evid, cmt, addl, ss, theta         , biovar[1]      , tlag_data,    rtol, atol, max_step);
  x = mixOde1CptModel_rk45(twoCptIndirectModelODE, nCmt, time, amt, rate, ii, evid, cmt, addl, ss, theta         , biovar[1]      , tlag_data[1], rtol, atol, max_step);
  x = mixOde1CptModel_rk45(twoCptIndirectModelODE, nCmt, time, amt, rate, ii, evid, cmt, addl, ss, theta         , biovar         , tlag_data[1], rtol, atol, max_step);

  /* theta: data                                                            , biovar: data   , tlag: param */
  x = mixOde1CptModel_rk45(twoCptIndirectModelODE, nCmt, time, amt, rate, ii, evid, cmt, addl, ss, theta_data    , biovar_data    , tlag,         rtol, atol, max_step);
  x = mixOde1CptModel_rk45(twoCptIndirectModelODE, nCmt, time, amt, rate, ii, evid, cmt, addl, ss, theta_data[1] , biovar_data    , tlag,         rtol, atol, max_step);  
  x = mixOde1CptModel_rk45(twoCptIndirectModelODE, nCmt, time, amt, rate, ii, evid, cmt, addl, ss, theta_data[1] , biovar_data[1] , tlag,         rtol, atol, max_step);
  x = mixOde1CptModel_rk45(twoCptIndirectModelODE, nCmt, time, amt, rate, ii, evid, cmt, addl, ss, theta_data[1] , biovar_data[1] , tlag[1],      rtol, atol, max_step);
  x = mixOde1CptModel_rk45(twoCptIndirectModelODE, nCmt, time, amt, rate, ii, evid, cmt, addl, ss, theta_data[1] , biovar_data    , tlag[1],      rtol, atol, max_step);
  x = mixOde1CptModel_rk45(twoCptIndirectModelODE, nCmt, time, amt, rate, ii, evid, cmt, addl, ss, theta_data    , biovar_data[1] , tlag,         rtol, atol, max_step);
  x = mixOde1CptModel_rk45(twoCptIndirectModelODE, nCmt, time, amt, rate, ii, evid, cmt, addl, ss, theta_data    , biovar_data[1] , tlag[1],      rtol, atol, max_step);
  x = mixOde1CptModel_rk45(twoCptIndirectModelODE, nCmt, time, amt, rate, ii, evid, cmt, addl, ss, theta_data    , biovar_data    , tlag[1],      rtol, atol, max_step);

  /* theta: data                                                            , biovar: param  , tlag: data */
  x = mixOde1CptModel_rk45(twoCptIndirectModelODE, nCmt, time, amt, rate, ii, evid, cmt, addl, ss, theta_data    , biovar         , tlag_data,    rtol, atol, max_step);
  x = mixOde1CptModel_rk45(twoCptIndirectModelODE, nCmt, time, amt, rate, ii, evid, cmt, addl, ss, theta_data[1] , biovar         , tlag_data,    rtol, atol, max_step);  
  x = mixOde1CptModel_rk45(twoCptIndirectModelODE, nCmt, time, amt, rate, ii, evid, cmt, addl, ss, theta_data[1] , biovar[1]      , tlag_data,    rtol, atol, max_step);
  x = mixOde1CptModel_rk45(twoCptIndirectModelODE, nCmt, time, amt, rate, ii, evid, cmt, addl, ss, theta_data[1] , biovar[1]      , tlag_data[1], rtol, atol, max_step);
  x = mixOde1CptModel_rk45(twoCptIndirectModelODE, nCmt, time, amt, rate, ii, evid, cmt, addl, ss, theta_data[1] , biovar         , tlag_data[1], rtol, atol, max_step);
  x = mixOde1CptModel_rk45(twoCptIndirectModelODE, nCmt, time, amt, rate, ii, evid, cmt, addl, ss, theta_data    , biovar[1]      , tlag_data,    rtol, atol, max_step);
  x = mixOde1CptModel_rk45(twoCptIndirectModelODE, nCmt, time, amt, rate, ii, evid, cmt, addl, ss, theta_data    , biovar[1]      , tlag_data[1], rtol, atol, max_step);
  x = mixOde1CptModel_rk45(twoCptIndirectModelODE, nCmt, time, amt, rate, ii, evid, cmt, addl, ss, theta_data    , biovar         , tlag_data[1], rtol, atol, max_step);

  /* theta: param                                                           , biovar: data   , tlag: data */
  x = mixOde1CptModel_rk45(twoCptIndirectModelODE, nCmt, time, amt, rate, ii, evid, cmt, addl, ss, theta         , biovar_data    , tlag_data,    rtol, atol, max_step);
  x = mixOde1CptModel_rk45(twoCptIndirectModelODE, nCmt, time, amt, rate, ii, evid, cmt, addl, ss, theta[1]      , biovar_data    , tlag_data,    rtol, atol, max_step);  
  x = mixOde1CptModel_rk45(twoCptIndirectModelODE, nCmt, time, amt, rate, ii, evid, cmt, addl, ss, theta[1]      , biovar_data[1] , tlag_data,    rtol, atol, max_step);
  x = mixOde1CptModel_rk45(twoCptIndirectModelODE, nCmt, time, amt, rate, ii, evid, cmt, addl, ss, theta[1]      , biovar_data[1] , tlag_data[1], rtol, atol, max_step);
  x = mixOde1CptModel_rk45(twoCptIndirectModelODE, nCmt, time, amt, rate, ii, evid, cmt, addl, ss, theta[1]      , biovar_data    , tlag_data[1], rtol, atol, max_step);
  x = mixOde1CptModel_rk45(twoCptIndirectModelODE, nCmt, time, amt, rate, ii, evid, cmt, addl, ss, theta         , biovar_data[1] , tlag_data,    rtol, atol, max_step);
  x = mixOde1CptModel_rk45(twoCptIndirectModelODE, nCmt, time, amt, rate, ii, evid, cmt, addl, ss, theta         , biovar_data[1] , tlag_data[1], rtol, atol, max_step);
  x = mixOde1CptModel_rk45(twoCptIndirectModelODE, nCmt, time, amt, rate, ii, evid, cmt, addl, ss, theta         , biovar_data    , tlag_data[1], rtol, atol, max_step);

  /* bdf */
  /* data args */
  x = mixOde1CptModel_bdf(twoCptIndirectModelODE, nCmt, time, amt, rate, ii, evid, cmt, addl, ss, theta_data    , biovar_data    , tlag_data,    rtol, atol, max_step);
  x = mixOde1CptModel_bdf(twoCptIndirectModelODE, nCmt, time, amt, rate, ii, evid, cmt, addl, ss, theta_data[1] , biovar_data    , tlag_data,    rtol, atol, max_step);
  x = mixOde1CptModel_bdf(twoCptIndirectModelODE, nCmt, time, amt, rate, ii, evid, cmt, addl, ss, theta_data[1] , biovar_data[1] , tlag_data,    rtol, atol, max_step);
  x = mixOde1CptModel_bdf(twoCptIndirectModelODE, nCmt, time, amt, rate, ii, evid, cmt, addl, ss, theta_data[1] , biovar_data[1] , tlag_data[1], rtol, atol, max_step);
  x = mixOde1CptModel_bdf(twoCptIndirectModelODE, nCmt, time, amt, rate, ii, evid, cmt, addl, ss, theta_data[1] , biovar_data    , tlag_data[1], rtol, atol, max_step);
  x = mixOde1CptModel_bdf(twoCptIndirectModelODE, nCmt, time, amt, rate, ii, evid, cmt, addl, ss, theta_data    , biovar_data[1] , tlag_data,    rtol, atol, max_step);
  x = mixOde1CptModel_bdf(twoCptIndirectModelODE, nCmt, time, amt, rate, ii, evid, cmt, addl, ss, theta_data    , biovar_data[1] , tlag_data[1], rtol, atol, max_step);
  x = mixOde1CptModel_bdf(twoCptIndirectModelODE, nCmt, time, amt, rate, ii, evid, cmt, addl, ss, theta_data    , biovar_data    , tlag_data,    rtol, atol, max_step);

  /* param args */
  x = mixOde1CptModel_bdf(twoCptIndirectModelODE, nCmt, time, amt, rate, ii, evid, cmt, addl, ss, theta         , biovar         , tlag,         rtol, atol, max_step);
  x = mixOde1CptModel_bdf(twoCptIndirectModelODE, nCmt, time, amt, rate, ii, evid, cmt, addl, ss, theta[1]      , biovar         , tlag,         rtol, atol, max_step);  
  x = mixOde1CptModel_bdf(twoCptIndirectModelODE, nCmt, time, amt, rate, ii, evid, cmt, addl, ss, theta[1]      , biovar[1]      , tlag,         rtol, atol, max_step);
  x = mixOde1CptModel_bdf(twoCptIndirectModelODE, nCmt, time, amt, rate, ii, evid, cmt, addl, ss, theta[1]      , biovar[1]      , tlag[1],      rtol, atol, max_step);
  x = mixOde1CptModel_bdf(twoCptIndirectModelODE, nCmt, time, amt, rate, ii, evid, cmt, addl, ss, theta[1]      , biovar         , tlag[1],      rtol, atol, max_step);
  x = mixOde1CptModel_bdf(twoCptIndirectModelODE, nCmt, time, amt, rate, ii, evid, cmt, addl, ss, theta         , biovar[1]      , tlag,         rtol, atol, max_step);
  x = mixOde1CptModel_bdf(twoCptIndirectModelODE, nCmt, time, amt, rate, ii, evid, cmt, addl, ss, theta         , biovar[1]      , tlag[1],      rtol, atol, max_step);
  x = mixOde1CptModel_bdf(twoCptIndirectModelODE, nCmt, time, amt, rate, ii, evid, cmt, addl, ss, theta         , biovar         , tlag[1],      rtol, atol, max_step);

  /* theta: data                                                            , biovar: param  , tlag: param */
  x = mixOde1CptModel_bdf(twoCptIndirectModelODE, nCmt, time, amt, rate, ii, evid, cmt, addl, ss, theta_data    , biovar         , tlag,         rtol, atol, max_step);
  x = mixOde1CptModel_bdf(twoCptIndirectModelODE, nCmt, time, amt, rate, ii, evid, cmt, addl, ss, theta_data[1] , biovar         , tlag,         rtol, atol, max_step);  
  x = mixOde1CptModel_bdf(twoCptIndirectModelODE, nCmt, time, amt, rate, ii, evid, cmt, addl, ss, theta_data[1] , biovar[1]      , tlag,         rtol, atol, max_step);
  x = mixOde1CptModel_bdf(twoCptIndirectModelODE, nCmt, time, amt, rate, ii, evid, cmt, addl, ss, theta_data[1] , biovar[1]      , tlag[1],      rtol, atol, max_step);
  x = mixOde1CptModel_bdf(twoCptIndirectModelODE, nCmt, time, amt, rate, ii, evid, cmt, addl, ss, theta_data[1] , biovar         , tlag[1],      rtol, atol, max_step);
  x = mixOde1CptModel_bdf(twoCptIndirectModelODE, nCmt, time, amt, rate, ii, evid, cmt, addl, ss, theta_data    , biovar[1]      , tlag,         rtol, atol, max_step);
  x = mixOde1CptModel_bdf(twoCptIndirectModelODE, nCmt, time, amt, rate, ii, evid, cmt, addl, ss, theta_data    , biovar[1]      , tlag[1],      rtol, atol, max_step);
  x = mixOde1CptModel_bdf(twoCptIndirectModelODE, nCmt, time, amt, rate, ii, evid, cmt, addl, ss, theta_data    , biovar         , tlag[1],      rtol, atol, max_step);

  /* theta: parm                                                            , biovar: data   , tlag: param */
  x = mixOde1CptModel_bdf(twoCptIndirectModelODE, nCmt, time, amt, rate, ii, evid, cmt, addl, ss, theta         , biovar_data    , tlag,         rtol, atol, max_step);
  x = mixOde1CptModel_bdf(twoCptIndirectModelODE, nCmt, time, amt, rate, ii, evid, cmt, addl, ss, theta[1]      , biovar_data    , tlag,         rtol, atol, max_step);  
  x = mixOde1CptModel_bdf(twoCptIndirectModelODE, nCmt, time, amt, rate, ii, evid, cmt, addl, ss, theta[1]      , biovar_data[1] , tlag,         rtol, atol, max_step);
  x = mixOde1CptModel_bdf(twoCptIndirectModelODE, nCmt, time, amt, rate, ii, evid, cmt, addl, ss, theta[1]      , biovar_data[1] , tlag[1],      rtol, atol, max_step);
  x = mixOde1CptModel_bdf(twoCptIndirectModelODE, nCmt, time, amt, rate, ii, evid, cmt, addl, ss, theta[1]      , biovar_data    , tlag[1],      rtol, atol, max_step);
  x = mixOde1CptModel_bdf(twoCptIndirectModelODE, nCmt, time, amt, rate, ii, evid, cmt, addl, ss, theta         , biovar_data[1] , tlag,         rtol, atol, max_step);
  x = mixOde1CptModel_bdf(twoCptIndirectModelODE, nCmt, time, amt, rate, ii, evid, cmt, addl, ss, theta         , biovar_data[1] , tlag[1],      rtol, atol, max_step);
  x = mixOde1CptModel_bdf(twoCptIndirectModelODE, nCmt, time, amt, rate, ii, evid, cmt, addl, ss, theta         , biovar_data    , tlag[1],      rtol, atol, max_step);

  /* theta: parm                                                            , biovar: parm   , tlag: data */
  x = mixOde1CptModel_bdf(twoCptIndirectModelODE, nCmt, time, amt, rate, ii, evid, cmt, addl, ss, theta         , biovar         , tlag_data,    rtol, atol, max_step);
  x = mixOde1CptModel_bdf(twoCptIndirectModelODE, nCmt, time, amt, rate, ii, evid, cmt, addl, ss, theta[1]      , biovar         , tlag_data,    rtol, atol, max_step);  
  x = mixOde1CptModel_bdf(twoCptIndirectModelODE, nCmt, time, amt, rate, ii, evid, cmt, addl, ss, theta[1]      , biovar[1]      , tlag_data,    rtol, atol, max_step);
  x = mixOde1CptModel_bdf(twoCptIndirectModelODE, nCmt, time, amt, rate, ii, evid, cmt, addl, ss, theta[1]      , biovar[1]      , tlag_data[1], rtol, atol, max_step);
  x = mixOde1CptModel_bdf(twoCptIndirectModelODE, nCmt, time, amt, rate, ii, evid, cmt, addl, ss, theta[1]      , biovar         , tlag_data[1], rtol, atol, max_step);
  x = mixOde1CptModel_bdf(twoCptIndirectModelODE, nCmt, time, amt, rate, ii, evid, cmt, addl, ss, theta         , biovar[1]      , tlag_data,    rtol, atol, max_step);
  x = mixOde1CptModel_bdf(twoCptIndirectModelODE, nCmt, time, amt, rate, ii, evid, cmt, addl, ss, theta         , biovar[1]      , tlag_data[1], rtol, atol, max_step);
  x = mixOde1CptModel_bdf(twoCptIndirectModelODE, nCmt, time, amt, rate, ii, evid, cmt, addl, ss, theta         , biovar         , tlag_data[1], rtol, atol, max_step);

  /* theta: data                                                            , biovar: data   , tlag: param */
  x = mixOde1CptModel_bdf(twoCptIndirectModelODE, nCmt, time, amt, rate, ii, evid, cmt, addl, ss, theta_data    , biovar_data    , tlag,         rtol, atol, max_step);
  x = mixOde1CptModel_bdf(twoCptIndirectModelODE, nCmt, time, amt, rate, ii, evid, cmt, addl, ss, theta_data[1] , biovar_data    , tlag,         rtol, atol, max_step);  
  x = mixOde1CptModel_bdf(twoCptIndirectModelODE, nCmt, time, amt, rate, ii, evid, cmt, addl, ss, theta_data[1] , biovar_data[1] , tlag,         rtol, atol, max_step);
  x = mixOde1CptModel_bdf(twoCptIndirectModelODE, nCmt, time, amt, rate, ii, evid, cmt, addl, ss, theta_data[1] , biovar_data[1] , tlag[1],      rtol, atol, max_step);
  x = mixOde1CptModel_bdf(twoCptIndirectModelODE, nCmt, time, amt, rate, ii, evid, cmt, addl, ss, theta_data[1] , biovar_data    , tlag[1],      rtol, atol, max_step);
  x = mixOde1CptModel_bdf(twoCptIndirectModelODE, nCmt, time, amt, rate, ii, evid, cmt, addl, ss, theta_data    , biovar_data[1] , tlag,         rtol, atol, max_step);
  x = mixOde1CptModel_bdf(twoCptIndirectModelODE, nCmt, time, amt, rate, ii, evid, cmt, addl, ss, theta_data    , biovar_data[1] , tlag[1],      rtol, atol, max_step);
  x = mixOde1CptModel_bdf(twoCptIndirectModelODE, nCmt, time, amt, rate, ii, evid, cmt, addl, ss, theta_data    , biovar_data    , tlag[1],      rtol, atol, max_step);

  /* theta: data                                                            , biovar: param  , tlag: data */
  x = mixOde1CptModel_bdf(twoCptIndirectModelODE, nCmt, time, amt, rate, ii, evid, cmt, addl, ss, theta_data    , biovar         , tlag_data,    rtol, atol, max_step);
  x = mixOde1CptModel_bdf(twoCptIndirectModelODE, nCmt, time, amt, rate, ii, evid, cmt, addl, ss, theta_data[1] , biovar         , tlag_data,    rtol, atol, max_step);  
  x = mixOde1CptModel_bdf(twoCptIndirectModelODE, nCmt, time, amt, rate, ii, evid, cmt, addl, ss, theta_data[1] , biovar[1]      , tlag_data,    rtol, atol, max_step);
  x = mixOde1CptModel_bdf(twoCptIndirectModelODE, nCmt, time, amt, rate, ii, evid, cmt, addl, ss, theta_data[1] , biovar[1]      , tlag_data[1], rtol, atol, max_step);
  x = mixOde1CptModel_bdf(twoCptIndirectModelODE, nCmt, time, amt, rate, ii, evid, cmt, addl, ss, theta_data[1] , biovar         , tlag_data[1], rtol, atol, max_step);
  x = mixOde1CptModel_bdf(twoCptIndirectModelODE, nCmt, time, amt, rate, ii, evid, cmt, addl, ss, theta_data    , biovar[1]      , tlag_data,    rtol, atol, max_step);
  x = mixOde1CptModel_bdf(twoCptIndirectModelODE, nCmt, time, amt, rate, ii, evid, cmt, addl, ss, theta_data    , biovar[1]      , tlag_data[1], rtol, atol, max_step);
  x = mixOde1CptModel_bdf(twoCptIndirectModelODE, nCmt, time, amt, rate, ii, evid, cmt, addl, ss, theta_data    , biovar         , tlag_data[1], rtol, atol, max_step);

  /* theta: param                                                           , biovar: data   , tlag: data */
  x = mixOde1CptModel_bdf(twoCptIndirectModelODE, nCmt, time, amt, rate, ii, evid, cmt, addl, ss, theta         , biovar_data    , tlag_data,    rtol, atol, max_step);
  x = mixOde1CptModel_bdf(twoCptIndirectModelODE, nCmt, time, amt, rate, ii, evid, cmt, addl, ss, theta[1]      , biovar_data    , tlag_data,    rtol, atol, max_step);  
  x = mixOde1CptModel_bdf(twoCptIndirectModelODE, nCmt, time, amt, rate, ii, evid, cmt, addl, ss, theta[1]      , biovar_data[1] , tlag_data,    rtol, atol, max_step);
  x = mixOde1CptModel_bdf(twoCptIndirectModelODE, nCmt, time, amt, rate, ii, evid, cmt, addl, ss, theta[1]      , biovar_data[1] , tlag_data[1], rtol, atol, max_step);
  x = mixOde1CptModel_bdf(twoCptIndirectModelODE, nCmt, time, amt, rate, ii, evid, cmt, addl, ss, theta[1]      , biovar_data    , tlag_data[1], rtol, atol, max_step);
  x = mixOde1CptModel_bdf(twoCptIndirectModelODE, nCmt, time, amt, rate, ii, evid, cmt, addl, ss, theta         , biovar_data[1] , tlag_data,    rtol, atol, max_step);
  x = mixOde1CptModel_bdf(twoCptIndirectModelODE, nCmt, time, amt, rate, ii, evid, cmt, addl, ss, theta         , biovar_data[1] , tlag_data[1], rtol, atol, max_step);
  x = mixOde1CptModel_bdf(twoCptIndirectModelODE, nCmt, time, amt, rate, ii, evid, cmt, addl, ss, theta         , biovar_data    , tlag_data[1], rtol, atol, max_step);
}

