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
  int nCmt = 3;  // number of compartments in model

  // Since we're not trying to evaluate the bio-variability (F) and 
  // the lag times, we declare them as data.
  real biovar_data[nt, nCmt];
  real tlag_data[nt, nCmt];

  real<lower = 0> CL_data;
  real<lower = 0> Q_data;
  real<lower = 0> V2_data;
  real<lower = 0> V3_data;
  real<lower = 0> ka_data;

  matrix[3, 3] K_data[nt];
  real k10_data;
  real k12_data;
  real k21_data;

  CL_data  = 5;
  ka_data  = 1.2;
  Q_data   = 8;
  V2_data  = 20;
  V3_data  = 70;
  k10_data = CL_data / V2_data;
  k12_data = Q_data  / V2_data;
  k21_data = Q_data  / V3_data;

  for (i in 1:nt) {
    K_data[i] = rep_matrix(0, 3, 3);
    K_data[i][1, 1] = -ka_data;
    K_data[i][2, 1] = ka_data;
    K_data[i][2, 2] = -(k10_data + k12_data);
    K_data[i][2, 3] = k21_data;
    K_data[i][3, 2] = k12_data;
    K_data[i][3, 3] = -k21_data;
  }

  for (i in 1:nt) {
    biovar_data[i] = { 1, 1, 1 };
    tlag_data[i] = { 0, 0, 0 };
  }

}

generated quantities {
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
  matrix[3, 3] K[nt];
  real k10;
  real k12;
  real k21;

  vector<lower = 0>[nObs] cHatDat1;
  vector<lower = 0>[nObs] cHatDat2;
  vector<lower = 0>[nObs] cHatDat3;
  vector<lower = 0>[nObs] cHatDat4;
  vector<lower = 0>[nObs] cHatDat5;
  vector<lower = 0>[nObs] cHatDat6;
  vector<lower = 0>[nObs] cHatDat7;
  vector<lower = 0>[nObs] cHatDat8;

  vector<lower = 0>[nObs] cHatObs1;
  vector<lower = 0>[nObs] cHatObs2;
  vector<lower = 0>[nObs] cHatObs3;
  vector<lower = 0>[nObs] cHatObs4;
  vector<lower = 0>[nObs] cHatObs5;
  vector<lower = 0>[nObs] cHatObs6;
  vector<lower = 0>[nObs] cHatObs7;
  vector<lower = 0>[nObs] cHatObs8;

  /* theta: data, biovar: param, tlag: param */
  vector<lower = 0>[nObs] cHatDPP1;
  vector<lower = 0>[nObs] cHatDPP2;
  vector<lower = 0>[nObs] cHatDPP3;
  vector<lower = 0>[nObs] cHatDPP4;
  vector<lower = 0>[nObs] cHatDPP5;
  vector<lower = 0>[nObs] cHatDPP6;
  vector<lower = 0>[nObs] cHatDPP7;
  vector<lower = 0>[nObs] cHatDPP8;

  /* theta: param, biovar: data, tlag: param */
  vector<lower = 0>[nObs] cHatPDP1;
  vector<lower = 0>[nObs] cHatPDP2;
  vector<lower = 0>[nObs] cHatPDP3;
  vector<lower = 0>[nObs] cHatPDP4;
  vector<lower = 0>[nObs] cHatPDP5;
  vector<lower = 0>[nObs] cHatPDP6;
  vector<lower = 0>[nObs] cHatPDP7;
  vector<lower = 0>[nObs] cHatPDP8;

  /* theta: param, biovar: param, tlag: data */
  vector<lower = 0>[nObs] cHatPPD1;
  vector<lower = 0>[nObs] cHatPPD2;
  vector<lower = 0>[nObs] cHatPPD3;
  vector<lower = 0>[nObs] cHatPPD4;
  vector<lower = 0>[nObs] cHatPPD5;
  vector<lower = 0>[nObs] cHatPPD6;
  vector<lower = 0>[nObs] cHatPPD7;
  vector<lower = 0>[nObs] cHatPPD8;

  /* theta: data, biovar: data, tlag: parm */
  vector<lower = 0>[nObs] cHatDDP1;
  vector<lower = 0>[nObs] cHatDDP2;
  vector<lower = 0>[nObs] cHatDDP3;
  vector<lower = 0>[nObs] cHatDDP4;
  vector<lower = 0>[nObs] cHatDDP5;
  vector<lower = 0>[nObs] cHatDDP6;
  vector<lower = 0>[nObs] cHatDDP7;
  vector<lower = 0>[nObs] cHatDDP8;

  /* theta: data, biovar: param, tlag: data */
  vector<lower = 0>[nObs] cHatDPD1;
  vector<lower = 0>[nObs] cHatDPD2;
  vector<lower = 0>[nObs] cHatDPD3;
  vector<lower = 0>[nObs] cHatDPD4;
  vector<lower = 0>[nObs] cHatDPD5;
  vector<lower = 0>[nObs] cHatDPD6;
  vector<lower = 0>[nObs] cHatDPD7;
  vector<lower = 0>[nObs] cHatDPD8;

  /* theta: param, biovar: data, tlag: data */
  vector<lower = 0>[nObs] cHatPDD1;
  vector<lower = 0>[nObs] cHatPDD2;
  vector<lower = 0>[nObs] cHatPDD3;
  vector<lower = 0>[nObs] cHatPDD4;
  vector<lower = 0>[nObs] cHatPDD5;
  vector<lower = 0>[nObs] cHatPDD6;
  vector<lower = 0>[nObs] cHatPDD7;
  vector<lower = 0>[nObs] cHatPDD8;

  CL    = 5;
  ka    = 1.2;
  Q     = 8;
  sigma = 0.01;
  V2    = 20; /* central cpt */
  V3    = 70;
  k10   = CL  / V2;
  k12   = Q   / V2;
  k21   = Q   / V3;

  for (i in 1:nt) {
    K[i] = rep_matrix(0, 3, 3);
    K[i][1, 1] = -ka;
    K[i][2, 1] = ka;
    K[i][2, 2] = -(k10 + k12);
    K[i][2, 3] = k21;
    K[i][3, 2] = k12;
    K[i][3, 3] = -k21;
  }

  for (i in 1:nt) {
    biovar[i] = { 1, 1, 1 };
    tlag[i] = { 0, 0, 0 };
  }

  /* data args */
  x = linOdeModel(time, amt, rate, ii, evid, cmt, addl, ss, K_data    , biovar_data    , tlag_data)    ;  cHat = col(x, 2) ./ V2;  cHatDat1 = cHat[iObs];
  x = linOdeModel(time, amt, rate, ii, evid, cmt, addl, ss, K_data[1] , biovar_data    , tlag_data)    ;  cHat = col(x, 2) ./ V2;  cHatDat2 = cHat[iObs];
  x = linOdeModel(time, amt, rate, ii, evid, cmt, addl, ss, K_data[1] , biovar_data[1] , tlag_data)    ;  cHat = col(x, 2) ./ V2;  cHatDat3 = cHat[iObs];
  x = linOdeModel(time, amt, rate, ii, evid, cmt, addl, ss, K_data[1] , biovar_data[1] , tlag_data[1]) ;  cHat = col(x, 2) ./ V2;  cHatDat4 = cHat[iObs];
  x = linOdeModel(time, amt, rate, ii, evid, cmt, addl, ss, K_data[1] , biovar_data    , tlag_data[1]) ;  cHat = col(x, 2) ./ V2;  cHatDat5 = cHat[iObs];
  x = linOdeModel(time, amt, rate, ii, evid, cmt, addl, ss, K_data    , biovar_data[1] , tlag_data)    ;  cHat = col(x, 2) ./ V2;  cHatDat6 = cHat[iObs];
  x = linOdeModel(time, amt, rate, ii, evid, cmt, addl, ss, K_data    , biovar_data[1] , tlag_data[1]) ;  cHat = col(x, 2) ./ V2;  cHatDat7 = cHat[iObs];
  x = linOdeModel(time, amt, rate, ii, evid, cmt, addl, ss, K_data    , biovar_data    , tlag_data)    ;  cHat = col(x, 2) ./ V2;  cHatDat8 = cHat[iObs];

  /* param args */
  x = linOdeModel(time, amt, rate, ii, evid, cmt, addl, ss, K         , biovar         , tlag)         ;  cHat = col(x, 2) ./ V2;  cHatObs1 = cHat[iObs];
  x = linOdeModel(time, amt, rate, ii, evid, cmt, addl, ss, K[1]      , biovar         , tlag)         ;  cHat = col(x, 2) ./ V2;  cHatObs2 = cHat[iObs];  
  x = linOdeModel(time, amt, rate, ii, evid, cmt, addl, ss, K[1]      , biovar[1]      , tlag)         ;  cHat = col(x, 2) ./ V2;  cHatObs3 = cHat[iObs];
  x = linOdeModel(time, amt, rate, ii, evid, cmt, addl, ss, K[1]      , biovar[1]      , tlag[1])      ;  cHat = col(x, 2) ./ V2;  cHatObs4 = cHat[iObs];
  x = linOdeModel(time, amt, rate, ii, evid, cmt, addl, ss, K[1]      , biovar         , tlag[1])      ;  cHat = col(x, 2) ./ V2;  cHatObs5 = cHat[iObs];
  x = linOdeModel(time, amt, rate, ii, evid, cmt, addl, ss, K         , biovar[1]      , tlag)         ;  cHat = col(x, 2) ./ V2;  cHatObs6 = cHat[iObs];
  x = linOdeModel(time, amt, rate, ii, evid, cmt, addl, ss, K         , biovar[1]      , tlag[1])      ;  cHat = col(x, 2) ./ V2;  cHatObs7 = cHat[iObs];
  x = linOdeModel(time, amt, rate, ii, evid, cmt, addl, ss, K         , biovar         , tlag[1])      ;  cHat = col(x, 2) ./ V2;  cHatObs8 = cHat[iObs];

  /* K: data                                                            , biovar: param  , tlag: param */
  x = linOdeModel(time, amt, rate, ii, evid, cmt, addl, ss, K_data    , biovar         , tlag)         ;  cHat = col(x, 2) ./ V2;  cHatDPP1 = cHat[iObs];
  x = linOdeModel(time, amt, rate, ii, evid, cmt, addl, ss, K_data[1] , biovar         , tlag)         ;  cHat = col(x, 2) ./ V2;  cHatDPP2 = cHat[iObs];  
  x = linOdeModel(time, amt, rate, ii, evid, cmt, addl, ss, K_data[1] , biovar[1]      , tlag)         ;  cHat = col(x, 2) ./ V2;  cHatDPP3 = cHat[iObs];
  x = linOdeModel(time, amt, rate, ii, evid, cmt, addl, ss, K_data[1] , biovar[1]      , tlag[1])      ;  cHat = col(x, 2) ./ V2;  cHatDPP4 = cHat[iObs];
  x = linOdeModel(time, amt, rate, ii, evid, cmt, addl, ss, K_data[1] , biovar         , tlag[1])      ;  cHat = col(x, 2) ./ V2;  cHatDPP5 = cHat[iObs];
  x = linOdeModel(time, amt, rate, ii, evid, cmt, addl, ss, K_data    , biovar[1]      , tlag)         ;  cHat = col(x, 2) ./ V2;  cHatDPP6 = cHat[iObs];
  x = linOdeModel(time, amt, rate, ii, evid, cmt, addl, ss, K_data    , biovar[1]      , tlag[1])      ;  cHat = col(x, 2) ./ V2;  cHatDPP7 = cHat[iObs];
  x = linOdeModel(time, amt, rate, ii, evid, cmt, addl, ss, K_data    , biovar         , tlag[1])      ;  cHat = col(x, 2) ./ V2;  cHatDPP8 = cHat[iObs];

  /* K: parm                                                            , biovar: data   , tlag: param */
  x = linOdeModel(time, amt, rate, ii, evid, cmt, addl, ss, K         , biovar_data    , tlag)         ;  cHat = col(x, 2) ./ V2;  cHatPDP1 = cHat[iObs];
  x = linOdeModel(time, amt, rate, ii, evid, cmt, addl, ss, K[1]      , biovar_data    , tlag)         ;  cHat = col(x, 2) ./ V2;  cHatPDP2 = cHat[iObs];  
  x = linOdeModel(time, amt, rate, ii, evid, cmt, addl, ss, K[1]      , biovar_data[1] , tlag)         ;  cHat = col(x, 2) ./ V2;  cHatPDP3 = cHat[iObs];
  x = linOdeModel(time, amt, rate, ii, evid, cmt, addl, ss, K[1]      , biovar_data[1] , tlag[1])      ;  cHat = col(x, 2) ./ V2;  cHatPDP4 = cHat[iObs];
  x = linOdeModel(time, amt, rate, ii, evid, cmt, addl, ss, K[1]      , biovar_data    , tlag[1])      ;  cHat = col(x, 2) ./ V2;  cHatPDP5 = cHat[iObs];
  x = linOdeModel(time, amt, rate, ii, evid, cmt, addl, ss, K         , biovar_data[1] , tlag)         ;  cHat = col(x, 2) ./ V2;  cHatPDP6 = cHat[iObs];
  x = linOdeModel(time, amt, rate, ii, evid, cmt, addl, ss, K         , biovar_data[1] , tlag[1])      ;  cHat = col(x, 2) ./ V2;  cHatPDP7 = cHat[iObs];
  x = linOdeModel(time, amt, rate, ii, evid, cmt, addl, ss, K         , biovar_data    , tlag[1])      ;  cHat = col(x, 2) ./ V2;  cHatPDP8 = cHat[iObs];

  /* K: parm                                                            , biovar: parm   , tlag: data */
  x = linOdeModel(time, amt, rate, ii, evid, cmt, addl, ss, K         , biovar         , tlag_data)    ;  cHat = col(x, 2) ./ V2;  cHatPPD1 = cHat[iObs];
  x = linOdeModel(time, amt, rate, ii, evid, cmt, addl, ss, K[1]      , biovar         , tlag_data)    ;  cHat = col(x, 2) ./ V2;  cHatPPD2 = cHat[iObs];  
  x = linOdeModel(time, amt, rate, ii, evid, cmt, addl, ss, K[1]      , biovar[1]      , tlag_data)    ;  cHat = col(x, 2) ./ V2;  cHatPPD3 = cHat[iObs];
  x = linOdeModel(time, amt, rate, ii, evid, cmt, addl, ss, K[1]      , biovar[1]      , tlag_data[1]) ;  cHat = col(x, 2) ./ V2;  cHatPPD4 = cHat[iObs];
  x = linOdeModel(time, amt, rate, ii, evid, cmt, addl, ss, K[1]      , biovar         , tlag_data[1]) ;  cHat = col(x, 2) ./ V2;  cHatPPD5 = cHat[iObs];
  x = linOdeModel(time, amt, rate, ii, evid, cmt, addl, ss, K         , biovar[1]      , tlag_data)    ;  cHat = col(x, 2) ./ V2;  cHatPPD6 = cHat[iObs];
  x = linOdeModel(time, amt, rate, ii, evid, cmt, addl, ss, K         , biovar[1]      , tlag_data[1]) ;  cHat = col(x, 2) ./ V2;  cHatPPD7 = cHat[iObs];
  x = linOdeModel(time, amt, rate, ii, evid, cmt, addl, ss, K         , biovar         , tlag_data[1]) ;  cHat = col(x, 2) ./ V2;  cHatPPD8 = cHat[iObs];

  /* K: data                                                            , biovar: data   , tlag: param */
  x = linOdeModel(time, amt, rate, ii, evid, cmt, addl, ss, K_data    , biovar_data    , tlag)         ;  cHat = col(x, 2) ./ V2;  cHatDDP1 = cHat[iObs];
  x = linOdeModel(time, amt, rate, ii, evid, cmt, addl, ss, K_data[1] , biovar_data    , tlag)         ;  cHat = col(x, 2) ./ V2;  cHatDDP2 = cHat[iObs];  
  x = linOdeModel(time, amt, rate, ii, evid, cmt, addl, ss, K_data[1] , biovar_data[1] , tlag)         ;  cHat = col(x, 2) ./ V2;  cHatDDP3 = cHat[iObs];
  x = linOdeModel(time, amt, rate, ii, evid, cmt, addl, ss, K_data[1] , biovar_data[1] , tlag[1])      ;  cHat = col(x, 2) ./ V2;  cHatDDP4 = cHat[iObs];
  x = linOdeModel(time, amt, rate, ii, evid, cmt, addl, ss, K_data[1] , biovar_data    , tlag[1])      ;  cHat = col(x, 2) ./ V2;  cHatDDP5 = cHat[iObs];
  x = linOdeModel(time, amt, rate, ii, evid, cmt, addl, ss, K_data    , biovar_data[1] , tlag)         ;  cHat = col(x, 2) ./ V2;  cHatDDP6 = cHat[iObs];
  x = linOdeModel(time, amt, rate, ii, evid, cmt, addl, ss, K_data    , biovar_data[1] , tlag[1])      ;  cHat = col(x, 2) ./ V2;  cHatDDP7 = cHat[iObs];
  x = linOdeModel(time, amt, rate, ii, evid, cmt, addl, ss, K_data    , biovar_data    , tlag[1])      ;  cHat = col(x, 2) ./ V2;  cHatDDP8 = cHat[iObs];

  /* K: data                                                            , biovar: param  , tlag: data */
  x = linOdeModel(time, amt, rate, ii, evid, cmt, addl, ss, K_data    , biovar         , tlag_data)    ;  cHat = col(x, 2) ./ V2;  cHatDPD1 = cHat[iObs];
  x = linOdeModel(time, amt, rate, ii, evid, cmt, addl, ss, K_data[1] , biovar         , tlag_data)    ;  cHat = col(x, 2) ./ V2;  cHatDPD2 = cHat[iObs];  
  x = linOdeModel(time, amt, rate, ii, evid, cmt, addl, ss, K_data[1] , biovar[1]      , tlag_data)    ;  cHat = col(x, 2) ./ V2;  cHatDPD3 = cHat[iObs];
  x = linOdeModel(time, amt, rate, ii, evid, cmt, addl, ss, K_data[1] , biovar[1]      , tlag_data[1]) ;  cHat = col(x, 2) ./ V2;  cHatDPD4 = cHat[iObs];
  x = linOdeModel(time, amt, rate, ii, evid, cmt, addl, ss, K_data[1] , biovar         , tlag_data[1]) ;  cHat = col(x, 2) ./ V2;  cHatDPD5 = cHat[iObs];
  x = linOdeModel(time, amt, rate, ii, evid, cmt, addl, ss, K_data    , biovar[1]      , tlag_data)    ;  cHat = col(x, 2) ./ V2;  cHatDPD6 = cHat[iObs];
  x = linOdeModel(time, amt, rate, ii, evid, cmt, addl, ss, K_data    , biovar[1]      , tlag_data[1]) ;  cHat = col(x, 2) ./ V2;  cHatDPD7 = cHat[iObs];
  x = linOdeModel(time, amt, rate, ii, evid, cmt, addl, ss, K_data    , biovar         , tlag_data[1]) ;  cHat = col(x, 2) ./ V2;  cHatDPD8 = cHat[iObs];

  /* K: param                                                           , biovar: data   , tlag: data */
  x = linOdeModel(time, amt, rate, ii, evid, cmt, addl, ss, K         , biovar_data    , tlag_data)    ;  cHat = col(x, 2) ./ V2;  cHatPDD1 = cHat[iObs];
  x = linOdeModel(time, amt, rate, ii, evid, cmt, addl, ss, K[1]      , biovar_data    , tlag_data)    ;  cHat = col(x, 2) ./ V2;  cHatPDD2 = cHat[iObs];  
  x = linOdeModel(time, amt, rate, ii, evid, cmt, addl, ss, K[1]      , biovar_data[1] , tlag_data)    ;  cHat = col(x, 2) ./ V2;  cHatPDD3 = cHat[iObs];
  x = linOdeModel(time, amt, rate, ii, evid, cmt, addl, ss, K[1]      , biovar_data[1] , tlag_data[1]) ;  cHat = col(x, 2) ./ V2;  cHatPDD4 = cHat[iObs];
  x = linOdeModel(time, amt, rate, ii, evid, cmt, addl, ss, K[1]      , biovar_data    , tlag_data[1]) ;  cHat = col(x, 2) ./ V2;  cHatPDD5 = cHat[iObs];
  x = linOdeModel(time, amt, rate, ii, evid, cmt, addl, ss, K         , biovar_data[1] , tlag_data)    ;  cHat = col(x, 2) ./ V2;  cHatPDD6 = cHat[iObs];
  x = linOdeModel(time, amt, rate, ii, evid, cmt, addl, ss, K         , biovar_data[1] , tlag_data[1]) ;  cHat = col(x, 2) ./ V2;  cHatPDD7 = cHat[iObs];
  x = linOdeModel(time, amt, rate, ii, evid, cmt, addl, ss, K         , biovar_data    , tlag_data[1]) ;  cHat = col(x, 2) ./ V2;  cHatPDD8 = cHat[iObs];

}

