functions {
  vector ode(real t,
             vector y,
             array[] real theta,
             array[] real x,
             array[] int x_int) {
    vector[2] dydt;
    return dydt;
  }
}

data {
  int<lower = 1> np;
  int<lower = 1> nt;
  int nTheta;
  array[nt] int<lower = 1> cmt;
  array[np] int len;
  array[nt] int evid;
  array[nt] int addl;
  array[nt] int ss;
  array[nt] real amt;
  array[nt] real time;
  array[nt] real rate;
  array[nt] real ii;

int<lower=1> T;
// real y0_d[2];
// real t0;
// real ts[T];
array[1] real theta_d;
// real x[0];
// int x_int[0];
}

transformed data {
  int nCmt = 2;
  array[nt, nTheta] real theta_data;
  array[nt, nCmt] real biovar_data;
  array[nt, nCmt] real tlag_data;
  matrix[nCmt, nt * np] x_data;

  array[nt, 2] real x_r;
  array[nt, 3] int x_i;

  /*****************************************************************
   pmx_solve_rk45/adams/bdf full sig
   *****************************************************************/
  x_data = pmx_solve_group_rk45  (ode, nCmt, len, time, amt, rate, ii, evid, cmt, addl, ss, theta_data    , biovar_data    , tlag_data    , 1e-8, 1e-8, 1e8);
  x_data = pmx_solve_group_bdf   (ode, nCmt, len, time, amt, rate, ii, evid, cmt, addl, ss, theta_data    , biovar_data    , tlag_data    , 1e-8, 1e-8, 1e8);
  x_data = pmx_solve_group_adams (ode, nCmt, len, time, amt, rate, ii, evid, cmt, addl, ss, theta_data    , biovar_data    , tlag_data    , 1e-8, 1e-8, 1e8);

  /*****************************************************************
   pmx_solve_rk45/adams/bdf default tlag
   *****************************************************************/

  x_data = pmx_solve_group_rk45  (ode, nCmt, len, time, amt, rate, ii, evid, cmt, addl, ss, theta_data    , biovar_data    , 1e-8, 1e-8, 1e8);
  x_data = pmx_solve_group_bdf   (ode, nCmt, len, time, amt, rate, ii, evid, cmt, addl, ss, theta_data    , biovar_data    , 1e-8, 1e-8, 1e8);
  x_data = pmx_solve_group_adams (ode, nCmt, len, time, amt, rate, ii, evid, cmt, addl, ss, theta_data    , biovar_data    , 1e-8, 1e-8, 1e8);

  /*****************************************************************
   pmx_solve_rk45/adams/bdf default F & tlag
   *****************************************************************/

  x_data = pmx_solve_group_rk45  (ode, nCmt, len, time, amt, rate, ii, evid, cmt, addl, ss, theta_data    , 1e-8, 1e-8, 1e8);
  x_data = pmx_solve_group_bdf   (ode, nCmt, len, time, amt, rate, ii, evid, cmt, addl, ss, theta_data    , 1e-8, 1e-8, 1e8);
  x_data = pmx_solve_group_adams (ode, nCmt, len, time, amt, rate, ii, evid, cmt, addl, ss, theta_data    , 1e-8, 1e-8, 1e8);

  /*****************************************************************
   pmx_solve_group_rk45/adams/bdf w/o ODE controls
   *****************************************************************/
  x_data = pmx_solve_group_rk45  (ode, nCmt, len, time, amt, rate, ii, evid, cmt, addl, ss, theta_data    , biovar_data    , tlag_data    );
  x_data = pmx_solve_group_bdf   (ode, nCmt, len, time, amt, rate, ii, evid, cmt, addl, ss, theta_data    , biovar_data    , tlag_data    );
  x_data = pmx_solve_group_adams (ode, nCmt, len, time, amt, rate, ii, evid, cmt, addl, ss, theta_data    , biovar_data    , tlag_data    );

  /*****************************************************************
   pmx_solve_group_rk45/adams/bdf w/o ODE controls or tlag
   *****************************************************************/
  x_data = pmx_solve_group_rk45  (ode, nCmt, len, time, amt, rate, ii, evid, cmt, addl, ss, theta_data    , biovar_data);
  x_data = pmx_solve_group_bdf   (ode, nCmt, len, time, amt, rate, ii, evid, cmt, addl, ss, theta_data    , biovar_data);
  x_data = pmx_solve_group_adams (ode, nCmt, len, time, amt, rate, ii, evid, cmt, addl, ss, theta_data    , biovar_data);

  /*****************************************************************
   pmx_solve_group_rk45/adams/bdf w/o ODE controls or tlag
   *****************************************************************/
  x_data = pmx_solve_group_rk45  (ode, nCmt, len, time, amt, rate, ii, evid, cmt, addl, ss, theta_data);
  x_data = pmx_solve_group_bdf   (ode, nCmt, len, time, amt, rate, ii, evid, cmt, addl, ss, theta_data);
  x_data = pmx_solve_group_adams (ode, nCmt, len, time, amt, rate, ii, evid, cmt, addl, ss, theta_data);
}

parameters {
  real y_p;

  array[2] real y0_p;
  array[1] real theta_p;
}

transformed parameters {
  array[nt, nTheta] real theta_parm;
  array[nt, nCmt] real biovar_parm;
  array[nt, nCmt] real tlag_parm;
  matrix[nCmt, nt * np] x_parm;

  /*****************************************************************
   pmx_solve_group_ode
   ****************************************************************/
  x_parm = pmx_solve_group_rk45  (ode, nCmt, len, time, amt, rate, ii, evid, cmt, addl, ss, theta_parm    , biovar_parm    , tlag_parm    , 1e-8, 1e-8, 1e8);
  x_parm = pmx_solve_group_bdf   (ode, nCmt, len, time, amt, rate, ii, evid, cmt, addl, ss, theta_parm    , biovar_parm    , tlag_parm    , 1e-8, 1e-8, 1e8);
  x_parm = pmx_solve_group_adams (ode, nCmt, len, time, amt, rate, ii, evid, cmt, addl, ss, theta_parm    , biovar_parm    , tlag_parm    , 1e-8, 1e-8, 1e8);
  x_parm = pmx_solve_group_rk45  (ode, nCmt, len, time, amt, rate, ii, evid, cmt, addl, ss, theta_parm    , biovar_parm    , 1e-8, 1e-8, 1e8);
  x_parm = pmx_solve_group_bdf   (ode, nCmt, len, time, amt, rate, ii, evid, cmt, addl, ss, theta_parm    , biovar_parm    , 1e-8, 1e-8, 1e8);
  x_parm = pmx_solve_group_adams (ode, nCmt, len, time, amt, rate, ii, evid, cmt, addl, ss, theta_parm    , biovar_parm    , 1e-8, 1e-8, 1e8);
  x_parm = pmx_solve_group_rk45  (ode, nCmt, len, time, amt, rate, ii, evid, cmt, addl, ss, theta_parm    , 1e-8, 1e-8, 1e8);
  x_parm = pmx_solve_group_bdf   (ode, nCmt, len, time, amt, rate, ii, evid, cmt, addl, ss, theta_parm    , 1e-8, 1e-8, 1e8);
  x_parm = pmx_solve_group_adams (ode, nCmt, len, time, amt, rate, ii, evid, cmt, addl, ss, theta_parm    , 1e-8, 1e-8, 1e8);


  /*****************************************************************
   pmx_solve_group_ode no ODE controls
   ****************************************************************/
  x_parm = pmx_solve_group_rk45  (ode, nCmt, len, time, amt, rate, ii, evid, cmt, addl, ss, theta_parm    , biovar_parm    , tlag_parm);
  x_parm = pmx_solve_group_bdf   (ode, nCmt, len, time, amt, rate, ii, evid, cmt, addl, ss, theta_parm    , biovar_parm    , tlag_parm);
  x_parm = pmx_solve_group_adams (ode, nCmt, len, time, amt, rate, ii, evid, cmt, addl, ss, theta_parm    , biovar_parm    , tlag_parm);
  x_parm = pmx_solve_group_rk45  (ode, nCmt, len, time, amt, rate, ii, evid, cmt, addl, ss, theta_parm    , biovar_parm);
  x_parm = pmx_solve_group_bdf   (ode, nCmt, len, time, amt, rate, ii, evid, cmt, addl, ss, theta_parm    , biovar_parm);
  x_parm = pmx_solve_group_adams (ode, nCmt, len, time, amt, rate, ii, evid, cmt, addl, ss, theta_parm    , biovar_parm);
  x_parm = pmx_solve_group_rk45  (ode, nCmt, len, time, amt, rate, ii, evid, cmt, addl, ss, theta_parm);
  x_parm = pmx_solve_group_bdf   (ode, nCmt, len, time, amt, rate, ii, evid, cmt, addl, ss, theta_parm);
  x_parm = pmx_solve_group_adams (ode, nCmt, len, time, amt, rate, ii, evid, cmt, addl, ss, theta_parm);

  /*****************************************************************
   pmx_solve_group_ode with data & w/o ODE controls
   ****************************************************************/
  x_parm = pmx_solve_group_rk45  (ode, nCmt, len, time, amt, rate, ii, evid, cmt, addl, ss, theta_parm    , biovar_parm    , tlag_parm, x_r, x_i);
  x_parm = pmx_solve_group_rk45  (ode, nCmt, len, time, amt, rate, ii, evid, cmt, addl, ss, theta_parm    , biovar_parm    , tlag_parm, x_r);
  x_parm = pmx_solve_group_bdf   (ode, nCmt, len, time, amt, rate, ii, evid, cmt, addl, ss, theta_parm    , biovar_parm    , tlag_parm, x_r, x_i);
  x_parm = pmx_solve_group_bdf   (ode, nCmt, len, time, amt, rate, ii, evid, cmt, addl, ss, theta_parm    , biovar_parm    , tlag_parm, x_r);
  x_parm = pmx_solve_group_adams (ode, nCmt, len, time, amt, rate, ii, evid, cmt, addl, ss, theta_parm    , biovar_parm    , tlag_parm, x_r, x_i);
  x_parm = pmx_solve_group_adams (ode, nCmt, len, time, amt, rate, ii, evid, cmt, addl, ss, theta_parm    , biovar_parm    , tlag_parm, x_r);

  x_parm = pmx_solve_group_rk45  (ode, nCmt, len, time, amt, rate, ii, evid, cmt, addl, ss, theta_parm    , biovar_parm    , tlag_parm, x_r, x_i, 1e-8, 1e-8, 1e8);
  x_parm = pmx_solve_group_rk45  (ode, nCmt, len, time, amt, rate, ii, evid, cmt, addl, ss, theta_parm    , biovar_parm    , tlag_parm, x_r, 1e-8, 1e-8, 1e8);
  x_parm = pmx_solve_group_bdf   (ode, nCmt, len, time, amt, rate, ii, evid, cmt, addl, ss, theta_parm    , biovar_parm    , tlag_parm, x_r, x_i, 1e-8, 1e-8, 1e8);
  x_parm = pmx_solve_group_bdf   (ode, nCmt, len, time, amt, rate, ii, evid, cmt, addl, ss, theta_parm    , biovar_parm    , tlag_parm, x_r, 1e-8, 1e-8, 1e8);
  x_parm = pmx_solve_group_adams (ode, nCmt, len, time, amt, rate, ii, evid, cmt, addl, ss, theta_parm    , biovar_parm    , tlag_parm, x_r, x_i, 1e-8, 1e-8, 1e8);
  x_parm = pmx_solve_group_adams (ode, nCmt, len, time, amt, rate, ii, evid, cmt, addl, ss, theta_parm    , biovar_parm    , tlag_parm, x_r, 1e-8, 1e-8, 1e8);
}

model {
	y_p ~ normal(0,1);
}
