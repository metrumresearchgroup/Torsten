functions {
  real[] ode(real t,
             real[] y,
             real[] theta,
             real[] x,
             int[] x_int) {
    real dydt[2];
    return dydt;
  }
}

data {
  int<lower = 1> np;
  int<lower = 1> nt;
  int nTheta;
  int<lower = 1> cmt[nt];
  int len[np];
  int evid[nt];
  int addl[nt];
  int ss[nt];
  real amt[nt];
  real time[nt];
  real rate[nt];
  real ii[nt];

int<lower=1> T;
real y0_d[2];
real t0;
real ts[T];
real theta_d[1];
real x[0];
int x_int[0];
}

transformed data {
  int nCmt = 2;
  real theta_data[nt, nTheta];
  real biovar_data[nt, nCmt];
  real tlag_data[nt, nCmt];
  matrix[nCmt, nt * np] x_data;
}

parameters {
  real y_p;
  real y0_p[2];
  real theta_p[1];
}

transformed parameters {
  real theta_parm[nt, nCmt];
  real biovar_parm;
  real tlag_parm[nt, nCmt];
  matrix[nCmt, nt * np] x_parm;

  /*****************************************************************
   pmx_solve_group_ode
   ****************************************************************/
  // RK45
  x_parm = pmx_solve_group_rk45(ode, nCmt, len, time, amt, rate, ii, evid, cmt, addl, ss, theta_parm, biovar_parm, tlag_parm);
}

model {
	y_p ~ normal(0,1);
}
