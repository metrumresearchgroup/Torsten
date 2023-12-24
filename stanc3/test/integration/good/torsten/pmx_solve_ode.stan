functions {
  vector ode_rhs(real t, vector x, array[] real parms, array[] real x_r, array[] int x_i){
    real CL = parms[1];
    real Q = parms[2];
    real V1 = parms[3];
    real V2 = parms[4];
    real ka = parms[5];
    
    real k10 = CL / V1;
    real k12 = Q / V1;
    real k21 = Q / V2;
    
    vector[3] y;

    y[1] = -ka*x[1];
    y[2] = ka*x[1] - (k10 + k12)*x[2] + k21*x[3];
    y[3] = k12*x[2] - k21*x[3];

    return y;
  }  
}
data{
  int<lower = 1> nt;  // number of events
  int<lower = 1> nObs;  // number of observation
  array[nObs] int<lower = 1> iObs;  // index of observation
  
  // NONMEM data
  array[nt] int<lower = 1> cmt;
  array[nt] int evid;
  array[nt] int addl;
  array[nt] int ss;
  array[nt] real amt;
  array[nt] real time;
  array[nt] real rate;
  array[nt] real ii;

  real rel_tol;
  real abs_tol;
  int max_num_steps;
  real as_rel_tol;
  real as_abs_tol;
  int as_max_num_steps;
}

parameters{
  real<lower = 0> CL;
  real<lower = 0> Q;
  real<lower = 0> V1;
  real<lower = 0> V2;
  real<lower = 0> ka;
  real<lower = 0> ke;
  real<lower = 0> sigma;
}

transformed parameters{
  array[5] real theta;  // ODE parameters
  array[3] real biovar;
  array[3] real tlag;
  array[nt, 5] real theta_t;
  array[nt, 3] real biovar_t;
  array[nt, 3] real tlag_t;
  row_vector<lower = 0>[nt] cHat;
  matrix<lower = 0>[3, nt] x;

  theta[1] = CL;
  theta[2] = Q;
  theta[3] = V1;
  theta[4] = V2;
  theta[5] = ka;

  x = pmx_solve_rk45(ode_rhs, 3, time, amt, rate, ii, evid, cmt, addl, ss, theta, rel_tol, abs_tol, max_num_steps, as_rel_tol, as_abs_tol, as_max_num_steps);
  x = pmx_solve_rk45(ode_rhs, 3, time, amt, rate, ii, evid, cmt, addl, ss, theta, biovar, rel_tol, abs_tol, max_num_steps, as_rel_tol, as_abs_tol, as_max_num_steps);
  x = pmx_solve_rk45(ode_rhs, 3, time, amt, rate, ii, evid, cmt, addl, ss, theta, biovar, tlag, rel_tol, abs_tol, max_num_steps, as_rel_tol, as_abs_tol, as_max_num_steps);
  x = pmx_solve_rk45(ode_rhs, 3, time, amt, rate, ii, evid, cmt, addl, ss, theta, rel_tol, abs_tol, max_num_steps);
  x = pmx_solve_rk45(ode_rhs, 3, time, amt, rate, ii, evid, cmt, addl, ss, theta, biovar, rel_tol, abs_tol, max_num_steps);
  x = pmx_solve_rk45(ode_rhs, 3, time, amt, rate, ii, evid, cmt, addl, ss, theta, biovar, tlag, rel_tol, abs_tol, max_num_steps);

  x = pmx_solve_rk45(ode_rhs, 3, time, amt, rate, ii, evid, cmt, addl, ss, theta_t);
  x = pmx_solve_rk45(ode_rhs, 3, time, amt, rate, ii, evid, cmt, addl, ss, theta_t, biovar_t);
  x = pmx_solve_rk45(ode_rhs, 3, time, amt, rate, ii, evid, cmt, addl, ss, theta_t, biovar_t, tlag_t);
  x = pmx_solve_rk45(ode_rhs, 3, time, amt, rate, ii, evid, cmt, addl, ss, theta_t, biovar);
  x = pmx_solve_rk45(ode_rhs, 3, time, amt, rate, ii, evid, cmt, addl, ss, theta_t, biovar, tlag_t);
  x = pmx_solve_rk45(ode_rhs, 3, time, amt, rate, ii, evid, cmt, addl, ss, theta, biovar_t);
  x = pmx_solve_rk45(ode_rhs, 3, time, amt, rate, ii, evid, cmt, addl, ss, theta, biovar_t, tlag_t);
  x = pmx_solve_rk45(ode_rhs, 3, time, amt, rate, ii, evid, cmt, addl, ss, theta_t, biovar_t);
  x = pmx_solve_rk45(ode_rhs, 3, time, amt, rate, ii, evid, cmt, addl, ss, theta_t, biovar_t, tlag);

  x = pmx_solve_bdf(ode_rhs, 3, time, amt, rate, ii, evid, cmt, addl, ss, theta, rel_tol, abs_tol, max_num_steps, as_rel_tol, as_abs_tol, as_max_num_steps);
  x = pmx_solve_bdf(ode_rhs, 3, time, amt, rate, ii, evid, cmt, addl, ss, theta, biovar, rel_tol, abs_tol, max_num_steps, as_rel_tol, as_abs_tol, as_max_num_steps);
  x = pmx_solve_bdf(ode_rhs, 3, time, amt, rate, ii, evid, cmt, addl, ss, theta, biovar, tlag, rel_tol, abs_tol, max_num_steps, as_rel_tol, as_abs_tol, as_max_num_steps);
  x = pmx_solve_bdf(ode_rhs, 3, time, amt, rate, ii, evid, cmt, addl, ss, theta, rel_tol, abs_tol, max_num_steps);
  x = pmx_solve_bdf(ode_rhs, 3, time, amt, rate, ii, evid, cmt, addl, ss, theta, biovar, rel_tol, abs_tol, max_num_steps);
  x = pmx_solve_bdf(ode_rhs, 3, time, amt, rate, ii, evid, cmt, addl, ss, theta, biovar, tlag, rel_tol, abs_tol, max_num_steps);

  x = pmx_solve_bdf(ode_rhs, 3, time, amt, rate, ii, evid, cmt, addl, ss, theta_t);
  x = pmx_solve_bdf(ode_rhs, 3, time, amt, rate, ii, evid, cmt, addl, ss, theta_t, biovar_t);
  x = pmx_solve_bdf(ode_rhs, 3, time, amt, rate, ii, evid, cmt, addl, ss, theta_t, biovar_t, tlag_t);
  x = pmx_solve_bdf(ode_rhs, 3, time, amt, rate, ii, evid, cmt, addl, ss, theta_t, biovar);
  x = pmx_solve_bdf(ode_rhs, 3, time, amt, rate, ii, evid, cmt, addl, ss, theta_t, biovar, tlag_t);
  x = pmx_solve_bdf(ode_rhs, 3, time, amt, rate, ii, evid, cmt, addl, ss, theta, biovar_t);
  x = pmx_solve_bdf(ode_rhs, 3, time, amt, rate, ii, evid, cmt, addl, ss, theta, biovar_t, tlag_t);
  x = pmx_solve_bdf(ode_rhs, 3, time, amt, rate, ii, evid, cmt, addl, ss, theta_t, biovar_t);
  x = pmx_solve_bdf(ode_rhs, 3, time, amt, rate, ii, evid, cmt, addl, ss, theta_t, biovar_t, tlag);

  x = pmx_solve_adams(ode_rhs, 3, time, amt, rate, ii, evid, cmt, addl, ss, theta, rel_tol, abs_tol, max_num_steps, as_rel_tol, as_abs_tol, as_max_num_steps);
  x = pmx_solve_adams(ode_rhs, 3, time, amt, rate, ii, evid, cmt, addl, ss, theta, biovar, rel_tol, abs_tol, max_num_steps, as_rel_tol, as_abs_tol, as_max_num_steps);
  x = pmx_solve_adams(ode_rhs, 3, time, amt, rate, ii, evid, cmt, addl, ss, theta, biovar, tlag, rel_tol, abs_tol, max_num_steps, as_rel_tol, as_abs_tol, as_max_num_steps);
  x = pmx_solve_adams(ode_rhs, 3, time, amt, rate, ii, evid, cmt, addl, ss, theta, rel_tol, abs_tol, max_num_steps);
  x = pmx_solve_adams(ode_rhs, 3, time, amt, rate, ii, evid, cmt, addl, ss, theta, biovar, rel_tol, abs_tol, max_num_steps);
  x = pmx_solve_adams(ode_rhs, 3, time, amt, rate, ii, evid, cmt, addl, ss, theta, biovar, tlag, rel_tol, abs_tol, max_num_steps);

  x = pmx_solve_adams(ode_rhs, 3, time, amt, rate, ii, evid, cmt, addl, ss, theta_t);
  x = pmx_solve_adams(ode_rhs, 3, time, amt, rate, ii, evid, cmt, addl, ss, theta_t, biovar_t);
  x = pmx_solve_adams(ode_rhs, 3, time, amt, rate, ii, evid, cmt, addl, ss, theta_t, biovar_t, tlag_t);
  x = pmx_solve_adams(ode_rhs, 3, time, amt, rate, ii, evid, cmt, addl, ss, theta_t, biovar);
  x = pmx_solve_adams(ode_rhs, 3, time, amt, rate, ii, evid, cmt, addl, ss, theta_t, biovar, tlag_t);
  x = pmx_solve_adams(ode_rhs, 3, time, amt, rate, ii, evid, cmt, addl, ss, theta, biovar_t);
  x = pmx_solve_adams(ode_rhs, 3, time, amt, rate, ii, evid, cmt, addl, ss, theta, biovar_t, tlag_t);
  x = pmx_solve_adams(ode_rhs, 3, time, amt, rate, ii, evid, cmt, addl, ss, theta_t, biovar_t);
  x = pmx_solve_adams(ode_rhs, 3, time, amt, rate, ii, evid, cmt, addl, ss, theta_t, biovar_t, tlag);

  cHat = x[2, :] ./ V1; // we're interested in the amount in the second compartment
}

model{
}
