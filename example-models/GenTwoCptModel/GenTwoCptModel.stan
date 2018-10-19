// GenTwoCptModelExample.stan
// Run two compartment model using numerical solution
// Heavily anotated to help new users

functions{
  
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
  int<lower = 1> nObs;  // number of observations
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
  
  vector<lower = 0>[nObs] cObs;  // observed concentration (dependent variable)
}

transformed data{
  vector[nObs] logCObs = log(cObs);
  int nTheta = 5;   // number of parameters
  int nCmt = 3;   // number of compartments
  real biovar[nCmt];
  real tlag[nCmt];

  for (i in 1:nCmt) {
    biovar[i] = 1;
    tlag[i] = 0;
  }
}

parameters{
  real<lower = 0> CL;
  real<lower = 0> Q;
  real<lower = 0> V1;
  real<lower = 0> V2;
  real<lower = 0> ka;
  real<lower = 0> sigma;
}

transformed parameters{
  real theta[nTheta];
  vector<lower = 0>[nt] cHat;
  vector<lower = 0>[nObs] cHatObs;
  matrix<lower = 0>[nt, 3] x; 

  theta[1] = CL;
  theta[2] = Q;
  theta[3] = V1;
  theta[4] = V2;
  theta[5] = ka;

  // generalCptModel takes in the ODE system, the number of compartment 
  // (here we have a two compartment model with first order absorption, so
  // three compartments), the parameters matrix, the NONEM data, and the tuning
  // parameters (relative tolerance, absolute tolerance, and maximum number of steps)
  // of the ODE integrator. The user can choose between the bdf and the rk45 integrator.
  // Returns a matrix with the predicted amount in each compartment 
  // at each event.

//  x = generalOdeModel_bdf(twoCptModelODE, 3,
//                          time, amt, rate, ii, evid, cmt, addl, ss,
//                          theta, biovar, tlag,
//                          1e-8, 1e-8, 1e8);

   x = generalOdeModel_rk45(twoCptModelODE, 3,
                           time, amt, rate, ii, evid, cmt, addl, ss,
                           theta, biovar, tlag,
                           1e-6, 1e-6, 1e6);

  cHat = col(x, 2) ./ V1;

  for(i in 1:nObs){
    cHatObs[i] = cHat[iObs[i]];  // predictions for observed data records
  }
}

model{
  // informative prior
  CL ~ lognormal(log(10), 0.25);
  Q ~ lognormal(log(15), 0.5);
  V1 ~ lognormal(log(35), 0.25);
  V2 ~ lognormal(log(105), 0.5);
  ka ~ lognormal(log(2.5), 1);
  sigma ~ cauchy(0, 1);

  logCObs ~ normal(log(cHatObs), sigma);
}

generated quantities{
  real cObsPred[nObs];

  for(i in 1:nObs){
      cObsPred[i] = exp(normal_rng(log(cHatObs[i]), sigma));
    }
			 
}
