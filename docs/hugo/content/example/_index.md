+++
title = "Examples"
author = ["Yi Zhang"]
date = 2021-06-25T00:00:00-07:00
lastmod = 2021-06-26T09:38:18-07:00
draft = false
weight = 1006
+++

All the PMX models in this chapter can be found in
`Torsten/example-models` directory:

-   `Torsten/example-models/pk2cpt` (Section [sec:pk2cpt](#sec:pk2cpt)).
-   `Torsten/example-models/pk2cpt_linode` (Section [sec:pk2cpt_linode](#sec:pk2cpt_linode)).
-   `Torsten/example-models/pk2cpt_ode` (Section [sec:pk2cpt_ode](#sec:pk2cpt_ode)).
-   `Torsten/example-models/FK_coupled` (Section [sec:fk_model](#sec:fk_model)).
-   `Torsten/example-models/twocpt_population` (Section [sec:twocpt_population](#sec:twocpt_population)).
-   `Torsten/example-models/lotka_volterra_ode_group_model` (Section [sec:lotka_volterra](#sec:lotka_volterra)).
-   `Torsten/example-models/effCpt` (Section [sec:effcpt_model](#sec:effcpt_model)).
-   `Torsten/example-models/FribergKarlsson` (Section [sec:fkpop_model](#sec:fkpop_model)).


## <span class="section-num">1</span> Two-compartment model for single patient {#two-compartment-model-for-single-patient}

\label{sec:pk2cpt}
  We model drug absorption in a single patient and simulate plasma drug concentrations:

-   Multiple Doses: 1250 mg, every 12 hours, for a total of 15 doses
-   PK measured at 0.083, 0.167, 0.25, 0.5, 0.75, 1, 1.5, 2, 4, 6,
    8, 10 and 12 hours after 1st, 2nd, and 15th dose. In addition, the
    PK is measured every 12 hours throughout the trial.

With the plasma concentration \\(\hat{c}\\) solved from
two-compartment ODEs in [{{< relref "two-cpt" >}}]({{< relref "two-cpt" >}}), we simulate \\(c\\) according to:

\begin{align\*}
  \log\left(c\right) &\sim N\left(\log\left(\widehat{c}\right),\sigma^2\right) \\\\\\
  \left(CL, Q, V\_2, V\_3, ka\right) &= \left(5\ {\rm L/h}, 8\  {\rm L/h}, 20\  {\rm L},  70\ {\rm L}, 1.2\ {\rm h^{-1}} \right) \\\\\\
  \sigma^2 &= 0.01
\end{align\*}

The data are generated using the R package `mrgsolve` <sup id="8dd98ac45050f8bfe813328322213083"><a href="#Baron000" title="Kyle Baron \&amp; Marc Gastonguay, Simulation from ODE-Based Population PK/PD and Systems Pharmacology Models in R with mrgsolve, {Journal of Pharmacokinetics and Pharmacodynamics}, v(W-23), S84--S85 (2015).">Baron000</a></sup>.

Code below shows how Torsten function `pmx_solve_twocpt` can be used to fit the above model.

```stan
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

  vector<lower = 0>[nObs] cObs;  // observed concentration (Dependent Variable)
}

transformed data{
  vector[nObs] logCObs = log(cObs);
  int nTheta = 5;  // number of ODE parameters in Two Compartment Model
  int nCmt = 3;  // number of compartments in model
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
  real theta[nTheta];  // ODE parameters
  row_vector<lower = 0>[nt] cHat;
  vector<lower = 0>[nObs] cHatObs;
  matrix<lower = 0>[nCmt, nt] x;

  theta[1] = CL;
  theta[2] = Q;
  theta[3] = V1;
  theta[4] = V2;
  theta[5] = ka;

  x = pmx_solve_twocpt(time, amt, rate, ii, evid, cmt, addl, ss, theta);

  cHat = x[2, :] ./ V1; // we're interested in the amount in the second compartment

  cHatObs = cHat'[iObs]; // predictions for observed data recors
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
```

Four MCMC chains of 2000 iterations (1000 warmup iterations and 1000
sampling iterations) are simulated. 1000 samples per chain were used for the subsequent analyses.
The MCMC history plots(Figure [1](#orgbe0d7f0))
suggest that the 4 chains have converged to common distributions for
all of the key model parameters. The fit to the plasma concentration
data (Figure [3](#org147efd5)) are in close agreement with the
data, which is not surprising since the fitted model is identical to
the one used to simulate the data. Similarly the parameter posterior
density can be examined in Figure [2](#orgffecbb9) and shows
consistency with the values used for simulation. Another way to
summarize the posterior is through `cmdstanr`'s `summary` method.

```r
## fit is a CmdStanMCMC object returned by sampling. See cmdstanr reference.
> pars = c("CL", "Q", "V1", "V2", "ka", "sigma")
> fit$summary(pars)
# A tibble: 6 x 10
  variable   mean median     sd    mad      q5    q95  rhat ess_bulk ess_tail
  <chr>     <dbl>  <dbl>  <dbl>  <dbl>   <dbl>  <dbl> <dbl>    <dbl>    <dbl>
1 CL        4.82   4.83  0.0910 0.0870  4.68    4.97   1.00    1439.    1067.
2 Q         7.56   7.55  0.588  0.586   6.61    8.56   1.00    1256.    1235.
3 V1       21.1   21.1   2.50   2.45   17.1    25.3    1.00    1057.    1177.
4 V2       76.1   76.1   5.33   4.93   67.5    84.9    1.01    1585.    1372.
5 ka        1.23   1.23  0.175  0.174   0.958   1.52   1.00    1070.    1122.
6 sigma     0.109  0.108 0.0117 0.0111  0.0911  0.130  1.01    1414.     905.
```

<a id="orgbe0d7f0"></a>

</ox-hugo/history.pdf>

<a id="orgffecbb9"></a>

</ox-hugo/density.pdf>

<a id="org147efd5"></a>

</ox-hugo/ppc_ribbon.pdf>


## <span class="section-num">2</span> Two-compartment model as a linear ODE model for single patient {#two-compartment-model-as-a-linear-ode-model-for-single-patient}

\label{sec:pk2cpt\_linode}
Using `pmx_solve_linode`, the following example fits a two-compartment model
with first order absorption. We omit `data` and
`model` block as they are identical to Sectiontion [sec:pk2cpt](#sec:pk2cpt).

```stan
transformed data{
  row_vector[nObs] logCObs = log(cObs);
  int nCmt = 3;
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
  matrix[3, 3] K;
  real k10 = CL / V1;
  real k12 = Q / V1;
  real k21 = Q / V2;
  row_vector<lower = 0>[nt] cHat;
  row_vector<lower = 0>[nObs] cHatObs;
  matrix<lower = 0>[3, nt] x;

  K = rep_matrix(0, 3, 3);

  K[1, 1] = -ka;
  K[2, 1] = ka;
  K[2, 2] = -(k10 + k12);
  K[2, 3] = k21;
  K[3, 2] = k12;
  K[3, 3] = -k21;

  x = pmx_solve_linode(time, amt, rate, ii, evid, cmt, addl, ss, K, biovar, tlag);

  cHat = row(x, 2) ./ V1;

  for(i in 1:nObs){
    cHatObs[i] = cHat[iObs[i]];  // predictions for observed data records
  }
}
```


## <span class="section-num">3</span> Two-compartment model solved by numerical integrator for single patient {#two-compartment-model-solved-by-numerical-integrator-for-single-patient}

\label{sec:pk2cpt\_ode}
Using `pmx_solve_rk45`, the following example fits a two-compartment model
with first order absorption. User-defined function
`ode_rhs` describes the RHS of the ODEs.

```stan
functions{
  vector ode_rhs(real t, vector x, real[] parms, real[] x_r, int[] x_i){
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
```

We omit `data` and
`model` block as they are identical to Section [sec:pk2cpt](#sec:pk2cpt).

```stan
transformed data {
  row_vector[nObs] logCObs = log(cObs);
  int nTheta = 5;   // number of parameters
  int nCmt = 3;   // number of compartments
}

parameters {
  real<lower = 0> CL;
  real<lower = 0> Q;
  real<lower = 0> V1;
  real<lower = 0> V2;
  real<lower = 0> ka;
  real<lower = 0> sigma;
}

transformed parameters {
  real theta[nTheta];
  row_vector<lower = 0>[nt] cHat;
  row_vector<lower = 0>[nObs] cHatObs;
  matrix<lower = 0>[3, nt] x;

  theta[1] = CL;
  theta[2] = Q;
  theta[3] = V1;
  theta[4] = V2;
  theta[5] = ka;

  x = pmx_solve_rk45(ode_rhs, 3, time, amt, rate, ii, evid, cmt, addl, ss, theta, 1e-5, 1e-8, 1e5);

  cHat = x[2, ] ./ V1;

  for(i in 1:nObs){
    cHatObs[i] = cHat[iObs[i]];  // predictions for observed data records
  }
}

model{
```


## <span class="section-num">4</span> Joint PK-PD model {#joint-pk-pd-model}

\label{sec:fk\_model}

Neutropenia is observed in patients receiving an ME-2 drug. Our goal
is to model the relation between neutrophil counts and drug
exposure. As shown in Figure [fig:FK_model](#fig:FK_model), the Friberg-Karlsson Semi-Mechanistic model <sup id="4698e09238445a27c4bd926a20f9846e"><a href="#friberg_mechanistic_2003" title="Friberg \&amp; Karlsson, Mechanistic {Models} for {Myelosuppression}, {Investigational New Drugs}, v(2), 183--194 (2003).">friberg_mechanistic_2003</a></sup> couples
a PK model with a PD
effect to describe a delayed feedback mechanism that keeps the
absolute neutrophil count (ANC) at the
baseline in a circulatory compartment (Circ), and
the drug's effect in
reducing the proliferation rate (prol).
The delay between prol and Circ is modeled using \\(n\\) transit
comparments with mean transit time MTT = \\((n + 1)/k\_{\text{tr}}\\),
with \\(k\_{\text{tr}}\\) the transit rate constant. In the current example, we use the two compartment model in section [{{< relref "two-cpt" >}}]({{< relref "two-cpt" >}}) for
PK model, and set \\(n = 3\\).

\begin{align}
  \log(\text{ANC})& \sim N(\log(y\_{\text{circ}}), \sigma^2\_{\text{ANC}}),  \\\\\\
  y\_{\text{circ}}& = f\_{\text{FK}}(\text{MTT}, \text{Circ}\_{0}, \alpha, \gamma, c),
\end{align}

  where \\(c\\) is the drug concentration calculated from the PK model, and function \\(f\_{\text{FK}}\\) represents solving the following
nonlinear ODE for \\(y\_{\text{circ}}\\)

\begin{subequations}
  \begin{align}
  \frac{dy\_\mathrm{prol}}{dt} &= k\_\mathrm{prol} y\_\mathrm{prol} (1 - E\_\mathrm{drug})\left(\frac{\text{Circ}\_0}{y\_\mathrm{circ}}\right)^\gamma - k\_\mathrm{tr}y\_\mathrm{prol}, \\\\\\
  \frac{dy\_\mathrm{trans1}}{dt} &= k\_\mathrm{tr} y\_\mathrm{prol} - k\_\mathrm{tr} y\_\mathrm{trans1}, \\\\\\
  \frac{dy\_\mathrm{trans2}}{dt} &= k\_\mathrm{tr} y\_\mathrm{trans1} - k\_\mathrm{tr} y\_\mathrm{trans2},  \\\\\\
  \frac{dy\_\mathrm{trans3}}{dt} &= k\_\mathrm{tr} y\_\mathrm{trans2} - k\_\mathrm{tr} y\_\mathrm{trans3},  \\\\\\
  \frac{dy\_\mathrm{circ}}{dt} &= k\_\mathrm{tr} y\_\mathrm{trans3} - k\_\mathrm{tr} y\_\mathrm{circ},
   \label{eq:FK}
  \end{align}
\end{subequations}

We use \\(E\_{\text{drug}} = \alpha c\\) to model the linear effect of drug
concentration in central compartment, with
\\(c=y\_{\text{cent}}/V\_{\text{cent}}\\) based on PK solutions.

Since the ODEs specifying the Two Compartment Model
(Equation \eqref{eq:twocpt}) do not depend on the PD ODEs
(Equation \eqref{eq:FK}) and can be solved analytically
using Torsten's `pmx_solve_twocpt` function
we can specify solve the system using a coupled solver function. We do not
expect our system to be stiff and use the Runge-Kutta 4th/5th order
integrator.

<a id="orgd46b903"></a>

{{< figure src="/ox-hugo/neutrophilModel.jpg" caption="Figure 4: Friberg-Karlsson semi-mechanistic Model." >}}

The model fitting is based on simulated data

\begin{align\*}
  (\text{MTT}, \text{Circ}\_{0}, \alpha, \gamma, k\_{\text{tr}})& = (125, 5.0, 3 \times 10^{-4}, 0.17) \\\\\\
  \sigma^2\_{\text{ANC}}& = 0.001.
\end{align\*}

```stan
functions{
  vector FK_ODE(real t, vector y, vector y_pk, real[] theta, real[] rdummy, int[] idummy){
    /* PK variables */
    real VC = theta[3];

    /* PD variable */
    real mtt      = theta[6];
    real circ0    = theta[7];
    real alpha    = theta[8];
    real gamma    = theta[9];
    real ktr      = 4.0 / mtt;
    real prol     = y[1] + circ0;
    real transit1 = y[2] + circ0;
    real transit2 = y[3] + circ0;
    real transit3 = y[4] + circ0;
    real circ     = fmax(machine_precision(), y[5] + circ0);
    real conc     = y_pk[2] / VC;
    real EDrug    = alpha * conc;

    vector[5] dydt;

    dydt[1] = ktr * prol * ((1 - EDrug) * ((circ0 / circ)^gamma) - 1);
    dydt[2] = ktr * (prol - transit1);
    dydt[3] = ktr * (transit1 - transit2);
    dydt[4] = ktr * (transit2 - transit3);
    dydt[5] = ktr * (transit3 - circ);

    return dydt;
  }
}

data{
  int<lower = 1> nt;
  int<lower = 1> nObsPK;
  int<lower = 1> nObsPD;
  int<lower = 1> iObsPK[nObsPK];
  int<lower = 1> iObsPD[nObsPD];
  real<lower = 0> amt[nt];
  int<lower = 1> cmt[nt];
  int<lower = 0> evid[nt];
  real<lower = 0> time[nt];
  real<lower = 0> ii[nt];
  int<lower = 0> addl[nt];
  int<lower = 0> ss[nt];
  real rate[nt];
  vector<lower = 0>[nObsPK] cObs;
  vector<lower = 0>[nObsPD] neutObs;

  real<lower = 0> circ0Prior;
  real<lower = 0> circ0PriorCV;
  real<lower = 0> mttPrior;
  real<lower = 0> mttPriorCV;
  real<lower = 0> gammaPrior;
  real<lower = 0> gammaPriorCV;
  real<lower = 0> alphaPrior;
  real<lower = 0> alphaPriorCV;
}

transformed data{
  int nOde = 5;
  vector[nObsPK] logCObs;
  vector[nObsPD] logNeutObs;

  int nTheta = 9; // number of parameters
  int nIIV = 7; // parameters with IIV

  int n = 8;                        /* ODE dimension */
  real rtol = 1e-8;
  real atol = 1e-8;;
  int max_step = 100000;

  logCObs = log(cObs);
  logNeutObs = log(neutObs);
}

parameters{

  real<lower = 0> CL;
  real<lower = 0> Q;
  real<lower = 0> VC;
  real<lower = 0> VP;
  real<lower = 0> ka;
  real<lower = 0> mtt;
  real<lower = 0> circ0;
  real<lower = 0> alpha;
  real<lower = 0> gamma;
  real<lower = 0> sigma;
  real<lower = 0> sigmaNeut;

  // IIV parameters
  cholesky_factor_corr[nIIV] L;
  vector<lower = 0>[nIIV] omega;
}

transformed parameters{
  row_vector[nt] cHat;
  vector<lower = 0>[nObsPK] cHatObs;
  row_vector[nt] neutHat;
  vector<lower = 0>[nObsPD] neutHatObs;
  real<lower = 0> theta[nTheta];
  matrix[nOde + 3, nt] x;
  real biovar[nTheta] = rep_array(1.0, nTheta);
  real tlag[nTheta] = rep_array(0.0, nTheta);

  theta[1] = CL;
  theta[2] = Q;
  theta[3] = VC;
  theta[4] = VP;
  theta[5] = ka;
  theta[6] = mtt;
  theta[7] = circ0;
  theta[8] = alpha;
  theta[9] = gamma;

  x = pmx_solve_twocpt_rk45(FK_ODE, nOde, time, amt, rate, ii, evid, cmt, addl, ss, theta, biovar, tlag, rtol, atol, max_step);

  cHat = x[2, ] / VC;
  neutHat = x[8, ] + circ0;

  for(i in 1:nObsPK) cHatObs[i]    = cHat[iObsPK[i]];
  for(i in 1:nObsPD) neutHatObs[i] = neutHat[iObsPD[i]];
}

model {
  // Priors
  CL    ~ normal(0, 20);
  Q     ~ normal(0, 20);
  VC    ~ normal(0, 100);
  VP    ~ normal(0, 1000);
  ka    ~ normal(0, 5);
  sigma ~ cauchy(0, 1);

  mtt       ~ lognormal(log(mttPrior), mttPriorCV);
  circ0     ~ lognormal(log(circ0Prior), circ0PriorCV);
  alpha     ~ lognormal(log(alphaPrior), alphaPriorCV);
  gamma     ~ lognormal(log(gammaPrior), gammaPriorCV);
  sigmaNeut ~ cauchy(0, 1);

  // Parameters for Matt's trick
  L ~ lkj_corr_cholesky(1);
  omega ~ cauchy(0, 1);

  // observed data likelihood
  logCObs ~ normal(log(cObs), sigma);
  logNeutObs ~ normal(log(neutObs), sigmaNeut);
}
```


## <span class="section-num">5</span> Two-compartment population model {#two-compartment-population-model}

\label{sec:twocpt\_population}
Using `pmx_solve_group_bdf`, the following example fits a
two-compartment population model.

```stan
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
  int<lower = 1> np;            /* population size */
  int<lower = 1> nt;  // number of events
  int<lower = 1> nObs;  // number of observations
  int<lower = 1> iObs[nObs];  // index of observation

  // NONMEM data
  int<lower = 1> cmt[np * nt];
  int evid[np * nt];
  int addl[np * nt];
  int ss[np * nt];
  real amt[np * nt];
  real time[np * nt];
  real rate[np * nt];
  real ii[np * nt];

  real<lower = 0> cObs[np*nObs];  // observed concentration (dependent variable)
}

transformed data {
  real logCObs[np*nObs];
  int<lower = 1> len[np];
  int<lower = 1> len_theta[np];
  int<lower = 1> len_biovar[np];
  int<lower = 1> len_tlag[np];

  int nTheta = 5;   // number of parameters
  int nCmt = 3;   // number of compartments
  real biovar[np * nt, nCmt];
  real tlag[np * nt, nCmt];

  logCObs = log(cObs);

  for (id in 1:np) {
    for (j in 1:nt) {
      for (i in 1:nCmt) {
        biovar[(id - 1) * nt + j, i] = 1;
        tlag[(id - 1) * nt + j, i] = 0;
      }
    }
    len[id] = nt;
    len_theta[id] = nt;
    len_biovar[id] = nt;
    len_tlag[id] = nt;
  }
}

parameters{
  real<lower = 0> CL[np];
  real<lower = 0> Q[np];
  real<lower = 0> V1[np];
  real<lower = 0> V2[np];
  real<lower = 0> ka[np];
  real<lower = 0> sigma[np];
}

transformed parameters{
  real theta[np * nt, nTheta];
  vector<lower = 0>[nt] cHat[np];
  real<lower = 0> cHatObs[np*nObs];
  matrix[3, nt * np] x;

  for (id in 1:np) {
    for (it in 1:nt) {
      theta[(id - 1) * nt + it, 1] = CL[id];
      theta[(id - 1) * nt + it, 2] = Q[id];
      theta[(id - 1) * nt + it, 3] = V1[id];
      theta[(id - 1) * nt + it, 4] = V2[id];
      theta[(id - 1) * nt + it, 5] = ka[id];
    }
  }

  x = pmx_solve_group_bdf(twoCptModelODE, 3, len,
                          time, amt, rate, ii, evid, cmt, addl, ss,
                          theta, biovar, tlag);

  for (id in 1:np) {
    for (j in 1:nt) {
      cHat[id][j] = x[2, (id - 1) * nt + j] ./ V1[id];
    }
  }

  for (id in 1:np) {
    for(i in 1:nObs){
      cHatObs[(id - 1)*nObs + i] = cHat[id][iObs[i]];  // predictions for observed data records
    }
  }
}

model{
  // informative prior
  for(id in 1:np){
    CL[id] ~ lognormal(log(10), 0.25);
    Q[id] ~ lognormal(log(15), 0.5);
    V1[id] ~ lognormal(log(35), 0.25);
    V2[id] ~ lognormal(log(105), 0.5);
    ka[id] ~ lognormal(log(2.5), 1);
    sigma[id] ~ cauchy(0, 1);

    for(i in 1:nObs){
      logCObs[(id - 1)*nObs + i] ~ normal(log(cHatObs[(id - 1)*nObs + i]), sigma[id]);
    }
  }
}
```

When the above model is compiled with MPI support(see Section
[sec:mpi_support](#sec:mpi_support)), one can run it with within-chain
parallelization:

```bash
mpiexec -n nproc ./twocpt_population sample data file=twocpt_population.data.R init=twocpt_population.init.R
```

Here `nproc` indicates the number of parallel
processes participating ODE solution. For example, with
`np=8` for a population of 8,
`nproc=4` indicates solving 8 subjects' ODEs in
parallel, with each process solving 2 subjects.


## <span class="section-num">6</span> Lotka-Volterra group model {#lotka-volterra-group-model}

\label{sec:lotka\_volterra}
Using `pmx_integrate_ode_group_rk45`, the following example fits
a Lotka-Volterra group model, based on [Stan's case study](https://mc-stan.org/users/documentation/case-studies/lotka-volterra-predator-prey.html).

```stan
functions {
  real[] dz_dt(real t,       // time
               real[] z,     // system state {prey, predator}
               real[] theta, // parameters
               real[] x_r,   // unused data
               int[] x_i) {
    real u = z[1];
    real v = z[2];

    real alpha = theta[1];
    real beta = theta[2];
    real gamma = theta[3];
    real delta = theta[4];

    real du_dt = (alpha - beta * v) * u;
    real dv_dt = (-gamma + delta * u) * v;
    return { du_dt, dv_dt };
  }
}
data {
  int<lower = 0> N_subj;      // number of subjects
  int<lower = 0> N;           // number of measurement times
  real ts_0[N];                 // measurement times > 0
  real y0_0[2];     // initial measured populations
  real<lower = 0> y_0[N, 2];    // measured populations
}
transformed data {
  int len[N_subj] = rep_array(N, N_subj);
  real y0[N_subj, 2] = rep_array(y0_0, N_subj);
  real y[N_subj, N, 2] = rep_array(y_0, N_subj);
  real ts[N_subj * N];
  for (i in 1:N_subj) {
    ts[((i-1)*N + 1) : (i*N)] = ts_0;
  }
}

parameters {
  real<lower = 0> theta[N_subj, 4];   // { alpha, beta, gamma, delta }
  real<lower = 0> z_init[N_subj, 2];  // initial population
  real<lower = 0> sigma[N_subj, 2];   // measurement errors
}
transformed parameters {
  matrix[2, N_subj * N] z;
  z = pmx_integrate_ode_group_rk45(dz_dt, z_init, 0, len, ts, theta, rep_array(rep_array(0.0, 0), N_subj), rep_array(rep_array(0, 0),N_subj));
}
model {
  for (isub in 1:N_subj) {
    theta[isub, {1, 3}] ~ normal(1, 0.5);
    theta[isub, {2, 4}] ~ normal(0.05, 0.05);
    sigma[isub] ~ lognormal(-1, 1);
    z_init[isub] ~ lognormal(10, 1);
    for (k in 1:2) {
      y0[isub, k] ~ lognormal(log(z_init[isub, k]), sigma[isub, k]);
      y[isub, , k] ~ lognormal(log(z[k, ((isub-1)*N + 1):(isub*N)]), sigma[isub, k]);
    }
  }
```


## <span class="section-num">7</span> Univariate integral of a quadratic function {#univariate-integral-of-a-quadratic-function}

integral of a quadratic function.
This example shows how to use `univariate_integral_rk45` to calculate the
integral of a quadratic function.

```stan
functions {
  real fun_ord2(real t, real[] theta, real[] x_r, int[] x_i) {
    real a = 2.3;
    real b = 2.0;
    real c = 1.5;
    real res;
    res = a + b * t + c * t * t;
    return res;
  }
}
data {
  real t0;
  real t1;
  real dtheta[2];
  real x_r[0];
  int x_i[0];
}
transformed data {
  real univar_integral;
  univar_integral = univariate_integral_rk45(func, t0, t1, dtheta,
                          x_r, x_i);
}
/* ... */
```


## <span class="section-num">8</span> Linear intepolation {#linear-intepolation}

This example illustrates how to use `linear_intepolationi`
to fit a piecewise linear function to a data set consisting
of \\((x, y)\\) pairs.

```stan
data{
  int nObs;
  real xObs[nObs];
  real yObs[nObs];
  int nx;
  int nPred;
  real xPred[nPred];
}

transformed data{
  real xmin = min(xObs);
  real xmax = max(xObs);
}

parameters{
  real y[nx];
  real<lower = 0> sigma;
  simplex[nx - 1] xSimplex;
}

transformed parameters{
  real yHat[nObs];
  real x[nx];

  x[1] = xmin;
  x[nx] = xmax;
  for(i in 2:(nx-1))
    x[i] = x[i-1] + xSimplex[i-1] * (xmax - xmin);

  yHat = linear_interpolation(xObs, x, y);
}

model{
  xSimplex ~ dirichlet(rep_vector(1, nx - 1));
  y ~ normal(0, 25);
  yObs ~ normal(yHat, sigma);
}

generated quantities{
  real yHatPred[nPred];
  real yPred[nPred];

  yHatPred = linear_interpolation(xPred, x, y);
  for(i in 1:nPred)
    yPred[i] = normal_rng(yHatPred[i], sigma);
}
```


## <span class="section-num">9</span> Effect Compartment Population Model {#effect-compartment-population-model}

\label{sec:effcpt\_model}
Here we expand the example in [{{< relref "two-cpt" >}}]({{< relref "two-cpt" >}}) to a population model fitted to the
combined data from phase I and phase IIa studies. The
parameters exhibit inter-individual variations (IIV), due to
both random effects and to the patients' body weight,
treated as a covariate and denoted \\(bw\\).


### <span class="section-num">9.1</span> Population Model for Plasma Drug Concentration \\(c\\) {#population-model-for-plasma-drug-concentration--c}

\begin{gather\*}
  \log\left(c\_{ij}\right) \sim N\left(\log\left(\widehat{c}\_{ij}\right),\sigma^2\right), \\\\\\
  \widehat{c}\_{ij} = f\_{2cpt}\left(t\_{ij},D\_j,\tau\_j,CL\_j,Q\_j,V\_{1j},V\_{2j},k\_{aj}\right), \\\\\\
  \log\left(CL\_j,Q\_j,V\_{ssj},k\_{aj}\right) \sim N\left(\log\left(\widehat{CL}\left(\frac{bw\_j}{70}\right)^{0.75},\widehat{Q}\left(\frac{bw\_j}{70}\right)^{0.75}, \widehat{V}\_{ss}\left(\frac{bw\_j}{70}\right),\widehat{k}\_a\right),\Omega\right), \\\\\\
  V\_{1j} = f\_{V\_1}V\_{ssj}, \\\\\\
  V\_{2j} = \left(1 - f\_{V\_1}\right)V\_{ssj}, \\\\\\
  \left(\widehat{CL},\widehat{Q},\widehat{V}\_{ss},\widehat{k}\_a, f\_{V\_1}\right) = \left(10\ {\rm L/h},15\  {\rm L/h},140\  {\rm L},2\ {\rm h^{-1}}, 0.25 \right), \\\\\\
  \Omega = \left(\begin{array}{cccc} 0.25^2 & 0 & 0 & 0 \\ 0 & 0.25^2 & 0 & 0 \\\\\\
                    0 & 0 & 0.25^2 & 0 \\ 0 & 0 & 0 & 0.25^2  \end{array}\right), \\\\\\
  \sigma = 0.1
\end{gather\*}

Furthermore we add a fourth compartment in which we measure
a PD effect(Figure [5](#org3c0289f)).

<a id="org3c0289f"></a>

{{< figure src="/ox-hugo/effCptModel.png" caption="Figure 5: Effect Compartment Model" >}}


### <span class="section-num">9.2</span> Effect Compartment Model for PD response \\(R\\). {#effect-compartment-model-for-pd-response--r--dot}

\begin{gather\*}
R\_{ij} \sim N\left(\widehat{R}\_{ij},\sigma\_{R}^2\right), \\\\\\
\widehat{R}\_{ij} = \frac{E\_{max}c\_{eij}}{EC\_{50j} + c\_{eij}}, \\\\\\
c\_{e\cdot j}^\prime = k\_{e0j}\left(c\_{\cdot j} - c\_{e\cdot j}\right), \\\\\\
\log\left(EC\_{50j}, k\_{e0j}\right) \sim N\left(\log\left(\widehat{EC}\_{50}, \widehat{k}\_{e0}\right),\Omega\_R\right), \\\\\\
\left(E\_{max}, \widehat{EC}\_{50},\widehat{k}\_{e0}\right) = \left(100, 100.7, 1\right), \\\\\\
\Omega\_R = \left(\begin{array}{cc} 0.2^2 & 0 \\ 0 & 0.25^2  \end{array}\right), \ \ \ \sigma\_R = 10.
\end{gather\*}

The PK and the PD data are simulated using the following
treatment.

-   Phase I study
    -   Single dose and multiple doses
    -   Parallel dose escalation design
    -   25 subjects per dose
    -   Single doses: 5, 10, 20, and 40 mg
    -   PK: plasma concentration of parent drug (\\(c\\))
    -   PD response: Emax function of effect compartment concentration (\\(R\\))
    -   PK and PD measured at 0.125, 0.25, 0.5, 0.75, 1, 2, 3, 4, 6, 8, 12, 18, and 24 hours
-   Phase IIa trial in patients
    -   100 subjects
    -   Multiple doses: 20 mg
    -   sparse PK and PD data (3-6 samples per patient)

The model is simultaneously fitted to the PK and the PD
data. For this effect compartment model, we construct a
constant rate matrix and use `pmx_solve_linode`. Correct use of
Torsten requires the user pass the entire event history
(observation and dosing events) for an individual to the
function. Thus the Stan model shows the call to `pmx_solve_linode`
within a loop over the individual subjects rather than over
the individual observations. Note that the correlation matrix \\(\rho\\) does not explicitly appear
in the model, but it is used to construct \\(\Omega\\), which parametrizes
the PK IIV.

```stan
data{
  int<lower = 1> nSubjects;
  int<lower = 1> nt;
  int<lower = 1> nObs;
  int<lower = 1> iObs[nObs];
  real<lower = 0> amt[nt];
  int<lower = 1> cmt[nt];
  int<lower = 0> evid[nt];
  int<lower = 1> start[nSubjects];
  int<lower = 1> end[nSubjects];
  real<lower = 0> time[nt];
  vector<lower = 0>[nObs] cObs;
  vector[nObs] respObs;
  real<lower = 0> weight[nSubjects];
  real<lower = 0> rate[nt];
  real<lower = 0> ii[nt];
  int<lower = 0> addl[nt];
  int<lower = 0> ss[nt];
}

transformed data{
  vector[nObs] logCObs = log(cObs);
  int<lower = 1> nRandom = 5;
  int nCmt = 4;
  real biovar[nCmt] = rep_array(1.0, nCmt);
  real tlag[nCmt] = rep_array(0.0, nCmt);
}

parameters{
  real<lower = 0> CLHat;
  real<lower = 0> QHat;
  real<lower = 0> V1Hat;
  real<lower = 0> V2Hat;
  //  real<lower = 0> kaHat;
  real<lower = (CLHat / V1Hat + QHat / V1Hat + QHat / V2Hat +
                sqrt((CLHat / V1Hat + QHat / V1Hat + QHat / V2Hat)^2 -
                     4 * CLHat / V1Hat * QHat / V2Hat)) / 2> kaHat; // ka > lambda_1
  real<lower = 0> ke0Hat;
  real<lower = 0> EC50Hat;
  vector<lower = 0>[nRandom] omega;
  corr_matrix[nRandom] rho;
  real<lower = 0> omegaKe0;
  real<lower = 0> omegaEC50;
  real<lower = 0> sigma;
  real<lower = 0> sigmaResp;

  // reparameterization
  vector[nRandom] logtheta_raw[nSubjects];
  real logKe0_raw[nSubjects];
  real logEC50_raw[nSubjects];
}

transformed parameters{
  vector<lower = 0>[nRandom] thetaHat;
  cov_matrix[nRandom] Omega;
  real<lower = 0> CL[nSubjects];
  real<lower = 0> Q[nSubjects];
  real<lower = 0> V1[nSubjects];
  real<lower = 0> V2[nSubjects];
  real<lower = 0> ka[nSubjects];
  real<lower = 0> ke0[nSubjects];
  real<lower = 0> EC50[nSubjects];
  matrix[nCmt, nCmt] K;
  real k10;
  real k12;
  real k21;
  row_vector<lower = 0>[nt] cHat;
  row_vector<lower = 0>[nObs] cHatObs;
  row_vector<lower = 0>[nt] respHat;
  row_vector<lower = 0>[nObs] respHatObs;
  row_vector<lower = 0>[nt] ceHat;
  matrix[nCmt, nt] x;

  matrix[nRandom, nRandom] L;
  vector[nRandom] logtheta[nSubjects];
  real logKe0[nSubjects];
  real logEC50[nSubjects];

  thetaHat[1] = CLHat;
  thetaHat[2] = QHat;
  thetaHat[3] = V1Hat;
  thetaHat[4] = V2Hat;
  thetaHat[5] = kaHat;

  Omega = quad_form_diag(rho, omega); // diag_matrix(omega) * rho * diag_matrix(omega)
  L = cholesky_decompose(Omega);

  for(j in 1:nSubjects){
    logtheta[j] = log(thetaHat) + L * logtheta_raw[j];
    logKe0[j] = log(ke0Hat) + logKe0_raw[j] * omegaKe0;
    logEC50[j] = log(EC50Hat) + logEC50_raw[j] * omegaEC50;

    CL[j] = exp(logtheta[j, 1]) * (weight[j] / 70)^0.75;
    Q[j] = exp(logtheta[j, 2]) * (weight[j] / 70)^0.75;
    V1[j] = exp(logtheta[j, 3]) * weight[j] / 70;
    V2[j] = exp(logtheta[j, 4]) * weight[j] / 70;
    ka[j] = exp(logtheta[j, 5]);
    ke0[j] = exp(logKe0[j]);
    EC50[j] = exp(logEC50[j]);

    k10 = CL[j] / V1[j];
    k12 = Q[j] / V1[j];
    k21 = Q[j] / V2[j];

    K = rep_matrix(0, nCmt, nCmt);

    K[1, 1] = -ka[j];
    K[2, 1] = ka[j];
    K[2, 2] = -(k10 + k12);
    K[2, 3] = k21;
    K[3, 2] = k12;
    K[3, 3] = -k21;
    K[4, 2] = ke0[j];
    K[4, 4] = -ke0[j];

    x[, start[j]:end[j] ] = pmx_solve_linode(time[start[j]:end[j]], amt[start[j]:end[j]],
                                             rate[start[j]:end[j]], ii[start[j]:end[j]],
                                             evid[start[j]:end[j]], cmt[start[j]:end[j]],
                                             addl[start[j]:end[j]], ss[start[j]:end[j]], K, biovar, tlag);

    cHat[start[j]:end[j]] = 1000 * x[2, start[j]:end[j]] ./ V1[j];
    ceHat[start[j]:end[j]] = 1000 * x[4, start[j]:end[j]] ./ V1[j];
    respHat[start[j]:end[j]] = 100 * ceHat[start[j]:end[j]] ./
       (EC50[j] + ceHat[start[j]:end[j]]);
  }

  cHatObs = cHat[iObs];
  respHatObs = respHat[iObs];
}

model{
  // Prior
  CLHat ~ lognormal(log(10), 0.2);
  QHat ~ lognormal(log(15), 0.2);
  V1Hat ~ lognormal(log(30), 0.2);
  V2Hat ~ lognormal(log(100), 0.2);
  kaHat ~ lognormal(log(5), 0.25);
  ke0Hat ~ lognormal(log(10), 0.25);
  EC50Hat ~ lognormal(log(1.0), 0.2);
  omega ~ normal(0, 0.2);
```


### <span class="section-num">9.3</span> Results {#results}

We use the same diagnosis tools as for the
previous examples. Table [effCptModelParms](#effCptModelParms) summarises the
statistics and diagnostics of the parameters. In particular, `rhat`
for all parameters being close to 1.0 indicates convergence of the 4
chains. Figure [effcpt_mcmc_density](#effcpt_mcmc_density) shows the posterior density of
the parameters.

Posterior prediction check (PPC) in Figure
[effcpt_ppc_5mg](#effcpt_ppc_5mg)-[effcpt_ppc_study_2_20mg](#effcpt_ppc_study_2_20mg) show that the fits to the plasma concentration
are in close agreement with the data, notably for the sparse data case (phase IIa study). The fits
to the PD data (Figure
[effcpt_ppc_resp_5mg](#effcpt_ppc_resp_5mg)-[effcpt_ppc_resp_study_2_20mg](#effcpt_ppc_resp_study_2_20mg)) look
reasonable considering data being more noisy so the model
produces larger credible intervals. Both the summary table and PPC
plots show that the estimated values of the
parameters are consistent with the values used to simulate the data.

<a id="table--effCptModelParms"></a>
<div class="table-caption">
  <span class="table-number"><a href="#table--effCptModelParms">Table 1</a></span>:
  Summary of the MCMC simulations of the marginal posterior distributions of the model parameters for the effect compartment model example.
</div>

| variable  | mean    | median  | sd    | mad   | q5     | q95     | rhat  | ess\_bulk | ess\_tail |
|-----------|---------|---------|-------|-------|--------|---------|-------|-----------|-----------|
| CLHat     | 10.121  | 10.120  | 0.195 | 0.192 | 9.797  | 10.445  | 1.007 | 319.942   | 630.619   |
| QHat      | 14.858  | 14.853  | 0.347 | 0.344 | 14.301 | 15.432  | 1.000 | 1106.126  | 1712.821  |
| V1Hat     | 34.493  | 34.516  | 1.004 | 0.995 | 32.814 | 36.086  | 1.004 | 671.777   | 1563.396  |
| V2Hat     | 103.269 | 103.291 | 2.876 | 2.878 | 98.568 | 108.019 | 1.002 | 1689.165  | 2580.382  |
| kaHat     | 1.968   | 1.969   | 0.076 | 0.074 | 1.843  | 2.087   | 1.001 | 1204.531  | 1747.427  |
| ke0Hat    | 1.102   | 1.100   | 0.046 | 0.045 | 1.030  | 1.180   | 1.001 | 4008.337  | 3167.030  |
| EC50Hat   | 99.512  | 99.542  | 2.124 | 2.098 | 95.981 | 102.987 | 1.000 | 2557.436  | 2773.519  |
| omega[1]  | 0.268   | 0.267   | 0.016 | 0.016 | 0.242  | 0.295   | 1.008 | 594.842   | 978.297   |
| omega[2]  | 0.229   | 0.228   | 0.021 | 0.021 | 0.195  | 0.264   | 1.002 | 1245.453  | 1966.911  |
| omega[3]  | 0.212   | 0.211   | 0.029 | 0.029 | 0.165  | 0.261   | 1.005 | 623.820   | 1692.248  |
| omega[4]  | 0.263   | 0.262   | 0.026 | 0.026 | 0.221  | 0.306   | 1.002 | 1396.611  | 2260.425  |
| omega[5]  | 0.272   | 0.271   | 0.036 | 0.035 | 0.217  | 0.335   | 1.008 | 293.132   | 728.867   |
| rho[1,2]  | 0.197   | 0.200   | 0.100 | 0.101 | 0.029  | 0.360   | 1.003 | 1322.261  | 1955.862  |
| rho[1,3]  | -0.161  | -0.161  | 0.122 | 0.121 | -0.361 | 0.042   | 1.001 | 1609.160  | 2270.515  |
| rho[1,4]  | -0.101  | -0.105  | 0.107 | 0.107 | -0.270 | 0.083   | 1.001 | 1685.591  | 2353.498  |
| rho[1,5]  | 0.016   | 0.015   | 0.128 | 0.128 | -0.192 | 0.226   | 1.000 | 2039.767  | 2939.988  |
| rho[2,3]  | 0.091   | 0.092   | 0.144 | 0.148 | -0.143 | 0.328   | 1.008 | 718.187   | 1550.836  |
| rho[2,4]  | 0.186   | 0.190   | 0.125 | 0.125 | -0.025 | 0.384   | 1.005 | 948.704   | 1819.199  |
| rho[2,5]  | 0.146   | 0.145   | 0.157 | 0.161 | -0.111 | 0.402   | 1.003 | 626.620   | 1546.157  |
| rho[3,4]  | 0.815   | 0.827   | 0.093 | 0.094 | 0.646  | 0.947   | 1.010 | 309.098   | 736.635   |
| rho[3,5]  | -0.318  | -0.323  | 0.219 | 0.228 | -0.678 | 0.055   | 1.016 | 200.806   | 607.958   |
| rho[4,5]  | -0.295  | -0.299  | 0.161 | 0.162 | -0.551 | -0.019  | 1.008 | 546.998   | 1151.092  |
| omegaKe0  | 0.265   | 0.265   | 0.047 | 0.047 | 0.188  | 0.346   | 1.001 | 1731.276  | 2049.892  |
| omegaEC50 | 0.216   | 0.216   | 0.020 | 0.020 | 0.182  | 0.249   | 1.001 | 1599.567  | 1844.056  |
| sigma     | 0.099   | 0.099   | 0.002 | 0.002 | 0.095  | 0.103   | 1.002 | 1726.283  | 2836.027  |
| sigmaResp | 10.165  | 10.166  | 0.198 | 0.198 | 9.844  | 10.495  | 1.002 | 4788.527  | 2923.203  |

<a id="org41d7191"></a>

</ox-hugo/density.pdf>

<a id="org24095d5"></a>

</ox-hugo/ppc_study_1_5mg.pdf>

<a id="org787a63e"></a>

</ox-hugo/ppc_study_1_10mg.pdf>

<a id="orgaae9f95"></a>

</ox-hugo/ppc_study_1_20mg.pdf>

<a id="org664a06f"></a>

</ox-hugo/ppc_study_1_40mg.pdf>

<a id="org9112827"></a>

</ox-hugo/ppc_study_2_20mg.pdf>

<a id="orgcb83f04"></a>

</ox-hugo/ppc_study_1_5mg_resp.pdf>

<a id="org2436410"></a>

</ox-hugo/ppc_study_1_10mg_resp.pdf>

<a id="orgd08227a"></a>

</ox-hugo/ppc_study_1_20mg_resp.pdf>

<a id="org096dec3"></a>

</ox-hugo/ppc_study_1_40mg_resp.pdf>

<a id="org8671abd"></a>

</ox-hugo/ppc_study_2_20mg_resp.pdf>


## <span class="section-num">10</span> Friberg-Karlsson Semi-Mechanistic Population Model {#friberg-karlsson-semi-mechanistic-population-model}

\label{sec:fkpop\_model}
We now return to the example in Section [sec:fk_model](#sec:fk_model) and extend
it to a population model. While we recommend using the coupled
solver, and this time we solve it using group solver. We leave it
as an exercise to the reader to rewrite the model with
coupled solver.


### <span class="section-num">10.1</span> Friberg-Karlsson Population Model for drug-induced myelosuppression (\\(ANC\\)) {#friberg-karlsson-population-model-for-drug-induced-myelosuppression--anc}

\begin{gather\*}
\log(ANC\_{ij}) \sim N(Circ\_{ij}, \sigma^2\_{ANC}), \\\\\\
\log\left(MTT\_j, Circ\_{0j}, \alpha\_j\right) \sim N\left(\log\left(\widehat{MTT}, \widehat{Circ\_0}, \widehat{\alpha}\right), \Omega\_{ANC}\right), \\\\\\
\left(\widehat{MTT}, \widehat{Circ}\_0,\widehat{\alpha}, \gamma \right) = \left(125, 5, 2, 0.17\right), \\\\\\
\Omega\_{ANC} = \left(\begin{array}{ccc} 0.2^2 & 0 & 0 \\ 0 & 0.35^2 & 0 \\ 0 & 0 & 0.2^2 \end{array}\right), \\\\\\
\sigma\_{ANC} = 0.1, \\\\\\
\Omega\_{PK} = \left(\begin{array}{ccccc} 0.25^2 & 0 &a 0 & 0 & 0 \\ 0 & 0.4^2 & 0 & 0 & 0 \\\\\\
0 & 0 & 0.25^2 & 0 & 0 \\ 0 & 0 & 0 & 0.4^2 & 0 \\ 0 & 0 & 0 & 0 & 0.25^2  \end{array}\right)
\end{gather\*}

The PK and the PD data are simulated using the following treatment.

-   Phase IIa trial in patients
    -   Multiple doses: 80,000 mg
    -   Parallel dose escalation design
    -   15 subjects
    -   PK: plasma concentration of parent drug (\\(c\\))
    -   PD response: Neutrophil count (\\(ANC\\))
    -   PK measured at 0.083, 0.167, 0.25, 0.5, 0.75, 1, 2, 3, 4, 6, 8, 12, 18, and 24 hours
    -   PD measured once every two days for 28 days.

Once again, we simultaneously fit the model to the PK and the PD
data. It pays off to construct informative priors. For instance, we could
fit the PK data first, as was done in  example 1, and get informative
priors on the PK parameters. The PD parameters are drug independent,
so we can use information from the neutropenia literature. In this
example, we choose to use strongly informative priors on both PK and PD
parameters.

The ODE is defined as

```stan
functions{
    vector twoCptNeutModelODE(real t, vector x, real[] parms, real[] rdummy, int[] idummy){
    real k10;
    real k12;
    real k21;
    real CL;
    real Q;
    real V1;
    real V2;
    real ka;
    real mtt;
    real circ0;
    real gamma;
    real alpha;
    real ktr;
    vector[8] dxdt;
    real conc;
    real EDrug;
    real transit1;
    real transit2;
    real transit3;
    real circ;
    real prol;

    CL = parms[1];
    Q = parms[2];
    V1 = parms[3];
    V2 = parms[4];
    ka = parms[5];
    mtt = parms[6];
    circ0 = parms[7];
    gamma = parms[8];
    alpha = parms[9];

    k10 = CL / V1;
    k12 = Q / V1;
    k21 = Q / V2;

    ktr = 4 / mtt;

    dxdt[1] = -ka * x[1];
    dxdt[2] = ka * x[1] - (k10 + k12) * x[2] + k21 * x[3];
    dxdt[3] = k12 * x[2] - k21 * x[3];
    conc = x[2]/V1;
    EDrug = alpha * conc;
    // x[4], x[5], x[6], x[7] and x[8] are differences from circ0.
    prol = x[4] + circ0;
    transit1 = x[5] + circ0;
    transit2 = x[6] + circ0;
    transit3 = x[7] + circ0;
    circ = fmax(machine_precision(), x[8] + circ0); // Device for implementing a modeled
                                                    // initial condition
    dxdt[4] = ktr * prol * ((1 - EDrug) * ((circ0 / circ)^gamma) - 1);
    dxdt[5] = ktr * (prol - transit1);
    dxdt[6] = ktr * (transit1 - transit2);
    dxdt[7] = ktr * (transit2 - transit3);
    dxdt[8] = ktr * (transit3 - circ);

    return dxdt;
  }
}
```

We use the `pmx_solve_group_rk45` function to
solve the entire population's events.

```stan
transformed parameters{
  row_vector[nt] cHat;
  vector[nObsPK] cHatObs;
  row_vector[nt] neutHat;
  vector[nObsPD] neutHatObs;
  matrix[8, nt] x;
  real<lower = 0> parms[nSubjects, nTheta]; // The [1] indicates the parameters are constant

  // variables for Matt's trick
  vector<lower = 0>[nIIV] thetaHat;
  matrix<lower = 0>[nSubjects, nIIV] thetaM;

  // Matt's trick to use unit scale
  thetaHat[1] = CLHat;
  thetaHat[2] = QHat;
  thetaHat[3] = V1Hat;
  thetaHat[4] = V2Hat;
  thetaHat[5] = mttHat;
  thetaHat[6] = circ0Hat;
  thetaHat[7] = alphaHat;
  thetaM = (rep_matrix(thetaHat, nSubjects) .*
             exp(diag_pre_multiply(omega, L * etaStd)))';

  for(i in 1:nSubjects) {
    parms[i, 1] = thetaM[i, 1] * (weight[i] / 70)^0.75; // CL
    parms[i, 2] = thetaM[i, 2] * (weight[i] / 70)^0.75; // Q
    parms[i, 3] = thetaM[i, 3] * (weight[i] / 70); // V1
    parms[i, 4] = thetaM[i, 4] * (weight[i] / 70); // V2
    parms[i, 5] = kaHat; // ka
    parms[i, 6] = thetaM[i, 5]; // mtt
    parms[i, 7] = thetaM[i, 6]; // circ0
    parms[i, 8] = gamma;
    parms[i, 9] = thetaM[i, 7]; // alpha
  }

  /* group solver */
  x = pmx_solve_group_rk45(twoCptNeutModelODE, 8, len,
                           time, amt, rate, ii, evid, cmt, addl, ss,
                           parms,
                           1e-6, 1e-6, 500);

  for(i in 1:nSubjects) {
    cHat[start[i]:end[i]] = x[2, start[i]:end[i]] / parms[i, 3]; // divide by V1
    neutHat[start[i]:end[i]] = x[8, start[i]:end[i]] + parms[i, 7]; // Add baseline
  }

  for(i in 1:nObsPK) cHatObs[i] = cHat[iObsPK[i]];
  for(i in 1:nObsPD) neutHatObs[i] = neutHat[iObsPD[i]];
}
```

This allows us to use within-chain paralleleisation to reduce
simulation time. When run from cmdstan, each MPI run generates one
chain, and we use 4 MPI runs to generate 4 chains.

```bash
# chain 1
mpiexec -n nproc ./FribergKarlsson sample adapt delta=0.95 data file=fribergkarlsson.data.R init=fribergkarlsson.init.R random seed=8765 id=1 output file=output.1.csv
# chain 2
mpiexec -n nproc ./FribergKarlsson sample adapt delta=0.95 data file=fribergkarlsson.data.R init=fribergkarlsson.init.R random seed=8765 id=2 output file=output.2.csv
# chain 3
mpiexec -n nproc ./FribergKarlsson sample adapt delta=0.95 data file=fribergkarlsson.data.R init=fribergkarlsson.init.R random seed=8765 id=3 output file=output.3.csv
# chain 4
mpiexec -n nproc ./FribergKarlsson sample adapt delta=0.95 data file=fribergkarlsson.data.R init=fribergkarlsson.init.R random seed=8765 id=4 output file=output.4.csv
```


### <span class="section-num">10.2</span> Results {#results}

Table [FkpopModelParms](#FkpopModelParms) summarizes the sampling and some diagnostics output.
estimation reflects the real value of the parameters (Table [FkpopModelParms](#FkpopModelParms) and Figure [fkpop_mcmc_density](#fkpop_mcmc_density).
Similar to the previous example, PPCs shown in Figure [fkpop_ppc_pk](#fkpop_ppc_pk)
and [fkpop_ppc_pd](#fkpop_ppc_pd) indicate the model is a good fit.

<a id="table--FkpopModelParms"></a>
<div class="table-caption">
  <span class="table-number"><a href="#table--FkpopModelParms">Table 2</a></span>:
  Summary of the MCMC simulations of the marginal posterior distributions of the model parameters for the Friberg-Karlsson population model example.
</div>

| variable  | mean    | median  | sd       | mad      | q5      | q95     | rhat  | ess\_bulk | ess\_tail |
|-----------|---------|---------|----------|----------|---------|---------|-------|-----------|-----------|
| CLHat     | 9.539   | 9.535   | 0.522    | 0.487    | 8.692   | 10.401  | 1.006 | 971.369   | 1655.449  |
| QHat      | 15.401  | 15.386  | 1.018    | 1.000    | 13.742  | 17.090  | 1.000 | 2263.843  | 2447.006  |
| V1Hat     | 37.396  | 37.360  | 2.244    | 2.228    | 33.762  | 41.058  | 1.001 | 1936.476  | 2372.815  |
| V2Hat     | 101.698 | 101.394 | 6.503    | 6.119    | 91.529  | 112.538 | 1.001 | 2580.227  | 2592.925  |
| kaHat     | 1.997   | 1.997   | 0.074    | 0.074    | 1.873   | 2.115   | 1.001 | 7056.877  | 2993.406  |
| mttHat    | 113.681 | 113.204 | 11.506   | 10.910   | 95.807  | 133.514 | 1.001 | 4255.900  | 3269.646  |
| circ0Hat  | 4.760   | 4.752   | 0.241    | 0.229    | 4.375   | 5.163   | 1.002 | 3774.920  | 2783.663  |
| omega[1]  | 0.223   | 0.217   | 0.047    | 0.042    | 0.160   | 0.307   | 1.000 | 1751.864  | 2235.607  |
| omega[2]  | 0.339   | 0.329   | 0.073    | 0.067    | 0.239   | 0.473   | 1.001 | 2363.843  | 2607.056  |
| omega[3]  | 0.264   | 0.256   | 0.057    | 0.051    | 0.186   | 0.367   | 1.002 | 2128.660  | 2018.425  |
| omega[4]  | 0.257   | 0.249   | 0.056    | 0.051    | 0.182   | 0.361   | 1.003 | 2293.877  | 2937.673  |
| omega[5]  | 0.177   | 0.169   | 0.112    | 0.118    | 0.019   | 0.376   | 1.000 | 1550.483  | 2045.025  |
| omega[6]  | 0.188   | 0.183   | 0.044    | 0.041    | 0.127   | 0.269   | 1.000 | 2377.698  | 2965.713  |
| omega[7]  | 0.409   | 0.394   | 0.256    | 0.259    | 0.045   | 0.865   | 1.003 | 1386.987  | 2015.873  |
| gamma     | 0.171   | 0.168   | 0.035    | 0.033    | 0.121   | 0.235   | 1.000 | 8809.668  | 3189.676  |
| sigma     | 0.097   | 0.096   | 0.003    | 0.003    | 0.093   | 0.101   | 1.002 | 5436.508  | 2899.706  |
| sigmaNeut | 0.106   | 0.105   | 0.012    | 0.011    | 0.088   | 0.127   | 1.000 | 2809.059  | 3031.605  |
| alphaHat  | 2.24e-4 | 2.19e-4 | 3.97e-05 | 3.80e-05 | 1.66e-4 | 2.96e-4 | 1.000 | 5138.105  | 2807.328  |

<a id="org017c598"></a>

</ox-hugo/density.pdf>

<a id="orgd630694"></a>

</ox-hugo/ppc_pk.pdf>

<a id="orge7cae4e"></a>

</ox-hugo/ppc_pd.pdf>

\appendix
