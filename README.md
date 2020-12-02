# Torsten


## Overview

Torsten is a collection of Stan functions to facilitate analysis of pharmacometric data using Stan. The current version includes:

-   Specific linear compartment models:
    -   One compartment model with first order absorption.
    -   Two compartment model with elimination from and first order absorption into central compartment
-   General linear compartment model described by a system of first-order \underline{linear} Ordinary Differential Equations (ODEs).
-   General compartment model described by a system of first order ODEs.
-   Mix compartment model with PK forcing function described by a linear one or two compartment model.

The models and data format are based on NONMEM\textregistered{}<sup><a id="fnr.1" class="footref" href="#fn.1">1</a></sup>/NMTRAN/PREDPP conventions including:

-   Recursive calculation of model predictions
    -   This permits piecewise constant covariate values
-   Bolus or constant rate inputs into any compartment
-   Handles single dose and multiple dose histories
-   Handles steady state dosing histories
    -   Note: The infusion time must be shorter than the inter-dose interval.
-   Implemented NMTRAN data items include: TIME, EVID, CMT, AMT, RATE, ADDL, II, SS

In general, all real variables may be passed as model parameters. A few exceptions apply /to functions which use a numerical integrator(i.e. the general and the mix compartment models). The below listed cases present technical difficulties, which we expect to overcome in Torsten's next release:

-   In the case of a multiple truncated infusion rate dosing regimen:
    -   The bioavailability (F) and the amount (AMT) must be fixed.

This library provides Stan language functions that calculate amounts in each compartment, given an event schedule and an ODE system.


## Installation

Currently Torsten is based on a forked version of Stan. The latest v0.88 is compatible with Stan v2.25.0. Torsten can be accessed from command line for cmdstan interface and =cmdstanr=(<https://mc-stan.org/cmdstanr/>) for R interface.


### Command line interface

After downloading the project

-   <https://github.com/metrumresearchgroup/Torsten>

The command line interface `cmdstan` is available to use without installation. The following command builds a Torsten model `model_name` in `model_path`

```sh
cd Torsten/cmdstan; make model_path/model_name
```


### R interface

After installing cmdstanr from <https://mc-stan.org/cmdstanr/>, use the following command to set path

```r
cmdstanr::set_cmdstan_path("Torsten/cmdstan")
```

Then one can follow

<https://mc-stan.org/cmdstanr/articles/cmdstanr.html>

to compile and run Torsten models.

-   MPI support

    Torsten's MPI support is of a different flavour than `reduce_sum` found in Stan. To be able to utilize MPI parallelisation, one first needs to ensure an MPI library such as
    
    -   <https://www.mpich.org/downloads/>
    -   <https://www.open-mpi.org/software/ompi/>
    
    is available. Torsen's implementation is tested on both `MPICH` and `OpenMPI`.
    
    To use MPI-supported population/group solvers, add/edit `make/local`
    
    ```sh
    TORSTEN_MPI=1
    
    # path to MPI headers
    CXXFLAGS += -isystem /usr/local/include
    # if you are using Metrum's metworx platform, add MPICH3's
    # headers with
    # CXXFLAGS += -isystem /usr/local/mpich3/include
    ```
    
    Note that currently `TORSTEN_MPI` and `STAN_MPI` flags conflict on processes management and cannot be used in a same Stan model, and MPI support is only available through `cmdstan` interface.


### Testing

Models in `example-models` directory are for tutorial and demonstration. The following shows how to build and run the two-compartment model using `cmdstanr`, and use `bayesplot` to examine posterior density of `CL`.

```r
library("cmdstanr")
set_cmdstan_path("Torsten/cmdstan")
file.dir <- file.path("Torsten", "example-models", "pk2cpt")
file  <- file.path(file.dir, "pk2cpt.stan")
model <- cmdstan_model(file)
fit <- model$sample(data = file.path(file.dir, "pk2cpt.data.R"),
                    init = file.path(file.dir, "pk2cpt.init.R"),
                    seed = 123,
                    chains = 4,
                    parallel_chains = 2,
                    refresh = 500)
bayesplot::mcmc_dens_overlay(fit$draws("CL"))
```


## Development plans

Our current plans for future development of Torsten include the following:

-   Build a system to easily share packages of Stan functions (written in C++ or in the Stan language)
-   Optimize Matrix exponential functions
    -   Function for the action of Matrix Exponential on a vector
    -   Hand-coded gradients
    -   Special algorithm for matrices with special properties
-   Develop new method for large-scale hierarchical models with costly ODE solving.
-   Fix issue that arises when computing the adjoint of the lag time parameter (in a dosing compartment) evaluated at \(t_{\text{lag}} = 0\).
-   Extend formal tests
    -   More unit tests and better CD/CI support.
    -   Comparison with simulations from the R package `mrgsolve` and the software NONMEM\textregistered{}
    -   Recruit non-developer users to conduct beta testing


## Changelog


### Version 0.88 <span class="timestamp-wrapper"><span class="timestamp">&lt;2020-12-01 Tue&gt;</span></span>

-   Added

    -   Experimental feature of cross-chain warmup
    -   Bioavailability, lag time, ODE real & integer data are optional in PMX function signatures.
    -   Support all EVID options from NM-TRAN and mrgsolve.

-   Changed

    -   More efficient memory management of COVDES implenmentation.
    -   Update of MPI framework to adapt multilevel paralleism.
    -   Update to Stan version 2.25.0.
    -   Use cmdstanr as R interface.
    -   Stop supporting rstan as R interface.


### Version 0.87 <span class="timestamp-wrapper"><span class="timestamp">&lt;2019-07-26 Fri&gt;</span></span>

-   Added

    -   MPI dynamic load balance for Torsten's population ODE integrators
        
        -   `pmx_integrate_ode_group_adams`
        -   `pmx_integrate_ode_group_bdf`
        -   `pmx_integrate_ode_group_rk45`
        
        To invoke dynamic load balance instead of default static balance for MPI, issue `TORSTEN_MPI=2` in `make/local`.
    -   Support `RATE` as parameter in `pmx_solve_rk45/bdf/adams` functions.

-   Changed

    -   Some fixes on steady-state solvers
    -   Update to rstan version 2.19.2.


### Version 0.86 <span class="timestamp-wrapper"><span class="timestamp">&lt;2019-05-15 Wed&gt;</span></span>

-   Added

    -   Torsten's ODE integrator functions
        
        -   `pmx_integrate_ode_adams`
        -   `pmx_integrate_ode_bdf`
        -   `pmx_integrate_ode_rk45`
        
        and their counterparts to solve a population/group of subjects governed by an ODE
        
        -   `pmx_integrate_ode_group_adams`
        -   `pmx_integrate_ode_group_bdf`
        -   `pmx_integrate_ode_group_rk45`
    -   Torsten's population PMX solver functions for general ODE models
        -   `pmx_solve_group_adams`
        -   `pmx_solve_group_bdf`
        -   `pmx_solve_group_rk45`
    -   Support time step `ts` as parameter in `pmx_integrate_ode_xxx` solvers.

-   Changed

    -   Renaming Torsten functions in previous releases, the old-new name mapping is
        
        -   `PKModelOneCpt` &rarr; `pmx_solve_onecpt`
        -   `PKModelTwoCpt` &rarr; `pmx_solve_onecpt`
        -   `linOdeModel` &rarr; `pmx_solve_linode`
        -   `generalOdeModel_adams` &rarr; `pmx_solve_adams`
        -   `generalOdeModel_bdf` &rarr; `pmx_solve_bdf`
        -   `generalOdeModel_rk45` &rarr; `pmx_solve_rk45`
        -   `mixOde1CptModel_bdf` &rarr; `pmx_solve_onecpt_bdf`
        -   `mixOde1CptModel_rk45` &rarr; `pmx_solve_onecpt_rk45`
        -   `mixOde2CptModel_bdf` &rarr; `pmx_solve_twocpt_bdf`
        -   `mixOde2CptModel_rk45` &rarr; `pmx_solve_twocpt_rk45`
        
        Note that the new version of the above functions return the *transpose* of the matrix returned by the old versions, in order to improve memory efficiency. The old version are retained but will be deprecated in the future.
    -   Update to Stan version 2.19.1.


### Version 0.85 <span class="timestamp-wrapper"><span class="timestamp">&lt;2018-12-04 Tue&gt;</span></span>

-   Added

    -   Dosing rate as parameter

-   Changed

    -   Update to Stan version 2.18.0.


### Version 0.84 <span class="timestamp-wrapper"><span class="timestamp">&lt;2018-02-24 Sat&gt;</span></span>

-   Added

    -   Piecewise linear interpolation function.
    -   Univariate integral functions.

-   Changed

    -   Update to Stan version 2.17.1.
    -   Minor revisions to User Manual.
    -   Bugfixes.


### Version 0.83 <span class="timestamp-wrapper"><span class="timestamp">&lt;2017-08-02 Wed&gt;</span></span>

-   Added

    -   Work with TorstenHeaders
    -   Each chain has a different initial estimate

-   Changed

    -   User manual
    -   Fix misspecification in ODE system for TwoCpt example.
    -   Other bugfixes


### Version 0.82 <span class="timestamp-wrapper"><span class="timestamp">&lt;2017-01-29 Sun&gt;</span></span>

-   Added

    -   Allow parameter arguments to be passed as 1D or 2D arrays
    -   More unit tests
    -   Unit tests check automatic differentiation against finite differentiation.

-   Changed

    -   Split the parameter argument into three arguments: pMatrix (parameters for the ODEs &#x2013; note: for `linOdeModel`, pMatrix is replaced by the constant rate matrix K), biovar (parameters for the biovariability), and tlag (parameters for the lag time).
    -   bugfixes


### Version 0.81 <span class="timestamp-wrapper"><span class="timestamp">&lt;2016-09-27 Tue&gt;</span></span>

-   Added

    linCptModel (linear compartmental model) function


### Version 0.80a <span class="timestamp-wrapper"><span class="timestamp">&lt;2016-09-21 Wed&gt;</span></span>

-   Added

    check_finite statements in pred_1 and pred_2 to reject metropolis proposal if initial conditions are not finite
    
    \bibliography{./docs/torsten}

## Footnotes

<sup><a id="fn.1" class="footnum" href="#fnr.1">1</a></sup> NONMEM\textregistered{} is licensed and distributed by ICON Development Solutions.
