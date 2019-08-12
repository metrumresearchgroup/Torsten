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

In general, all real variables may be passed as Stan parameters. A few exceptions apply /to functions which use a numerical integrator/(i.e. the general and the mix compartment models). The below listed cases present technical difficulties, which we expect to overcome in Torsten's next release:

-   In the case of a multiple truncated infusion rate dosing regimen:
    -   The bioavailability (F) and the amount (AMT) must be fixed.

This library provides Stan language functions that calculate amounts in each compartment, given an event schedule and an ODE system.


## Installation

We are working with Stan development team to create a system to add and share Stan packages. In the mean time, the current repo contains forked version of Stan with Torsten. The latest version of Torsten (v0.87) is compatible with Stan v2.19.1. Torsten is agnostic to which Stan interface you use. Here we provide command line and R interfaces.


### Command line interface

After downloading the project

-   <https://github.com/metrumresearchgroup/Torsten>

The command line interface `cmdstan` is available to use without installation. The following command builds a Torsten model `model_name` in `model_path`

```sh
cd $TORSTEN_PATH/cmdstan; make model_path/model_name
```

Currently MPI support is only available through `cmdstan` interface. To use MPI-supported population/group solvers, add/edit `make/local`

```sh
TORSTEN_MPI=1

# path to MPI headers
CXXFLAGS += -isystem /usr/local/include
# if you are using Metrum's metworx platform, add MPICH3's
# headers with
# CXXFLAGS += -isystem /usr/local/mpich3/include
```

Note that currently `TORSTEN_MPI` and `STAN_MPI` flags conflict on processes management and cannot be used in a same Stan model.


### R interface

The R interface is based on [rstan](https://cran.r-project.org/web/packages/rstan/index.html), the Stan's interface for R. To install R version of Torsten, at `$TORSTEN_PATH`, in R

```R
source('install.R')
```

Please ensure the R toolchain includes a C++ compiler with C++14 support. In particular, R 3.4.0 and later is recommended as it contains toolchain based on gcc 4.9.3. On Windows platform, such a toolchain can be found in Rtools34 and later.

Please ensure `.R/Makevars` constains the following flags

```sh
CXX14 = g++ -fPIC               # or CXX14 = clang++ -fPIC

CXXFLAGS=-O3 -std=c++1y -mtune=native -march=native -Wno-unused-variable -Wno-unused-function
CXXFLAGS += -DBOOST_MPL_CFG_NO_PREPROCESSED_HEADERS -DBOOST_MPL_LIMIT_LIST_SIZE=30

CXX14FLAGS=-O3 -std=c++1y -mtune=native -march=native -Wno-unused-variable -Wno-unused-function
CXX14FLAGS += -DBOOST_MPL_CFG_NO_PREPROCESSED_HEADERS -DBOOST_MPL_LIMIT_LIST_SIZE=30
```

Fore more information of setting up `makevar` and its functionality, see

-   <http://dirk.eddelbuettel.com/code/rcpp/Rcpp-package.pdf>

For more information of installation troubleshooting, please consult [rstan wiki](https://github.com/stan-dev/rstan/wiki).


### Testing

With project in `torsten_path`, set the envionment variable `TORSTEN_PATH` as

```sh
# in bash
export TORSTEN_PATH=torsten_path
# in csh
setenv TORSTEN_PATH torsten_path
```

To test the installation, run

```sh
./test-torsten.sh --unit        # math unit test
./test-torsten.sh --signature   # stan function # signature test
./test-torsten.sh --model       # R model test, takes long time to finish
```


## Development plans

Our current plans for future development of Torsten include the following:

-   Build a system to easily share packages of Stan functions (written in C++ or in the Stan language)
-   Allow numerical methods to handle bioavailability fraction (F) as parameters in all cases.
-   Optimize Matrix exponential functions
    -   Function for the action of Matrix Exponential on a vector
    -   Hand-coded gradients
    -   Special algorithm for matrices with special properties
-   Fix issue that arises when computing the adjoint of the lag time parameter (in a dosing compartment) evaluated at \(t_{\text{lag}} = 0\).
-   Extend formal tests
    -   We want more C++ Google unit tests to address cases users may encounter
    -   Comparison with simulations from the R package `mrgsolve` and the software NONMEM\textregistered{}
    -   Recruit non-developer users to conduct beta testing


## Changelog


### 0.87 <span class="timestamp-wrapper"><span class="timestamp">&lt;2019-07-26 Fri&gt;</span></span>

-   Added

    -   MPI dynamic load balance for Torsten's population ODE integrators
        
        -   `pmx_integrate_ode_group_adams`
        -   `pmx_integrate_ode_group_bdf`
        -   `pmx_integrate_ode_group_rk45`
        
        To invoke dynamic load balance instead of default static balance for MPI, issue `TORSTEN_MPI=2` in `make/local`.
    -   Support `time` as parameter in `pmx_solve_rk45/bdf/adams` functions.

-   Changed

    -   Some fixes on steady-state solvers
    -   Update to rstan version 2.19.2.


### 0.86 <span class="timestamp-wrapper"><span class="timestamp">&lt;2019-05-15 Wed&gt;</span></span>

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


### 0.85 <span class="timestamp-wrapper"><span class="timestamp">&lt;2018-12-04 Tue&gt;</span></span>

-   Added

    -   Dosing rate as parameter

-   Changed

    -   Update to Stan version 2.18.0.


### 0.84 <span class="timestamp-wrapper"><span class="timestamp">&lt;2018-02-24 Sat&gt;</span></span>

-   Added

    -   Piecewise linear interpolation function.
    -   Univariate integral functions.

-   Changed

    -   Update to Stan version 2.17.1.
    -   Minor revisions to User Manual.
    -   Bugfixes.


### 0.83 <span class="timestamp-wrapper"><span class="timestamp">&lt;2017-08-02 Wed&gt;</span></span>

-   Added

    -   Work with TorstenHeaders
    -   Each chain has a different initial estimate

-   Changed

    -   User manual
    -   Fix misspecification in ODE system for TwoCpt example.
    -   Other bugfixes


### 0.82 <span class="timestamp-wrapper"><span class="timestamp">&lt;2017-01-29 Sun&gt;</span></span>

-   Added

    -   Allow parameter arguments to be passed as 1D or 2D arrays
    -   More unit tests
    -   Unit tests check automatic differentiation against finite differentiation.

-   Changed

    -   Split the parameter argument into three arguments: pMatrix (parameters for the ODEs &#x2013; note: for `linOdeModel`, pMatrix is replaced by the constant rate matrix K), biovar (parameters for the biovariability), and tlag (parameters for the lag time).
    -   bugfixes


### 0.81 <span class="timestamp-wrapper"><span class="timestamp">&lt;2016-09-27 Tue&gt;</span></span>

-   Added

    linCptModel (linear compartmental model) function


### 0.80a <span class="timestamp-wrapper"><span class="timestamp">&lt;2016-09-21 Wed&gt;</span></span>

-   Added

    check_finite statements in pred_1 and pred_2 to reject metropolis proposal if initial conditions are not finite
    
    \bibliography{./doc/torsten}

## Footnotes

<sup><a id="fn.1" class="footnum" href="#fnr.1">1</a></sup> NONMEM\textregistered{} is licensed and distributed by ICON Development Solutions.
