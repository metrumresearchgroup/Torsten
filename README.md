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

We are working with Stan development team to create a system to add and share Stan packages. In the mean time, the current repo contains forked version of Stan with Torsten. The latest version of Torsten (v0.85) is compatible with Stan v2.18.1. Torsten is agnostic to which Stan interface you use. Here we provide command line and R interfaces.

After downloading the project

-   <https://github.com/metrumresearchgroup/Torsten>

to `torsten_path`, set the envionment variable `TORSTEN_PATH` as

```sh
# in bash
export TORSTEN_PATH=torsten_path
# in csh
setenv TORSTEN_PATH torsten_path
```

-   Command line interface

    The command line interface `cmdstan` does not require installation. The following command builds a Torsten model `model_name` in `model_path`
    
    ```sh
    cd $TORSTEN_PATH/cmdstan; make model_path/model_name
    ```

-   R interface

    The R interface is based on [rstan](https://cran.r-project.org/web/packages/rstan/index.html), the Stan's interface for R. To install R version of Torsten, at `$TORSTEN_PATH`, in R
    
    ```R
    source('install.R')
    ```
    
    Please ensure the R toolchain includes a C++ compiler with C++14 support. In particular, R 3.4.0 and later is recommended as it contains toolchain based on gcc 4.9.3. On Windows platform, such a toolchain can be found in Rtools34 and later.
    
    Please ensure `.R/Makevars` constains the following flags
    
    ```sh
    CXXFLAGS=-O3 -std=c++1y -mtune=native -march=native -Wno-unused-variable -Wno-unused-function
    CXXFLAGS += -DBOOST_MPL_CFG_NO_PREPROCESSED_HEADERS -DBOOST_MPL_LIMIT_LIST_SIZE=30
    
    CXX14FLAGS=-O3 -std=c++1y -mtune=native -march=native -Wno-unused-variable -Wno-unused-function
    CXX14FLAGS += -DBOOST_MPL_CFG_NO_PREPROCESSED_HEADERS -DBOOST_MPL_LIMIT_LIST_SIZE=30
    ```
    
    For more information of installation troubleshooting, please consult [rstan wiki](https://github.com/stan-dev/rstan/wiki).

-   Testing

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
    -   Comparison with simulations from the R package *mrgsolve* and the software NONMEM\textregistered{}
    -   Recruit non-developer users to conduct beta testing


## Changelog


### 0.85 <span class="timestamp-wrapper"><span class="timestamp">&lt;2018-10-20 Sat&gt;</span></span>

-   Added

    -   Dosing rate as parameter

-   Changed

    -   Update with Stan version 2.18.0.


### 0.84 <span class="timestamp-wrapper"><span class="timestamp">&lt;2018-02-24 Sat&gt;</span></span>

-   Added

    -   Piecewise linear interpolation function.
    -   Univariate integral functions.

-   Changed

    -   Update with Stan version 2.17.1.
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

    -   Split the parameter argument into three arguments: pMatrix (parameters for the ODEs &#x2013; note: for linOdeModel, pMatrix is replaced by the constant rate matrix K), biovar (parameters for the biovariability), and tlag (parameters for the lag time).
    -   bugfixes


### 0.81 <span class="timestamp-wrapper"><span class="timestamp">&lt;2016-09-27 Tue&gt;</span></span>

-   Added

    linCptModel (linear compartmental model) function


### 0.80a <span class="timestamp-wrapper"><span class="timestamp">&lt;2016-09-21 Wed&gt;</span></span>

-   Added

    check<sub>finite</sub> statements in pred<sub>1</sub> and pred<sub>2</sub> to reject metropolis proposal if initial conditions are not finite
    
    \bibliography{./doc/torsten}

## Footnotes

<sup><a id="fn.1" class="footnum" href="#fnr.1">1</a></sup> NONMEM\textregistered{} is licensed and distributed by ICON Development Solutions.
