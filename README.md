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

-   The RATE and TIME arguments must be fixed
-   In the case of a multiple truncated infusion rate dosing regimen:
    -   The bioavailability (F) and the amount (AMT) must be fixed.

This library provides Stan language functions that calculate amounts in each compartment, given an event schedule and an ODE system.


## Installation

Installation files are available on GitHub

-   <https://github.com/metrumresearchgroup/Torsten>

We are working with Stan development team to create a system to add and share Stan packages. In the mean time, the current repo contains forked version of Stan with Torsten. The latest version of Torsten (v0.84) is compatible with Stan v2.17.1. Torsten is agnostic to which Stan interface you use. Here we provide command line and R interfaces.

First, download the project. The root path of the project is your `torsten_path`. Set the envionment variable `TORSTEN_PATH` to `torsten_path`. In bash this is

```sh
export TORSTEN_PATH=torsten_path
```


### Command line version

There is no need to install command line version. To build a Stan model with `model_name` in `model_path` using command line, in `torsten_path/cmdstan`, do

```sh
make model_path/model_name
```


### R version

The R version is based on [rstan](https://cran.r-project.org/web/packages/rstan/index.html), the Stan's interface for R. To install R version of Torsten, at `torsten_path`, in R

```R
source('install.R')
```

Please ensure the R toolchain includes a C++ compiler with C++11 support. In particular, R 3.3.0 and later is recommended as it contains toolchain based on gcc 4.9.3. On Windows platform, such a toolchain can be found in Rtools33 and later.

Please ensure CXXFLAGS in .R/Makevars constains flag -std=c++11.

For more information of installation troubleshooting, please consult [rstan wiki](https://github.com/stan-dev/rstan/wiki).


### Testing

To test after installation, run

```sh
./test-torsten.sh --unit --signature --model
```


## Development plans

Our current plans for future development of Torsten include the following:

-   Build a system to easily share packages of Stan functions (written in C++ or in the Stan language)
-   Allow numerical methods to handle RATE, AMT, TIME, and the bioavailability fraction (F) as parameters in all cases.
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
