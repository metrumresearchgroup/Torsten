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
-   Single dose and multiple dose events
-   Handles steady state dosing events
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

\bibliography{./docs/torsten}

## Footnotes

<sup><a id="fn.1" class="footnum" href="#fnr.1">1</a></sup> NONMEM\textregistered{} is licensed and distributed by ICON Development Solutions.
