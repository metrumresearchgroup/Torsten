+++
title = "Installation"
author = ["Yi Zhang"]
date = 2021-06-25T00:00:00-07:00
lastmod = 2021-06-29T11:51:41-07:00
draft = false
weight = 1005
[menu.main]
  weight = 1005
  identifier = "installation"
+++

Currently Torsten is based on a forked version of Stan and hosted on GitHub

-   <https://github.com/metrumresearchgroup/Torsten>

The latest v0.89rc is
compatible with Stan v2.27.0. Torsten can be accessed from
command line for cmdstan interface and `cmdstanr`
(<https://mc-stan.org/cmdstanr/>) for R interface. It requires
a modern C++11 compiler as well as a Make utility. See ([Team 2020](#org22b1af0)) for details of installation and
required toolchain. In particular, we recommend the folowing versions
of C++ compilers:

-   Linux: g++ >=7.5 or clang >=8.0,
-   macOS: the XCode version of clang,
-   Windows: g++ 8.1 (available with RTools 4.0).

On windows, the Make utility `mingw32-make` can be installed as part
of RTools.


### Command line interface {#command-line-interface}

The command line interface `cmdstan` is available along with Torsten
and can be found at `Torsten/cmdstan`.

After installation, one can use the following command to build a Torsten model `model_name` in `model_path`

```sh
cd Torsten/cmdstan
make model_path/model_name # replace "make" with "mingw32-make" on Windows platform
```


### R interface {#r-interface}

After installing cmdstanr from <https://mc-stan.org/cmdstanr/>, use the
following command to set path

```r
cmdstanr::set_cmdstan_path("Torsten/cmdstan")
```

Then one can follow <https://mc-stan.org/cmdstanr/articles/cmdstanr.html> to compile
and run Torsten models.


### MPI support {#mpi-support}

Torsten's MPI support is of a different flavour than
`reduce_sum` found in Stan. To be able to utilize MPI
parallelisation, one first needs to ensure an MPI library
such as

-   <https://www.mpich.org/downloads/>
-   <https://www.open-mpi.org/software/ompi/>

is available. Torsen's implementation is tested on
both `MPICH` and `OpenMPI`.

To use MPI-supported population/group solvers,
add/edit `make/local`

```sh
TORSTEN_MPI=1

# path to MPI headers
CXXFLAGS += -isystem /usr/local/include
# if you are using Metrum's metworx platform, add MPICH3's
# headers with
# CXXFLAGS += -isystem /usr/local/mpich3/include
```

Note that currently `TORSTEN_MPI` and `STAN_MPI` flags
conflict on processes management and cannot be used in a
same Stan model, and MPI support is only available through `cmdstan`
interface.


### Testing {#testing}

Models in `example-models` directory are for tutorial and demonstration.
The following shows how to build and run the two-compartment model
using `cmdstanr`, and use `bayesplot` to examine posterior density of `CL`.

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


## Bibliography {#bibliography}

<a id="org22b1af0"></a>Team, Stan Development. 2020. _CmdStan User’s Guide_. <https://mc-stan.org/docs/2%5F26/cmdstan-guide/index.html>.
