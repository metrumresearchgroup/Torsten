+++
title = "Torsten: A Pharmacokinetic/Pharmacodynamic Model Library for Stan"
author = ["Yi Zhang"]
draft = false
+++

### Development team {#development-team}


#### Bill Gillespie {#bill-gillespie}

[`billg@metrumrg.com`](mailto:billg@metrumrg.com),
Metrum Research Group, LLC


#### Yi Zhang {#yi-zhang}

[`yiz@metrumrg.com`](mailto:yiz@metrumrg.com), Metrum Research Group, LLC


#### Charles Margossian {#charles-margossian}

[`charles.margossian@columbia.edu`](mailto:charles.margossian@columbia.edu), Columbia University, Department of Statistics(formerly Metrum Research Group, LLC)


## Acknowledgements {#acknowledgements}


### Institutions {#institutions}

We thank Metrum Research Group, Columbia University, and AstraZeneca.


### Funding {#funding}

This work was funded in part by the following organizations:


#### Office of Naval Research (ONR) contract N00014-16-P-2039 {#office-of-naval-research--onr--contract-n00014-16-p-2039}

provided as part of the Small Business Technology Transfer (STTR)
program. The content of the information presented in this document
does not necessarily reflect the position or policy of the
Government and no official endorsement should be inferred.


#### Bill & Melinda Gates Foundation. {#bill-and-melinda-gates-foundation-dot}


### Individuals {#individuals}

We thank the Stan Development Team for giving us guidance on how to
create new Stan functions and adding features to Stan's core language
that facilitate building ODE-based models.

We also thank Kyle Baron and Hunter Ford for helpful advice on coding
in C++ and using GitHub, Curtis Johnston for reviewing the User
Manual, and Yaming Su for using Torsten and giving us feedback.


## Introduction {#introduction}

Stan is an open source probabilistic programing language designed
primarily to do Bayesian data analysis
<sup id="db2e464c2efb7b0f512669a129827ca4"><a href="#carpenter17_stan" title="Carpenter, Gelman, , Hoffman, Lee, , Goodrich, Betancourt, , Brubaker, Guo, Li, Peter \&amp; Riddell, Stan: {A} {Probabilistic} {Programming}  {Language}, {Journal of Statistical software}, v(), (2017).">carpenter17_stan</a></sup>. Several of its
features make it a powerful tool to specify and fit complex
models. First, its language is very expressive and
flexible. Secondly, it implements a variant of No U-Turn
Sampler(NUTS), an adaptative Hamiltonian Monte Carlo
algorithm that was proven more efficient than commonly used Monte Carlo Markov Chains
(MCMC) samplers for complex high dimensional problems <sup id="b10c573e0d702948c57bd33da6e3717e"><a href="#hoffman_no-u-turn_2011" title="Hoffman \&amp; Gelman, The {No}-{U}-{Turn} {Sampler}: {Adaptively} {Setting} {Path} {Lengths} in {Hamiltonian} {Monte} {Carlo}, {arXiv:1111.4246 [cs, stat]}, v(), (2011).">hoffman_no-u-turn_2011</a></sup><sup>,</sup><sup id="fe3285d2c36f9e6ca95ea6664e70a82f"><a href="#betancourt_hmc_2018" title="Betancourt, A Conceptual Introduction to Hamiltonian Monte Carlo, v(), (2018).">betancourt_hmc_2018</a></sup>. Our
goal is to harness these innovations and make Stan a better
software for pharmacometrics modeling. Our efforts are twofold:

1.  We contribute to the development of features, such as functions that support differential equations based models, and implement them directly into Stan's core language.
2.  We develop Torsten, an extension with specialized pharmacometrics functions.

Throughout the process, we work very closely with the Stan Development
Team. We have benefited immensely from their mentorship, advice, and
feedback. Just like Stan, Torsten is an open source project that
fosters collaborative work. Interested in contributing?
Comment at Torsten repository

[https://github.com/metrumresearchgroup/Torsten](https://github.com/metrumresearchgroup/Torsten)

or email us ([billg@metrumrg.com](mailto:billg@metrumrg.com), [yiz@metrumrg.com](mailto:yz@yizh.org)) and we will help you help us!

Torsten is licensed under the BSD 3-clause license.

<div class="mdframed">
  <div></div>

**WARNING:** The current version of Torsten is a _prototype_. It
is being released for review and comment, and to support limited
research applications. It has not been rigorously tested and should
not be used for critical applications without further testing or
cross-checking by comparison with other methods.

We encourage interested users to try Torsten out and are happy to
assist. Please report issues, bugs, and feature requests on
[our GitHub page](https://github.com/metrumresearchgroup/stan).

</div>


### Overview {#overview}

Torsten is a collection of Stan functions to facilitate analysis of
pharmacometric data using Stan. The current version
includes:

-   Specific linear compartment models:
    -   One compartment model with first order absorption.
    -   Two compartment model with elimination from and first order absorption into central compartment
-   General linear compartment model described by a system of first-order \underline{linear} Ordinary Differential Equations (ODEs).
-   General compartment model described by a system of first order ODEs.
-   Mix compartment model with PK forcing function described by a linear one or two compartment model.

The models and data format are based on
NONMEM\textregistered{}[^fn:1]/NMTRAN/PREDPP
conventions including:

-   Recursive calculation of model predictions
    -   This permits piecewise constant covariate values
-   Bolus or constant rate inputs into any compartment
-   Single dose and multiple dose events
-   Handles steady state dosing events
-   Implemented NMTRAN data items include: TIME, EVID, CMT, AMT, RATE, ADDL, II, SS

In general, all real variables may be passed as model parameters. A
few exceptions apply /to functions which use a numerical
integrator(i.e. the general and the mix compartment
models). The below listed cases present technical difficulties, which we expect to
overcome in Torsten's next release:

-   In the case of a multiple truncated infusion rate dosing regimen:
    -   The bioavailability (F) and the amount (AMT) must be fixed.

This library provides Stan language functions that calculate amounts
in each compartment, given an event schedule and an ODE system.


### Implementation summary {#imp_details}

-   Current Torsten v0.89rc is based on Stan v2.27.0.
-   All functions are programmed in C++ and are compatible
    with the Stan math automatic differentiation library <sup id="5edcf88ae76ae7132583b12fffe35055"><a href="#carpenter15_stan_math_librar" title="Carpenter, Hoffman, , Brubaker, Lee, Li, Peter \&amp; Betancourt, The {Stan} {Math} {Library}:  {Reverse}-{Mode} {Automatic}  {Differentiation} in {C}++, {arXiv:1509.07164 [cs]}, v(), (2015).">carpenter15_stan_math_librar</a></sup>
-   One and two compartment models are based on analytical solutions of governing ODEs.
-   General linear compartment models are based on semi-analytical solutions using the built-in matrix exponential function
-   General compartment models are solved numerically using built-in ODE integrators in Stan. The tuning parameters of the solver are adjustable. The steady state solution is calculated using a numerical algebraic solver.
-   A mix compartment model's PK forcing function is solved analytically, and its forced ODE system is solved numerically.


### Development plans {#dev_plans}

Our current plans for future development of Torsten include the
following:

-   Build a system to easily share packages of Stan functions
    (written in C++ or in the Stan language)
-   Optimize Matrix exponential functions
    -   Function for the action of Matrix Exponential on a vector
    -   Hand-coded gradients
    -   Special algorithm for matrices with special properties
-   Develop new method for large-scale hierarchical models with costly
    ODE solving.
-   Fix issue that arises when computing the adjoint of the lag time
    parameter (in a dosing compartment) evaluated at \\(t\_{\text{lag}} = 0\\).
-   Extend formal tests
    -   More unit tests and better CD/CI support.
    -   Comparison with simulations from the R package
        `mrgsolve` and the software NONMEM\textregistered{}
    -   Recruit non-developer users to conduct beta testing


### Changelog {#changelog}


#### Version 0.89 <span class="timestamp-wrapper"><span class="timestamp">&lt;2021-06-15 Tue&gt;</span></span> {#version-0-dot-89}

<!--list-separator-->

-  Changed

    -   New backend for ODE events solvers.
    -   Use vector instead of array as ODE function state & return type.
    -   Simplified ODE integrator naming,
        e.g.  and  .
    -   Update to Stan version 2.27.0.


#### Version 0.88 <span class="timestamp-wrapper"><span class="timestamp">&lt;2020-12-18 Fri&gt;</span></span> {#version-0-dot-88}

<!--list-separator-->

-  Added

    -   Bioavailability, lag time, ODE real & integer data are optional in PMX function signatures.
    -   Support all EVID options from NM-TRAN and mrgsolve.
    -   Support steady-state infusion through multiple interdose intervals.

<!--list-separator-->

-  Changed

    -   More efficient memory management of COVDES implenmentation.
    -   Update of MPI framework to adapt multilevel paralleism.
    -   Update to Stan version 2.25.0.
    -   Use cmdstanr as R interface.
    -   Stop supporting rstan as R interface.


#### Version 0.87 <span class="timestamp-wrapper"><span class="timestamp">&lt;2019-07-26 Fri&gt;</span></span> {#version-0-dot-87}

<!--list-separator-->

-  Added

    -   MPI dynamic load balance for Torsten's population ODE integrators

        -   `pmx_integrate_ode_group_adams`
        -   `pmx_integrate_ode_group_bdf`
        -   `pmx_integrate_ode_group_rk45`

        To invoke dynamic load balance instead of default static
        balance for MPI, issue `TORSTEN_MPI=2` in `make/local`.
    -   Support `RATE` as parameter in `pmx_solve_rk45/bdf/adams`
        functions.

<!--list-separator-->

-  Changed

    -   Some fixes on steady-state solvers
    -   Update to rstan version 2.19.2.


#### Version 0.86 <span class="timestamp-wrapper"><span class="timestamp">&lt;2019-05-15 Wed&gt;</span></span> {#version-0-dot-86}

<!--list-separator-->

-  Added

    -   Torsten's ODE integrator functions

        -   `pmx_integrate_ode_adams`
        -   `pmx_integrate_ode_bdf`
        -   `pmx_integrate_ode_rk45`

        and their counterparts to solve a population/group of
        subjects governed by an ODE

        -   `pmx_integrate_ode_group_adams`
        -   `pmx_integrate_ode_group_bdf`
        -   `pmx_integrate_ode_group_rk45`
    -   Torsten's population PMX solver functions for general
        ODE models
        -   `pmx_solve_group_adams`
        -   `pmx_solve_group_bdf`
        -   `pmx_solve_group_rk45`
    -   Support time step `ts` as parameter in `pmx_integrate_ode_xxx`
        solvers.

<!--list-separator-->

-  Changed

    -   Renaming Torsten functions in previous releases, the
        old-new name mapping is

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

        Note that the new version of the above functions return
        the _transpose_ of the matrix returned by the old
        versions, in order to improve memory efficiency. The old version are retained but will be
        deprecated in the future.
    -   Update to Stan version 2.19.1.


#### Version 0.85 <span class="timestamp-wrapper"><span class="timestamp">&lt;2018-12-04 Tue&gt;</span></span> {#version-0-dot-85}

<!--list-separator-->

-  Added

    -   Dosing rate as parameter

<!--list-separator-->

-  Changed

    -   Update to Stan version 2.18.0.


#### Version 0.84 <span class="timestamp-wrapper"><span class="timestamp">&lt;2018-02-24 Sat&gt;</span></span> {#version-0-dot-84}

<!--list-separator-->

-  Added

    -   Piecewise linear interpolation function.
    -   Univariate integral functions.

<!--list-separator-->

-  Changed

    -   Update to Stan version 2.17.1.
    -   Minor revisions to User Manual.
    -   Bugfixes.


#### Version 0.83 <span class="timestamp-wrapper"><span class="timestamp">&lt;2017-08-02 Wed&gt;</span></span> {#version-0-dot-83}

<!--list-separator-->

-  Added

    -   Work with TorstenHeaders
    -   Each chain has a different initial estimate

<!--list-separator-->

-  Changed

    -   User manual
    -   Fix misspecification in ODE system for TwoCpt example.
    -   Other bugfixes


#### Version 0.82 <span class="timestamp-wrapper"><span class="timestamp">&lt;2017-01-29 Sun&gt;</span></span> {#0-82-added}

<!--list-separator-->

-  Added

    -   Allow parameter arguments to be passed as 1D or 2D arrays
    -   More unit tests
    -   Unit tests check automatic differentiation against finite differentiation.

<!--list-separator-->

-  Changed

    -   Split the parameter argument into three arguments: pMatrix
        (parameters for the ODEs -- note: for `linOdeModel`, pMatrix
        is replaced by the constant rate matrix K), biovar
        (parameters for the biovariability), and tlag (parameters
        for the lag time).
    -   bugfixes


#### Version 0.81 <span class="timestamp-wrapper"><span class="timestamp">&lt;2016-09-27 Tue&gt;</span></span> {#version-0-dot-81}

<!--list-separator-->

-  Added

    linCptModel (linear compartmental model) function


#### Version 0.80a <span class="timestamp-wrapper"><span class="timestamp">&lt;2016-09-21 Wed&gt;</span></span> {#version-0-dot-80a}

<!--list-separator-->

-  Added

    check_finite statements in pred_1 and pred_2 to reject metropolis proposal if initial conditions are not finite


## Installation {#installation}

Currently Torsten is based on a forked version of Stan and hosted on GitHub

-   <https://github.com/metrumresearchgroup/Torsten>

The latest v0.89rc is
compatible with Stan v2.27.0. Torsten can be accessed from
command line for cmdstan interface and `cmdstanr`
(<https://mc-stan.org/cmdstanr/>) for R interface. It requires
a modern C++11 compiler as well as a Make utility. See <sup id="8aa45ba469168c91a8faa4436edebb97"><a href="#cmdstan_team_2020" title="@manual{cmdstan_team_2020,
    Author = {Stan Development Team},
    Title = {CmdStan User's Guide},
    Year = {2020},
    url = {https://mc-stan.org/docs/2_26/cmdstan-guide/index.html},
    }">cmdstan_team_2020</a></sup> for details of installation and
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

\label{sec:mpi\_support}
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


## Using Torsten {#using-torsten}

<a id="org1fe3fd5"></a>

The reader should have a basic understanding of how Stan works before
reading this chapter. There are excellent resources online to get
started with Stan ([http://mc-stan.org/documentation](http://mc-stan.org/documentation)).
In this section we go through the different functions Torsten adds to
Stan. The code for the examples can be found at the `example-models` folder.

Torsten's functions are prefixed with .
For some of their arguments we adopt NM-TRAN format for events
specification(Table [tab:event_args](#tab:event_args)).

<a id="table--tab:event-args"></a>
<div class="table-caption">
  <span class="table-number"><a href="#table--tab:event-args">Table 1</a></span>:
  NM-TRAN compatible event specification arguments. All arrays should have the same length corresponding to the number of events.
</div>

| Argument Name | Definition                  | Stan data type |
|---------------|-----------------------------|----------------|
|               | time                        |                |
|               | amount                      |                |
|               | infusion rate               |                |
|               | interdose interval          |                |
|               | event ID                    |                |
|               | event compartment           |                |
|               | additionial identical doses |                |
|               | steady-state dosing flag    |                |

All the `real[]` arguments above are allowed to
be `parameters` in a Stan model.
In addtion, Torsten functions
support optional arguments and overloaded signatures.
Optional arguments are indicated by surrounding square bracket `[]`.
Table below shows three commonly used PMX model arguments that support
overloading. In the rest of this document we assume this convention unless indicated otherwise.

<a id="table--tab:event-params"></a>
<div class="table-caption">
  <span class="table-number"><a href="#table--tab:event-params">Table 2</a></span>:
  PMX model parameter overloadings. One can use 1d array <code class="src src-stan"><span style="color: #b58900;">real</span>[]</code> to indicate constants of all events, or 2d array <code class="src src-stan"><span style="color: #b58900;">real</span>[ , ]</code> so that the \(i\)th row of the array describes the model arguments for time interval \((t_{i-1}, t_i)\), and the number of the rows equals to the size of .
</div>

| Argument Name | Definition               | Stan data type | Optional           |
|---------------|--------------------------|----------------|--------------------|
|               | model parameters         |  or            | N                  |
|               | bioavailability fraction |  or            | Y (default to 1.0) |
|               | lag time                 |  or            | Y (default to 0.0) |


### One Compartment Model {#one-compartment-model}

\label{sec:onecpt}


#### Description {#description}

Function  solves a one-compartment PK
model (Figure [1](#org3446f9f)). The model obtains plasma concentrations of parent drug \\(c=y\_2/V\_2\\)
by solving for the mass of drug in the central compartment
\\(y\_2\\) from ordinary differential equations(ODEs)

\begin{subequations}
\begin{align}
  y\_1' &= -k\_a y\_1, \\\\\\
  y\_2' &= k\_a y\_1 - \left(\frac{CL}{V\_2} + \frac{Q}{V\_2}\right) y\_2.
\end{align}
\label{eq:onecpt}
\end{subequations}

<a id="org3446f9f"></a>

{{< figure src="/ox-hugo/cptModels.png" caption="Figure 1: One and two compartment models with first order absorption implemented in Torsten." >}}


#### Usage {#usage}

```stan
matrix = pmx_solve_onecpt(time, amt, rate, ii, evid, cmt, addl, ss, theta [, biovar, tlag ] )
```


#### Arguments {#arguments}

See Table [tab:event_args](#tab:event_args) and Table [tab:event_params](#tab:event_params).


#### Return value {#return-value}

An `ncmt`-by-`nt` matrix, where `nt` is the number of time steps and `ncmt=2` is the number of compartments.


#### Note {#note}

-   ODE Parameters  should consist of \\(CL\\), \\(V\_2\\), \\(k\_a\\), in that order.
-   and  are optional, so that the following are allowed:

<!--listend-->

```stan
pmx_solve_onecpt(..., theta);
pmx_solve_onecpt(..., theta, biovar);
pmx_solve_onecpt(..., theta, biovar, tlag);
```

-   Setting \\(k\_a = 0\\) eliminates the first-order absorption.


### Two Compartment Model {#two-compartment-model}

<a id="orga5326f3"></a>


#### Description {#description}

Function  solves a two-compartment PK
model (Figure [1](#org3446f9f)). The model obtains plasma concentrations of parent drug \\(c=y\_2/V\_2\\)
by solving for the mass of drug in the central compartment
\\(y\_2\\) from ordinary differential equations(ODEs)

\begin{subequations}
  \begin{align} \label{eq:twocpt}
    y\_1' &= -k\_a y\_1 \\\\\\
    y\_2' &= k\_a y\_1 - \left(\frac{CL}{V\_2} + \frac{Q}{V\_2}\right) y\_2 +  \frac{Q}{V\_3}  y\_3  \\\\\\
    y\_3' &= \frac{Q}{V\_2} y\_2 - \frac{Q}{V\_3} y\_3
  \end{align}
\end{subequations}


#### Usage {#usage}

```stan
matrix = pmx_solve_twocpt(time, amt, rate, ii, evid, cmt, addl, ss, theta [, biovar, tlag ] )
```


#### Arguments {#arguments}

See Table [tab:event_args](#tab:event_args) and Table [tab:event_params](#tab:event_params).


#### Return value {#return-value}

An `ncmt`-by-`nt` matrix, where `nt` is the number of time steps and `ncmt=3` is the number of compartments.


#### Note {#note}

-   ODE Parameters  consists of \\(CL\\), \\(Q\\), \\(V\_2\\), \\(V\_3\\), \\(k\_a\\).
-   and  are optional, so that the following are allowed:

<!--listend-->

```stan
pmx_solve_twocpt(..., theta);
pmx_solve_twocpt(..., theta, biovar);
pmx_solve_twocpt(..., theta, biovar, tlag);
```

-   Setting \\(k\_a = 0\\) eliminates the first-order absorption.


### General Linear ODE Model Function {#general-linear-ode-model-function}


#### Description {#description}

Function  solves a (piecewise) linear ODEs model with coefficients
in form of matrix \\(K\\)

\begin{equation}
y^\prime\left(t\right) = Ky\left(t\right)
\end{equation}

For example, in a two-compartment model with first order absorption, \\(K\\) is

\begin{equation}
  K = \left[\begin{array}{ccc}
              -k\_a & 0 & 0 \\\\\\
              k\_a & -\left(k\_{10} + k\_{12}\right) & k\_{21} \\\\\\
              0 & k\_{12} & -k\_{21}
            \end{array}\right]
\end{equation}

where \\(k\_{10}=CL/V\_2\\), \\(k\_{12}=Q/V\_2\\), and \\(k\_{21}=Q/V\_3\\).


#### Usage {#usage}

```stan
matrix = pmx_solve_linode(time, amt, rate, ii, evid, cmt, addl, ss, K, biovar, tlag )
```


#### Arguments {#arguments}

<!--list-separator-->

-

    System parameters.  can be either

    -   a `matrix` for constant parameters in all events, or
    -   an array of matrices `matrix K[nt]` so that the \\(i\\)th entry of the array describes
        the model parameters for time interval \\((t\_{i-1}, t\_i)\\),
        and the number of the rows equals to the number of event time .

<!--list-separator-->

-  See Table [tab:event_args](#tab:event_args) and Table [tab:event_params](#tab:event_params) for the rest of arguments.


#### Return value {#return-value}

An `n`-by-`nt` matrix, where `nt` is the number of time steps and `n` is the number of rows(columns) of square matrix .


### General ODE Model Function {#general-ode-model-function}

\label{sec:general\_ode}


#### Description {#description}

Function , , and  solve a first-order ODE system
specified by user-specified right-hand-side function  \\(f\\)

\begin{equation\*}
y'(t) = f(t, y(t))
\end{equation\*}

In the case where the  vector \\(r\\) is non-zero, this equation becomes:

\begin{equation\*}
y'(t) = f(t, y(t)) + r
\end{equation\*}


#### Usage {#usage}

```stan
matrix pmx_solve_[adams || rk45 || bdf](ODE_rhs, int nCmt, time, amt, rate, ii, evid, cmt, addl, ss, theta, [ biovar, tlag, real[,] x_r, int [,] x_i, real rel_tol, real abs_tol, int max_step, real as_rel_tol, real as_abs_tol, int as_max_step ] );
```

Here `[adams || rk45 || bdf]` indicates the
function name can use any of the three suffixes. See below.


#### Arguments {#arguments}

<!--list-separator-->

-

    ODE right-hand-side \\(f\\). It should be defined in
    `functions` block and has the following format

    ```stan
    vector = f(real t, vector y, real[] param, real[] dat_r, int[] dat_i) {...}
    ```

    Here  is time,  the unknowns of ODE,  the parameters,  the real data,
    the integer data. ,
    , and  are from
    the entry of , ,
    and  corresponding to
    , respectively.
    \\(f\\) should not include dosing rates in its
    definition, as Torsten automatically update \\(f\\)
    when the corresponding event indicates infusion dosage.

<!--list-separator-->

-

    The number of compartments, equivalently, the dimension of the ODE system.

<!--list-separator-->

-

    2d arary real data to be passed to ODE RHS. If specified, its 1st
    dimension should have the same size as .

<!--list-separator-->

-

    2d arary integer data to be passed to ODE RHS. If specified, its 1st
    dimension should have the same size as .

<!--list-separator-->

-

    The relative tolerance for numerical integration, default to 1.0E-6.

<!--list-separator-->

-

    The absolute tolerance for numerical integration, default to 1.0E-6.

<!--list-separator-->

-

    The maximum number of steps in numerical integration, default to \\(10^6\\).

<!--list-separator-->

-

    The relative tolerance for algebra solver for steady state solution, default to 1.0E-6.

<!--list-separator-->

-

    The absolute tolerance for algebra solver for steady state solution, default to 1.0E-6.

<!--list-separator-->

-

    The maximum number of interations in algebra solver for steady state solution, default to \\(10^2\\).

<!--list-separator-->

-  See Table [tab:event_args](#tab:event_args) and Table [tab:event_params](#tab:event_params) for the rest of arguments.


#### Return value {#return-value}

An `nCmt`-by-`nt` matrix, where `nt` is the size of .


#### Note {#note}

-   See section [sec:ode_func_note](#sec:ode_func_note) for different types of integrator and general guidance.
-   See section [sec:ode_func_note](#sec:ode_func_note) for comments on accuracy and tolerance.
-   The default values of ,
    , and  are
    based on a limited amount of PKPD test problems and should not be considered as
    universally applicable. We strongly recommend user to set these values
    according to physical intuition and numerical tests. See also Secion
    [sec:ode_func_note](#sec:ode_func_note).
-   With optional arguments indicated by square bracket, the following calls are allowed:

<!--listend-->

```stan
pmx_solve_[adams || rk45 || bdf](..., theta);
pmx_solve_[adams || rk45 || bdf](..., theta, rel_tol, abs_tol, max_step);
pmx_solve_[adams || rk45 || bdf](..., theta, rel_tol, abs_tol, max_step, as_rel_tol, as_abs_tol, as_max_step);
pmx_solve_[adams || rk45 || bdf](..., theta, biovar);
pmx_solve_[adams || rk45 || bdf](..., theta, biovar, rel_tol, abs_tol, max_step);
pmx_solve_[adams || rk45 || bdf](..., theta, biovar, rel_tol, abs_tol, max_step, as_rel_tol, as_abs_tol, as_max_step);
pmx_solve_[adams || rk45 || bdf](..., theta, biovar, tlag);
pmx_solve_[adams || rk45 || bdf](..., theta, biovar, tlag, rel_tol, abs_tol, max_step);
pmx_solve_[adams || rk45 || bdf](..., theta, biovar, tlag, rel_tol, abs_tol, max_step, as_rel_tol, as_abs_tol, as_max_step);
pmx_solve_[adams || rk45 || bdf](..., theta, biovar, tlag, x_r);
pmx_solve_[adams || rk45 || bdf](..., theta, biovar, tlag, x_r, rel_tol, abs_tol, max_step);
pmx_solve_[adams || rk45 || bdf](..., theta, biovar, tlag, x_r, rel_tol, abs_tol, max_step, as_rel_tol, as_abs_tol, as_max_step);
pmx_solve_[adams || rk45 || bdf](..., theta, biovar, tlag, x_r, x_i);
pmx_solve_[adams || rk45 || bdf](..., theta, biovar, tlag, x_r, x_i, rel_tol, abs_tol, max_step);
pmx_solve_[adams || rk45 || bdf](..., theta, biovar, tlag, x_r, x_i, rel_tol, abs_tol, max_step, as_rel_tol, as_abs_tol, as_max_step);
```


### Coupled ODE Model Function {#coupled-ode-model-function}



#### Description {#description}

When the ODE system consists of two subsystems in form of

\begin{align\*}
  y\_1^\prime &= f\_1(t, y\_1), \\\\\\
  y\_2^\prime &= f\_2(t, y\_1, y\_2),
\end{align\*}

with \\(y\_1\\), \\(y\_2\\), \\(f\_1\\), and \\(f\_2\\) being vector-valued functions, and
\\(y\_1^\prime\\) independent of \\(y\_2\\), the solution can be
accelerated if \\(y\_1\\) admits an analytical solution which can
be introduced into the ODE for \\(y\_2\\) for numerical
integration. This structure arises in PK/PD
models, where \\(y\_1\\) describes a forcing PK function and \\(y\_2\\) the PD
effects. In the example of a Friberg-Karlsson
semi-mechanistic model(see below), we observe an average speedup of
\\(\sim 47 \pm 18 \%\\) when using the mix solver in lieu of the numerical
integrator. In the context, currently the couple solver supports one-
& two-compartment for PK model, and  &
 integrator for nonlinear PD model.


#### Usage {#usage}

```stan
matrix pmx_solve_onecpt_[ rk45 || bdf ](reduced_ODE_system, int nOde, time, amt, rate, ii, evid, cmt, addl, ss, theta, biovar, tlag [, real rel_tol, real abs_tol, int max_step, real as_rel_tol, real as_abs_tol, int as_max_step ] );
matrix pmx_solve_twocpt_[ rk45 || bdf ](reduced_ODE_system, int nOde, time, amt, rate, ii, evid, cmt, addl, ss, theta, biovar, tlag [, real rel_tol, real abs_tol, int max_step, real as_rel_tol, real as_abs_tol, int as_max_step ] );
```


#### Arguments {#arguments}

<!--list-separator-->

-

    The system  numerically solve (\\(y\_2\\) in the above discussion, also called the
    _reduced system_ and  the number of equations in
    the \underline{reduced} system. The function that defines a reduced
    system has an almost identical signature to that used for a full
    system, but takes one additional argument: \\(y\_1\\), the PK states,
    i.e. solution to the PK ODEs.

    ```stan
    vector reduced_ODE_rhs(real t, vector y2, vector y1, real[] theta, real[] x_r, int[] x_i)
    ```

<!--list-separator-->

-

    The number of compartments. Equivalently, the dimension of the ODE system.

<!--list-separator-->

-

    The relative tolerance for numerical integration, default to 1.0E-6.

<!--list-separator-->

-

    The absolute tolerance for numerical integration, default to 1.0E-6.

<!--list-separator-->

-

    The maximum number of steps in numerical integration, default to \\(10^6\\).

<!--list-separator-->

-  See Table [tab:event_args](#tab:event_args) and Table [tab:event_params](#tab:event_params) for the rest of arguments.


#### Return value {#return-value}

An `nPk + nOde`-by-`nt` matrix, where `nt` is the size of
, and `nPk` equals to 2 in
 functions
and 3 in  functions.


### General ODE-based Population Model Function {#general-ode-based-population-model-function}


#### Description {#description}

All the preivous functions solves for a single sunject. Torsten also
provides population modeling counterparts for ODE solutions. The
functions solve for a population that share an ODE model but with
subject-level parameters and event specifications and have similar
signatures to single-subject functions, except that now events
arguments , , , ,
, ,
,  specifies the entire
population, one subject after another.


#### Usage {#usage}

```stan
matrix pmx_solve_group_[adams || rk45 || bdf](ODE_rhs, int nCmt, int[] len, time, amt, rate, ii, evid, cmt, addl, ss, theta, [ biovar, tlag, real[,] x_r, int [,] x_i, real rel_tol, real abs_tol, int max_step, real as_rel_tol, real as_abs_tol, int as_max_step ] );
```

Here `[adams || rk45 || bdf]` indicates the
function name can be of any of the three suffixes. See section [sec:ode_func_note](#sec:ode_func_note).


#### Arguments {#arguments}

<!--list-separator-->

-

    Same as in Section [sec:general_ode](#sec:general_ode).

<!--list-separator-->

-  , , , , , , ,

    2d-array arguments that describe data record for the
    entire population (see also Table [tab:event_args](#tab:event_args) and Table [tab:event_params](#tab:event_params)). They must have same size in the first
    dimension. Take  for example. Let \\(N\\) be the
    population size, then  to
     specifies events ID for subject 1,
     to
     for subject 2, etc. With \\(n\_i\\)
    being the number of events for subject \\(i\\), \\(i=1, 2, \dots, N\\), the
    size of 's first dimension is \\(\sum\_{i}n\_i\\).

<!--list-separator-->

-

    The length of data for each subject within
    the above events arrays. The size of  equals
    to population size \\(N\\).

<!--list-separator-->

-

    The number of compartments. Equivalently, the dimension of the ODE system.

<!--list-separator-->

-

    2d arary real data to be passed to ODE RHS. If specified, its 1st
    dimension should have the same size as .

<!--list-separator-->

-

    2d arary integer data to be passed to ODE RHS. If specified, its 1st
    dimension should have the same size as .

<!--list-separator-->

-

    The relative tolerance for numerical integration, default to 1.0E-6.

<!--list-separator-->

-

    The absolute tolerance for numerical integration, default to 1.0E-6.

<!--list-separator-->

-

    The maximum number of steps in numerical integration, default to \\(10^6\\).

<!--list-separator-->

-

    The relative tolerance for algebra solver for steady state solution, default to 1.0E-6.

<!--list-separator-->

-

    The absolute tolerance for algebra solver for steady state solution, default to 1.0E-6.

<!--list-separator-->

-

    The maximum number of interations in algebra solver for steady state solution, default to \\(10^2\\).


#### Return value {#return-value}

An `nCmt`-by-`nt` matrix, where `nt` is the total size of
events \\(\sum\_{i}n\_i\\).


#### Note {#note}

\label{sec:ode\_group\_note}

-   Similar to single-subject solvers, three numerical integrator are provided:
    -   : Adams-Moulton method,
    -   : Backward-differentiation formular,
    -   : Runge-Kutta 4/5 method.
-   With optional arguments indicated by square bracket, the following calls are allowed:

<!--listend-->

```stan
pmx_solve_group_[adams || rk45 || bdf](..., theta);
pmx_solve_group_[adams || rk45 || bdf](..., theta, rel_tol, abs_tol, max_step);
pmx_solve_group_[adams || rk45 || bdf](..., theta, rel_tol, abs_tol, max_step, as_rel_tol, as_abs_tol, as_max_step);
pmx_solve_group_[adams || rk45 || bdf](..., theta, biovar);
pmx_solve_group_[adams || rk45 || bdf](..., theta, biovar, rel_tol, abs_tol, max_step);
pmx_solve_group_[adams || rk45 || bdf](..., theta, biovar, rel_tol, abs_tol, max_step, as_rel_tol, as_abs_tol, as_max_step);
pmx_solve_group_[adams || rk45 || bdf](..., theta, biovar, tlag);
pmx_solve_group_[adams || rk45 || bdf](..., theta, biovar, tlag, rel_tol, abs_tol, max_step);
pmx_solve_group_[adams || rk45 || bdf](..., theta, biovar, tlag, rel_tol, abs_tol, max_step, as_rel_tol, as_abs_tol, as_max_step);
pmx_solve_group_[adams || rk45 || bdf](..., theta, biovar, tlag, x_r);
pmx_solve_group_[adams || rk45 || bdf](..., theta, biovar, tlag, x_r, rel_tol, abs_tol, max_step);
pmx_solve_group_[adams || rk45 || bdf](..., theta, biovar, tlag, x_r, rel_tol, abs_tol, max_step, as_rel_tol, as_abs_tol, as_max_step);
pmx_solve_group_[adams || rk45 || bdf](..., theta, biovar, tlag, x_r, x_i);
pmx_solve_group_[adams || rk45 || bdf](..., theta, biovar, tlag, x_r, x_i, rel_tol, abs_tol, max_step);
pmx_solve_group_[adams || rk45 || bdf](..., theta, biovar, tlag, x_r, x_i, rel_tol, abs_tol, max_step, as_rel_tol, as_abs_tol, as_max_step);
```

-   The group solvers support paralleisation through Message Passing
    Interface(MPI). One can access this feature through `cmdstan` or
    `cmdstanr` interface.

<!--listend-->

```bash
# cmdstan interface user need to add "TORSTEN_MPI=1" and
# "TBB_CXX_TYPE=gcc" in "cmdstan/make/local" file. In linux & macos
# this can be done as
echo "TORSTEN_MPI=1" > cmdstan/make/local
echo "TBB_CXX_TYPE=gcc" > cmdstan/make/local # "gcc" should be replaced by user's C compiler
make path-to-model/model_name
mpiexec -n number_of_processes model_name sample... # additional cmdstan options
```

```r
library("cmdstanr")
cmdstan_make_local(cpp_options = list("TORSTEN_MPI" = "1", "TBB_CXX_TYPE"="gcc"))  # "gcc" should be replaced by user's C compiler
rebuild_cmdstan()
mod <- cmdstan_model(path-to-model-file, quiet=FALSE, force_recompile=TRUE)
f <- mod$sample_mpi(data = ..., chains = 1, mpi_args = list("n" = number_of_processes), refresh = 200)
```

Here \\(n\\) denotes number of MPI processes, so that \\(N\\)
ODE systems (each specified by a same RHS function and
subject-dependent events) are distributed to and solved by \\(n\\)
processes evenly. Note that to access this feature user must have
MPI installed, and some MPI installation may require set additional
compiler arguments, such as `CXXLFAGS` and `LDFLAGS`.


### ODE  integrator Function {#ode-integrator-function}


#### Description {#description}

Torsten provides its own implementation of ODE solvers that solves

\begin{equation\*}
  y'(t) = f(t, y(t)), \quad y(t\_0) = y\_0
\end{equation\*}

for \\(y\\). These solvers
are customized for Torsten applications and different from those found
in Stan. The general ODE PMX solvers in previous sections are internally powered
by these functions.


#### Usage {#usage}

```stan
real[ , ] pmx_integrate_ode_[ adams || bdf || rk45 ](ODE_rhs, real[] y0, real t0, real[] ts, real[] theta, real[] x_r, int[] x_i [ , real rtol, real atol, int max_step ]);
```


#### Arguments {#arguments}

\label{sec:ode\_func\_args}

<!--list-separator-->

-

    Function that specifies the right-hand-side \\(f\\).
    It should be defined in
    `functions` block and has the following format

    ```stan
    vector = f(real t, vector y, real[] param, real[] dat_r, int[] dat_i) {...}
    ```

    Here  is time,  the unknowns of ODE,  the parameters,  the real data,
    the integer data.

<!--list-separator-->

-

    Initial condition \\(y\_0\\).

<!--list-separator-->

-

    Initial time \\(t\_0\\).

<!--list-separator-->

-

    Output time when solution is seeked.

<!--list-separator-->

-

    Parameters to be passed to  function.

<!--list-separator-->

-

    Real data to be passed to  function.

<!--list-separator-->

-

    Integer data to be passed to  function.

<!--list-separator-->

-

    Relative tolerance, default to 1.e-6() and 1.e-8( and ).

<!--list-separator-->

-

    Absolute tolerance, default to 1.e-6() and 1.e-8( and ).

<!--list-separator-->

-

    Maximum number of steps allowed between neighboring time in ,
    default to 100000.


#### Return value {#return-value}

An `n`-by-`nd` 2d-array, where `n` is the size of
and `nd` the dimension of the system.


#### Note {#note}

\label{sec:ode\_func\_note}

-   Three numerical integrator are provided:

    -   : Adams-Moulton method,
    -   : Backward-differentiation formular,
    -   : Runge-Kutta 4/5 method.

    When not equipped with further understanding of the ODE system, as a
    rule of thumb we suggest user try
     integrator first,
    integrator when the system is suspected to be stiff, and
     when a non-stiff system needs to be solved
    with higher accuracy/smaller tolerance.

-   All three integrators support adaptive stepping. To achieve
    that, at step \\(i\\) estimated error \\(e\_i\\) is calculated and
    compared with given tolerance so that

    \begin{equation}
      e\_i < \Vert\text{rtol} \times \tilde{y} + \text{atol}\Vert
    \end{equation}

    Here \\(\tilde{y}\\) is the numerical solution of \\(y\\) at current
    step and \\(\Vert \cdot \Vert\\) indicates certain norm. When the above check fails, the solver attempts
    to reduce step size and retry. The default values of ,
    , and  are
    based on Stan's ODE functions and should not be considered as
    optimal. User should make problem-dependent
    decision on  and ,
    according to estimated scale of the unknowns, so that the error
    would not affect inference on statistical variance of quantities
    that enter the Stan model. In particular, when an unknown can be neglected
    below certain threshold without affecting the rest of
    the dynamic system, setting
     greater than that threshold will avoid
    spurious and error-prone computation. See
    <sup id="22a46696f36cafa68984a1b79da368aa"><a href="#hindmarsh_cvodes_2020" title="@manual{hindmarsh_cvodes_2020,
        Author = {Hindmarsh, Alan C. and Serban, Radu and Balos, Cody
                      J. and Gardner, David J. and Woodward, Carol S. and
                      Reynolds, Daniel R.},
        Title = {User Documentation for cvodes v5.4.0},
        Year = {2020},
        }">hindmarsh_cvodes_2020</a></sup> and
    1.4 of <sup id="5b42fed7488ad403e679e0ba2797e56a"><a href="#shampine_solving_2003" title="Shampine, Gladwell, Shampine \&amp; Thompson, Solving {ODEs} with {MATLAB}, Cambridge University Press (2003).">shampine_solving_2003</a></sup> for details.

-   With optional arguments indicated by square bracket, the following calls are allowed:

<!--listend-->

```stan
pmx_integrate_ode_[adams || rk45 || bdf](..., x_i);
pmx_integrate_ode_[adams || rk45 || bdf](..., x_i, rel_tol, abs_tol, max_step);
```


### ODE group  integrator Function {#ode-group-integrator-function}



#### Description {#description}

All the preivous functions solves for a single ODE system. Torsten also
provides group modeling counterparts for ODE integrators. The
functions solve for a group of ODE systems that share an ODE RHS but with
different parameters. They have similar
signatures to single-ODE integration functions.


#### Usage {#usage}

```stan
matrix pmx_integrate_ode_group_[adams || rk45 || bdf](ODE_system, real[ , ] y0, real t0, int[] len, real[] ts, real[ , ] theta, real[ , ] x_r, int[ , ] x_i, [ real rtol, real atol, int max_step ] );
```

Here `[adams || rk45 || bdf]` indicates the
function name can be of any of the three suffixes. See section [sec:ode_func_note](#sec:ode_func_note).

<!--list-separator-->

-

    Function that specifies the right-hand-side \\(f\\). See Section [sec:ode_func_args](#sec:ode_func_args).

<!--list-separator-->

-

    Initial condition \\(y\_0\\) for each subsystem in the group. The
    first dimension equals to the size of the group.

<!--list-separator-->

-

    Initial time \\(t\_0\\).

<!--list-separator-->

-

    A vector that contains the number of output time points for each
    subsystem. The lenght of the vector equals to the size of the group.

<!--list-separator-->

-

    Output time when solution is seeked, consisting of
     of each subsystem concatenated.

<!--list-separator-->

-

    2d-array parameters to be passed to
    function. Each row corresponds to one subsystem.

<!--list-separator-->

-

    2d-array real data to be passed to  function.
    Each row corresponds to one subsystem.

<!--list-separator-->

-

    2d-array integer data to be passed to  function.
    Each row corresponds to one subsystem.

<!--list-separator-->

-

    Relative tolerance.

<!--list-separator-->

-

    Absolute tolerance.

<!--list-separator-->

-

    Maximum number of steps allowed between neighboring time in .


#### Return value {#return-value}

An `n`-by-`nd` matrix, where `n` is the size of
and `nd` the dimension of the system.


#### Note {#note}

\label{sec:ode\_integ\_group\_note}

-   With optional arguments indicated by square bracket, the following calls are allowed:

<!--listend-->

```stan
pmx_integrate_group_[adams || rk45 || bdf](..., x_i);
pmx_integrate_group_[adams || rk45 || bdf](..., x_i, rel_tol, abs_tol, max_step);
```

-   The group integrators support paralleisation through Message Passing
    Interface(MPI). One can access this feature through `cmdstan` or
    `cmdstanr` interface.

<!--listend-->

```bash
# cmdstan interface user need to add "TORSTEN_MPI=1" and
# "TBB_CXX_TYPE=gcc" in "cmdstan/make/local" file. In linux & macos
# this can be done as
echo "TORSTEN_MPI=1" > cmdstan/make/local
echo "TBB_CXX_TYPE=gcc" > cmdstan/make/local # "gcc" should be replaced by user's C compiler
make path-to-model/model_name
mpiexec -n number_of_processes model_name sample... # additional cmdstan options
```

```r
library("cmdstanr")
cmdstan_make_local(cpp_options = list("TORSTEN_MPI" = "1", "TBB_CXX_TYPE"="gcc"))  # "gcc" should be replaced by user's C compiler
rebuild_cmdstan()
mod <- cmdstan_model(path-to-model-file, quiet=FALSE, force_recompile=TRUE)
f <- mod$sample_mpi(data = ..., chains = 1, mpi_args = list("n" = number_of_processes), refresh = 200)
```

Here \\(n\\) denotes number of MPI processes, so that \\(N\\)
ODE systems are distributed to and solved by \\(n\\)
processes evenly. Note that to access this feature user must have
MPI installed, and some MPI installation may require set additional
compiler arguments, such as `CXXLFAGS` and `LDFLAGS`.


### Univariate integral {#univariate-integral}

```stan
real univariate_integral_rk45(f, t0, t1, theta, x_r, x_i)
```

```stan
real univariate_integral_bdf(f, t0, t1, theta, x_r, x_i)
```

Based on the ODE solver capability in Stan, Torsten provides functions
calculating the integral of a univariate function. The integrand function \\(f\\) must follow the signature

```stan
  real f(real t, real[] theta, real[] x_r, int[] x_i) {
    /* ... */
}
```


### Piecewise linear interpolation {#piecewise-linear-interpolation}

```stan
real linear_interpolation(real xout, real[] x, real[] y)
```

```stan
real[] linear_interpolation(real[] xout, real[] x, real[] y)
```

Torsten also provides function `linear_interpolation` for piecewise linear interpolation over a
set of x, y pairs. It returns the values of a piecewise linear
function at specified values `xout` of the first function argument. The
function is specified in terms of a set of x, y
pairs. Specifically, `linear_interpolation` implements the following function

\begin{align\*}
  y\_{\text{out}} = \left\\{\begin{array}{ll}
                 y\_1, & x\_{\text{out}} < x\_1 \\\\\\
                 y\_i + \frac{y\_{i+1} - y\_i}{x\_{i+1} - x\_i}
                 \left(x\_{\text{out}} - x\_i\right), & x\_{\text{out}} \in [x\_i, x\_{i+1}) \\\\\\
                 y\_n, & x\_{\text{out}} \ge x\_n
                          \end{array}\right.
\end{align\*}

-   The x values must be in increasing order, i.e. \\(x\_i < x\_{i+1}\\).
-   All three arguments may be data or parameters.


## Examples {#examples}

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


### Two-compartment model for single patient {#two-compartment-model-for-single-patient}

\label{sec:pk2cpt}
  We model drug absorption in a single patient and simulate plasma drug concentrations:

-   Multiple Doses: 1250 mg, every 12 hours, for a total of 15 doses
-   PK measured at 0.083, 0.167, 0.25, 0.5, 0.75, 1, 1.5, 2, 4, 6,
    8, 10 and 12 hours after 1st, 2nd, and 15th dose. In addition, the
    PK is measured every 12 hours throughout the trial.

With the plasma concentration \\(\hat{c}\\) solved from
two-compartment ODEs in , we simulate \\(c\\) according to:

\begin{align\*}
  \log\left(c\right) &\sim N\left(\log\left(\widehat{c}\right),\sigma^2\right) \\\\\\
  \left(CL, Q, V\_2, V\_3, ka\right) &= \left(5\ {\rm L/h}, 8\  {\rm L/h}, 20\  {\rm L},  70\ {\rm L}, 1.2\ {\rm h^{-1}} \right) \\\\\\
  \sigma^2 &= 0.01
\end{align\*}

The data are generated using the R package `mrgsolve` <sup id="8dd98ac45050f8bfe813328322213083"><a href="#Baron000" title="Kyle Baron \&amp; Marc Gastonguay, Simulation from ODE-Based Population PK/PD and Systems Pharmacology Models in R with mrgsolve, {Journal of Pharmacokinetics and Pharmacodynamics}, v(W-23), S84--S85 (2015).">Baron000</a></sup>.

Code below shows how Torsten function  can be used to fit the above model.

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
The MCMC history plots(Figure [2](#org69b1c2f))
suggest that the 4 chains have converged to common distributions for
all of the key model parameters. The fit to the plasma concentration
data (Figure [4](#org8a71d25)) are in close agreement with the
data, which is not surprising since the fitted model is identical to
the one used to simulate the data. Similarly the parameter posterior
density can be examined in Figure [3](#org142a44f) and shows
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

<a id="org69b1c2f"></a>

</ox-hugo/history.pdf>

<a id="org142a44f"></a>

</ox-hugo/density.pdf>

<a id="org8a71d25"></a>

</ox-hugo/ppc_ribbon.pdf>


### Two-compartment model as a linear ODE model for single patient {#two-compartment-model-as-a-linear-ode-model-for-single-patient}

\label{sec:pk2cpt\_linode}
Using , the following example fits a two-compartment model
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


### Two-compartment model solved by numerical integrator for single patient {#two-compartment-model-solved-by-numerical-integrator-for-single-patient}

\label{sec:pk2cpt\_ode}
Using , the following example fits a two-compartment model
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


### Joint PK-PD model {#joint-pk-pd-model}

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
with \\(k\_{\text{tr}}\\) the transit rate constant. In the current example, we use the two compartment model in section for
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

<a id="orga717f09"></a>

{{< figure src="/ox-hugo/neutrophilModel.jpg" caption="Figure 5: Friberg-Karlsson semi-mechanistic Model." >}}

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


### Two-compartment population model {#two-compartment-population-model}

\label{sec:twocpt\_population}
Using , the following example fits a
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


### Lotka-Volterra group model {#lotka-volterra-group-model}

\label{sec:lotka\_volterra}
Using , the following example fits
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


### Univariate integral of a quadratic function {#univariate-integral-of-a-quadratic-function}

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


### Linear intepolation {#linear-intepolation}

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


### Effect Compartment Population Model {#effect-compartment-population-model}

\label{sec:effcpt\_model}
Here we expand the example in to a population model fitted to the
combined data from phase I and phase IIa studies. The
parameters exhibit inter-individual variations (IIV), due to
both random effects and to the patients' body weight,
treated as a covariate and denoted \\(bw\\).


#### Population Model for Plasma Drug Concentration \\(c\\) {#population-model-for-plasma-drug-concentration--c}

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
a PD effect(Figure [6](#org04ec6cd)).

<a id="org04ec6cd"></a>

{{< figure src="/ox-hugo/effCptModel.png" caption="Figure 6: Effect Compartment Model" >}}


#### Effect Compartment Model for PD response \\(R\\). {#effect-compartment-model-for-pd-response--r--dot}

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
constant rate matrix and use . Correct use of
Torsten requires the user pass the entire event history
(observation and dosing events) for an individual to the
function. Thus the Stan model shows the call to
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


#### Results {#results}

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
  <span class="table-number"><a href="#table--effCptModelParms">Table 3</a></span>:
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

<a id="orgb94feb3"></a>

</ox-hugo/density.pdf>

<a id="orge672b43"></a>

</ox-hugo/ppc_study_1_5mg.pdf>

<a id="orgd031097"></a>

</ox-hugo/ppc_study_1_10mg.pdf>

<a id="org06b4275"></a>

</ox-hugo/ppc_study_1_20mg.pdf>

<a id="orgc51152c"></a>

</ox-hugo/ppc_study_1_40mg.pdf>

<a id="org214d057"></a>

</ox-hugo/ppc_study_2_20mg.pdf>

<a id="orgf85364a"></a>

</ox-hugo/ppc_study_1_5mg_resp.pdf>

<a id="org97d3ef4"></a>

</ox-hugo/ppc_study_1_10mg_resp.pdf>

<a id="org8806ca3"></a>

</ox-hugo/ppc_study_1_20mg_resp.pdf>

<a id="org314f462"></a>

</ox-hugo/ppc_study_1_40mg_resp.pdf>

<a id="orge7901bb"></a>

</ox-hugo/ppc_study_2_20mg_resp.pdf>


### Friberg-Karlsson Semi-Mechanistic Population Model {#friberg-karlsson-semi-mechanistic-population-model}

\label{sec:fkpop\_model}
We now return to the example in Section [sec:fk_model](#sec:fk_model) and extend
it to a population model. While we recommend using the coupled
solver, and this time we solve it using group solver. We leave it
as an exercise to the reader to rewrite the model with
coupled solver.


#### Friberg-Karlsson Population Model for drug-induced myelosuppression (\\(ANC\\)) {#friberg-karlsson-population-model-for-drug-induced-myelosuppression--anc}

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

We use the  function to
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


#### Results {#results}

Table [FkpopModelParms](#FkpopModelParms) summarizes the sampling and some diagnostics output.
estimation reflects the real value of the parameters (Table [FkpopModelParms](#FkpopModelParms) and Figure [fkpop_mcmc_density](#fkpop_mcmc_density).
Similar to the previous example, PPCs shown in Figure [fkpop_ppc_pk](#fkpop_ppc_pk)
and [fkpop_ppc_pd](#fkpop_ppc_pd) indicate the model is a good fit.

<a id="table--FkpopModelParms"></a>
<div class="table-caption">
  <span class="table-number"><a href="#table--FkpopModelParms">Table 4</a></span>:
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

<a id="orgadce6dd"></a>

</ox-hugo/density.pdf>

<a id="org31be84a"></a>

</ox-hugo/ppc_pk.pdf>

<a id="org1110753"></a>

</ox-hugo/ppc_pd.pdf>

\appendix


## Compiling constants {#compiling-constants}

Several constants are used in Torsten's makefile. These constants can
be used in `cmdstan/make/local` file, or use `set_make_local` command
in `cmdstanr`.

-   turns on within-chain parallelsation of MPI-enable
    functions. To use this option one must also point
     to proper C compiler. See also section [sec:ode_group_note](#sec:ode_group_note) and [sec:ode_integ_group_note](#sec:ode_integ_group_note).
-   makes BDF and Adams
    integrator use Stan's automatic differentiation to calculate
    Jacobian matrix in nonlinear solver <sup id="22a46696f36cafa68984a1b79da368aa"><a href="#hindmarsh_cvodes_2020" title="@manual{hindmarsh_cvodes_2020,
        Author = {Hindmarsh, Alan C. and Serban, Radu and Balos, Cody
                      J. and Gardner, David J. and Woodward, Carol S. and
                      Reynolds, Daniel R.},
        Title = {User Documentation for cvodes v5.4.0},
        Year = {2020},
        }">hindmarsh_cvodes_2020</a></sup>.

\printindex

\backmatter

\bibliography{torsten}

[^fn:1]: NONMEM\textregistered{} is licensed and distributed by ICON Development Solutions.
