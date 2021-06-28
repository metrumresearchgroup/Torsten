+++
title = "Introduction"
author = ["Yi Zhang"]
date = 2021-06-25T00:00:00-07:00
lastmod = 2021-06-28T14:19:24-07:00
draft = false
weight = 1003
+++

Stan is an open source probabilistic programing language designed
primarily to do Bayesian data analysis
([Carpenter et al. 2017](#org535d69f)). Several of its
features make it a powerful tool to specify and fit complex
models. First, its language is very expressive and
flexible. Secondly, it implements a variant of No U-Turn
Sampler(NUTS), an adaptative Hamiltonian Monte Carlo
algorithm that was proven more efficient than commonly used Monte Carlo Markov Chains
(MCMC) samplers for complex high dimensional problems ([Hoffman and Gelman 2011](#orgca0d9d2); [Betancourt 2018](#org6c7e7f4)). Our
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


### <span class="section-num">0.1</span> Overview {#overview}

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


### <span class="section-num">0.2</span> Implementation summary {#implementation-summary}

-   Current Torsten v0.89rc is based on Stan v2.27.0.
-   All functions are programmed in C++ and are compatible
    with the Stan math automatic differentiation library ([Carpenter et al. 2015](#org177c092))
-   One and two compartment models are based on analytical solutions of governing ODEs.
-   General linear compartment models are based on semi-analytical solutions using the built-in matrix exponential function
-   General compartment models are solved numerically using built-in ODE integrators in Stan. The tuning parameters of the solver are adjustable. The steady state solution is calculated using a numerical algebraic solver.
-   A mix compartment model's PK forcing function is solved analytically, and its forced ODE system is solved numerically.


### <span class="section-num">0.3</span> Development plans {#development-plans}

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


### <span class="section-num">0.4</span> Changelog {#changelog}


#### <span class="section-num">0.4.1</span> Version 0.89 <span class="timestamp-wrapper"><span class="timestamp">&lt;2021-06-15 Tue&gt;</span></span> {#version-0-dot-89}

<!--list-separator-->

-  Changed

    -   New backend for ODE events solvers.
    -   Use vector instead of array as ODE function state & return type.
    -   Simplified ODE integrator naming,
        e.g. `pmx_ode_bdf` and  `pmx_ode_bdf_ctrl`.
    -   Update to Stan version 2.27.0.


#### <span class="section-num">0.4.2</span> Version 0.88 <span class="timestamp-wrapper"><span class="timestamp">&lt;2020-12-18 Fri&gt;</span></span> {#version-0-dot-88}

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


#### <span class="section-num">0.4.3</span> Version 0.87 <span class="timestamp-wrapper"><span class="timestamp">&lt;2019-07-26 Fri&gt;</span></span> {#version-0-dot-87}

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


#### <span class="section-num">0.4.4</span> Version 0.86 <span class="timestamp-wrapper"><span class="timestamp">&lt;2019-05-15 Wed&gt;</span></span> {#version-0-dot-86}

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


#### <span class="section-num">0.4.5</span> Version 0.85 <span class="timestamp-wrapper"><span class="timestamp">&lt;2018-12-04 Tue&gt;</span></span> {#version-0-dot-85}

<!--list-separator-->

-  Added

    -   Dosing rate as parameter

<!--list-separator-->

-  Changed

    -   Update to Stan version 2.18.0.


#### <span class="section-num">0.4.6</span> Version 0.84 <span class="timestamp-wrapper"><span class="timestamp">&lt;2018-02-24 Sat&gt;</span></span> {#version-0-dot-84}

<!--list-separator-->

-  Added

    -   Piecewise linear interpolation function.
    -   Univariate integral functions.

<!--list-separator-->

-  Changed

    -   Update to Stan version 2.17.1.
    -   Minor revisions to User Manual.
    -   Bugfixes.


#### <span class="section-num">0.4.7</span> Version 0.83 <span class="timestamp-wrapper"><span class="timestamp">&lt;2017-08-02 Wed&gt;</span></span> {#version-0-dot-83}

<!--list-separator-->

-  Added

    -   Work with TorstenHeaders
    -   Each chain has a different initial estimate

<!--list-separator-->

-  Changed

    -   User manual
    -   Fix misspecification in ODE system for TwoCpt example.
    -   Other bugfixes


#### <span class="section-num">0.4.8</span> Version 0.82 <span class="timestamp-wrapper"><span class="timestamp">&lt;2017-01-29 Sun&gt;</span></span> {#0-82-added}

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


#### <span class="section-num">0.4.9</span> Version 0.81 <span class="timestamp-wrapper"><span class="timestamp">&lt;2016-09-27 Tue&gt;</span></span> {#version-0-dot-81}

<!--list-separator-->

-  Added

    linCptModel (linear compartmental model) function


#### <span class="section-num">0.4.10</span> Version 0.80a <span class="timestamp-wrapper"><span class="timestamp">&lt;2016-09-21 Wed&gt;</span></span> {#version-0-dot-80a}

<!--list-separator-->

-  Added

    check_finite statements in pred_1 and pred_2 to reject metropolis proposal if initial conditions are not finite


## <span class="section-num">1</span> Bibliography {#bibliography}

<a id="org6c7e7f4"></a>Betancourt, Michael. 2018. “A Conceptual Introduction to Hamiltonian Monte Carlo.” <https://arxiv.org/abs/1701.02434>.

<a id="org535d69f"></a>Carpenter, Bob, Andrew Gelman, Matthew D. Hoffman, Daniel Lee, Ben Goodrich, Michael Betancourt, Marcus Brubaker, Jiqiang Guo, Peter Li, and Allen Riddell. 2017. “Stan: A Probabilistic Programming Language.” _Journal of Statistical Software_ 76.

<a id="org177c092"></a>Carpenter, Bob, Matthew D. Hoffman, Marcus Brubaker, Daniel Lee, Peter Li, and Michael Betancourt. 2015. “The Stan Math Library: Reverse-Mode Automatic Differentiation in C++.” _arXiv:1509.07164 [Cs]_, September. <http://arxiv.org/abs/1509.07164>.

<a id="orgca0d9d2"></a>Hoffman, Matthew D., and Andrew Gelman. 2011. “The No-U-Turn Sampler: Adaptively Setting Path Lengths in Hamiltonian Monte Carlo.” _arXiv:1111.4246 [Cs, Stat]_, November. <http://arxiv.org/abs/1111.4246>.

[^fn:1]: NONMEM\textregistered{} is licensed and distributed by ICON Development Solutions.
