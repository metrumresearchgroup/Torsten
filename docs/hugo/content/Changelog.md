+++
title = "Changelog"
author = ["Yi Zhang"]
date = 2021-06-25T00:00:00-07:00
lastmod = 2021-06-26T00:13:54-07:00
draft = false
weight = 2004
+++

## <span class="section-num">1</span> Version 0.89 <span class="timestamp-wrapper"><span class="timestamp">&lt;2021-06-15 Tue&gt;</span></span> {#version-0-dot-89}


### <span class="section-num">1.1</span> Changed {#0-85-changed}

-   New backend for ODE events solvers.
-   Use vector instead of array as ODE function state & return type.
-   Simplified ODE integrator naming,
    e.g. `pmx_ode_bdf` and  `pmx_ode_bdf_ctrl`.
-   Update to Stan version 2.27.0.


## <span class="section-num">2</span> Version 0.88 <span class="timestamp-wrapper"><span class="timestamp">&lt;2020-12-18 Fri&gt;</span></span> {#version-0-dot-88}


### <span class="section-num">2.1</span> Added {#0-85-added}

-   Bioavailability, lag time, ODE real & integer data are optional in PMX function signatures.
-   Support all EVID options from NM-TRAN and mrgsolve.
-   Support steady-state infusion through multiple interdose intervals.


### <span class="section-num">2.2</span> Changed {#0-85-changed}

-   More efficient memory management of COVDES implenmentation.
-   Update of MPI framework to adapt multilevel paralleism.
-   Update to Stan version 2.25.0.
-   Use cmdstanr as R interface.
-   Stop supporting rstan as R interface.


## <span class="section-num">3</span> Version 0.87 <span class="timestamp-wrapper"><span class="timestamp">&lt;2019-07-26 Fri&gt;</span></span> {#version-0-dot-87}


### <span class="section-num">3.1</span> Added {#0-85-added}

-   MPI dynamic load balance for Torsten's population ODE integrators

    -   `pmx_integrate_ode_group_adams`
    -   `pmx_integrate_ode_group_bdf`
    -   `pmx_integrate_ode_group_rk45`

    To invoke dynamic load balance instead of default static
    balance for MPI, issue `TORSTEN_MPI=2` in `make/local`.
-   Support `RATE` as parameter in `pmx_solve_rk45/bdf/adams`
    functions.


### <span class="section-num">3.2</span> Changed {#0-85-changed}

-   Some fixes on steady-state solvers
-   Update to rstan version 2.19.2.


## <span class="section-num">4</span> Version 0.86 <span class="timestamp-wrapper"><span class="timestamp">&lt;2019-05-15 Wed&gt;</span></span> {#version-0-dot-86}


### <span class="section-num">4.1</span> Added {#0-85-added}

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


### <span class="section-num">4.2</span> Changed {#0-85-changed}

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


## <span class="section-num">5</span> Version 0.85 <span class="timestamp-wrapper"><span class="timestamp">&lt;2018-12-04 Tue&gt;</span></span> {#version-0-dot-85}


### <span class="section-num">5.1</span> Added {#0-85-added}

-   Dosing rate as parameter


### <span class="section-num">5.2</span> Changed {#0-85-changed}

-   Update to Stan version 2.18.0.


## <span class="section-num">6</span> Version 0.84 <span class="timestamp-wrapper"><span class="timestamp">&lt;2018-02-24 Sat&gt;</span></span> {#version-0-dot-84}


### <span class="section-num">6.1</span> Added {#0-84-added}

-   Piecewise linear interpolation function.
-   Univariate integral functions.


### <span class="section-num">6.2</span> Changed {#0-84-changed}

-   Update to Stan version 2.17.1.
-   Minor revisions to User Manual.
-   Bugfixes.


## <span class="section-num">7</span> Version 0.83 <span class="timestamp-wrapper"><span class="timestamp">&lt;2017-08-02 Wed&gt;</span></span> {#version-0-dot-83}


### <span class="section-num">7.1</span> Added {#0-83-added}

-   Work with TorstenHeaders
-   Each chain has a different initial estimate


### <span class="section-num">7.2</span> Changed {#0-83-changed}

-   User manual
-   Fix misspecification in ODE system for TwoCpt example.
-   Other bugfixes


## <span class="section-num">8</span> Version 0.82 <span class="timestamp-wrapper"><span class="timestamp">&lt;2017-01-29 Sun&gt;</span></span> {#0-82-added}


### <span class="section-num">8.1</span> Added {#added}

-   Allow parameter arguments to be passed as 1D or 2D arrays
-   More unit tests
-   Unit tests check automatic differentiation against finite differentiation.


### <span class="section-num">8.2</span> Changed {#0-82-changed}

-   Split the parameter argument into three arguments: pMatrix
    (parameters for the ODEs -- note: for `linOdeModel`, pMatrix
    is replaced by the constant rate matrix K), biovar
    (parameters for the biovariability), and tlag (parameters
    for the lag time).
-   bugfixes


## <span class="section-num">9</span> Version 0.81 <span class="timestamp-wrapper"><span class="timestamp">&lt;2016-09-27 Tue&gt;</span></span> {#version-0-dot-81}


### <span class="section-num">9.1</span> Added {#0-81-added}

linCptModel (linear compartmental model) function


## <span class="section-num">10</span> Version 0.80a <span class="timestamp-wrapper"><span class="timestamp">&lt;2016-09-21 Wed&gt;</span></span> {#version-0-dot-80a}


### <span class="section-num">10.1</span> Added {#0-80a-added}

check_finite statements in pred_1 and pred_2 to reject metropolis proposal if initial conditions are not finite
