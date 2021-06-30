:EXPORT<sub>HUGO</sub><sub>MENU</sub>: :menu "main" :EXPORT<sub>OPTIONS</sub>: toc:t num:nil :EXPORT<sub>DATE</sub>: <span class="timestamp-wrapper"><span class="timestamp">&lt;2021-06-25 Fri&gt;</span></span>

-   Version 0.89 <span class="timestamp-wrapper"><span class="timestamp">&lt;2021-06-15 Tue&gt;</span></span>
    -   Changed
        -   New backend for ODE events solvers.
        -   Use vector instead of array as ODE function state & return type.
        -   Simplified ODE integrator naming, e.g. `pmx_ode_bdf[_ctrl]`.
        -   Update to Stan version 2.27.0.
-   Version 0.88 <span class="timestamp-wrapper"><span class="timestamp">&lt;2020-12-18 Fri&gt;</span></span>
    
    -   Added
        -   Bioavailability, lag time, ODE real & integer data are optional in PMX function signatures.
        -   Support all EVID options from NM-TRAN and mrgsolve.
        -   Support steady-state infusion through multiple interdose intervals.
    
    -   Changed
        -   More efficient memory management of COVDES implenmentation.
        -   Update of MPI framework to adapt multilevel paralleism.
        -   Update to Stan version 2.25.0.
        -   Use cmdstanr as R interface.
        -   Stop supporting rstan as R interface.
-   Version 0.87 <span class="timestamp-wrapper"><span class="timestamp">&lt;2019-07-26 Fri&gt;</span></span>
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
-   Version 0.86 <span class="timestamp-wrapper"><span class="timestamp">&lt;2019-05-15 Wed&gt;</span></span>
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

-   Version 0.85 <span class="timestamp-wrapper"><span class="timestamp">&lt;2018-12-04 Tue&gt;</span></span>
    -   Added
        -   Dosing rate as parameter
    -   Changed
    -   Update to Stan version 2.18.0.

-   Version 0.84 <span class="timestamp-wrapper"><span class="timestamp">&lt;2018-02-24 Sat&gt;</span></span>
    -   Added
        -   Piecewise linear interpolation function.
        -   Univariate integral functions.
    -   Changed
        -   Update to Stan version 2.17.1.
        -   Minor revisions to User Manual.
        -   Bugfixes.
-   Version 0.83 <span class="timestamp-wrapper"><span class="timestamp">&lt;2017-08-02 Wed&gt;</span></span>
    -   Added
        -   Work with TorstenHeaders
        -   Each chain has a different initial estimate
    -   Changed
        -   User manual
        -   Fix misspecification in ODE system for TwoCpt example.
        -   Other bugfixes
-   Version 0.82 <span class="timestamp-wrapper"><span class="timestamp">&lt;2017-01-29 Sun&gt;</span></span>
    -   Added
        -   Allow parameter arguments to be passed as 1D or 2D arrays
        -   More unit tests
        -   Unit tests check automatic differentiation against finite differentiation.
    -   Changed
        -   Split the parameter argument into three arguments: pMatrix (parameters for the ODEs &#x2013; note: for `linOdeModel`, pMatrix is replaced by the constant rate matrix K), biovar (parameters for the biovariability), and tlag (parameters for the lag time).
        -   bugfixes
-   Version 0.81 <span class="timestamp-wrapper"><span class="timestamp">&lt;2016-09-27 Tue&gt;</span></span>
    -   Added
        -   linCptModel (linear compartmental model) function
-   Version 0.80a <span class="timestamp-wrapper"><span class="timestamp">&lt;2016-09-21 Wed&gt;</span></span>
    -   Added
        -   check_finite statements in pred_1 and pred_2 to reject metropolis proposal if initial conditions are not finite
