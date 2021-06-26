+++
title = "General ODE-based Population Model Function"
author = ["Yi Zhang"]
date = 2021-06-25T00:00:00-07:00
lastmod = 2021-06-26T09:38:15-07:00
draft = false
weight = 2006
+++

<style>
  .ox-hugo-toc ul {
    list-style: none;
  }
</style>
<div class="ox-hugo-toc toc">
<div></div>

<div class="heading">Table of Contents</div>

- <span class="section-num">1</span> [Description](#description)
- <span class="section-num">2</span> [Usage](#usage)
- <span class="section-num">3</span> [Arguments](#arguments)
- <span class="section-num">4</span> [Return value](#return-value)
- <span class="section-num">5</span> [Note](#note)

</div>
<!--endtoc-->


## <span class="section-num">1</span> Description {#description}

All the preivous functions solves for a single sunject. Torsten also
provides population modeling counterparts for ODE solutions. The
functions solve for a population that share an ODE model but with
subject-level parameters and event specifications and have similar
signatures to single-subject functions, except that now events
arguments `time`, `amt`, `rate`, `ii`,
`evid`, `cmt`,
`addl`, `ss` specifies the entire
population, one subject after another.


## <span class="section-num">2</span> Usage {#usage}

```stan
matrix pmx_solve_group_[adams || rk45 || bdf](ODE_rhs, int nCmt, int[] len, time, amt, rate, ii, evid, cmt, addl, ss, theta, [ biovar, tlag, real[,] x_r, int [,] x_i, real rel_tol, real abs_tol, int max_step, real as_rel_tol, real as_abs_tol, int as_max_step ] );
```

Here `[adams || rk45 || bdf]` indicates the
function name can be of any of the three suffixes. See section [sec:ode_func_note](#sec:ode_func_note).


## <span class="section-num">3</span> Arguments {#arguments}

-   `ODE_rhs`
    Same as in Section [sec:general_ode](#sec:general_ode).
-   `time`, `amt`, `rate`, `ii`, `evid`, `cmt`, `addl`, `ss`
    2d-array arguments that describe data record for the
    entire population (see also Table [tab:event_args](#tab:event_args) and Table [tab:event_params](#tab:event_params)). They must have same size in the first
    dimension. Take `evid` for example. Let \\(N\\) be the
    population size, then `evid[1,]` to
    `evid[n1,]` specifies events ID for subject 1,
    `evid[n1 + 1,]` to
    `evid[n1 + n2,]` for subject 2, etc. With \\(n\_i\\)
    being the number of events for subject \\(i\\), \\(i=1, 2, \dots, N\\), the
    size of `evid`'s first dimension is \\(\sum\_{i}n\_i\\).
-   `len`
    The length of data for each subject within
    the above events arrays. The size of `len` equals
    to population size \\(N\\).
-   `nCmt`
    The number of compartments. Equivalently, the dimension of the ODE system.
-   `x_r`
    2d arary real data to be passed to ODE RHS. If specified, its 1st
    dimension should have the same size as `time`.
-   `x_i`
    2d arary integer data to be passed to ODE RHS. If specified, its 1st
    dimension should have the same size as `time`.
-   `rel_tol`
    The relative tolerance for numerical integration, default to 1.0E-6.
-   `abs_tol`
    The absolute tolerance for numerical integration, default to 1.0E-6.
-   `max_step`
    The maximum number of steps in numerical integration, default to \\(10^6\\).
-   `as_rel_tol`
    The relative tolerance for algebra solver for steady state solution, default to 1.0E-6.
-   `as_abs_tol`
    The absolute tolerance for algebra solver for steady state solution, default to 1.0E-6.
-   `as_max_step`
    The maximum number of interations in algebra solver for steady state solution, default to \\(10^2\\).


## <span class="section-num">4</span> Return value {#return-value}

An `nCmt`-by-`nt` matrix, where `nt` is the total size of
events \\(\sum\_{i}n\_i\\).


## <span class="section-num">5</span> Note {#note}

\label{sec:ode\_group\_note}

-   Similar to single-subject solvers, three numerical integrator are provided:
    -   `pmx_solve_group_adams`: Adams-Moulton method,
    -   `pmx_solve_group_bdf`: Backward-differentiation formular,
    -   `pmx_solve_group_rk45`: Runge-Kutta 4/5 method.
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
