+++
title = "ODE group  integrator Function"
author = ["Yi Zhang"]
date = 2021-06-25T00:00:00-07:00
lastmod = 2021-06-28T20:40:58-07:00
draft = false
weight = 2008
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
- <span class="section-num">3</span> [Return value](#return-value)
- <span class="section-num">4</span> [Note](#note)

</div>
<!--endtoc-->



## <span class="section-num">1</span> Description {#description}

All the preivous functions solves for a single ODE system. Torsten also
provides group modeling counterparts for ODE integrators. The
functions solve for a group of ODE systems that share an ODE RHS but with
different parameters. They have similar
signatures to single-ODE integration functions.


## <span class="section-num">2</span> Usage {#usage}

```stan
matrix pmx_integrate_ode_group_[adams || rk45 || bdf](ODE_system, real[ , ] y0, real t0, int[] len, real[] ts, real[ , ] theta, real[ , ] x_r, int[ , ] x_i, [ real rtol, real atol, int max_step ] );
```

Here `[adams || rk45 || bdf]` indicates the
function name can be of any of the three suffixes. See section [sec:ode_func_note](#sec:ode_func_note).

-   `ODE_rhs`
    Function that specifies the right-hand-side \\(f\\). See Section [sec:ode_func_args](#sec:ode_func_args).
-   `y0`
    Initial condition \\(y\_0\\) for each subsystem in the group. The
    first dimension equals to the size of the group.
-   `t0`
    Initial time \\(t\_0\\).
-   `len`
    A vector that contains the number of output time points for each
    subsystem. The lenght of the vector equals to the size of the group.
-   `ts`
    Output time when solution is seeked, consisting of
    `ts` of each subsystem concatenated.
-   `theta`
    2d-array parameters to be passed to `ODE_rhs`
    function. Each row corresponds to one subsystem.
-   `x_r`
    2d-array real data to be passed to `ODE_rhs` function.
    Each row corresponds to one subsystem.
-   `x_i`
    2d-array integer data to be passed to `ODE_rhs` function.
    Each row corresponds to one subsystem.
-   `rtol`
    Relative tolerance.
-   `atol`
    Absolute tolerance.
-   `max_step`
    Maximum number of steps allowed between neighboring time in `ts`.


## <span class="section-num">3</span> Return value {#return-value}

An `n`-by-`nd` matrix, where `n` is the size of `ts`
and `nd` the dimension of the system.


## <span class="section-num">4</span> Note {#note}

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
