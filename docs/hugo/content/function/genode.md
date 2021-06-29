+++
title = "General ODE Model Function"
author = ["Yi Zhang"]
date = 2021-06-25T00:00:00-07:00
lastmod = 2021-06-29T14:35:26-07:00
draft = false
weight = 2005
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

Function `pmx_solve_adams`, `pmx_solve_bdf`, and `pmx_solve_rk45` solve a first-order ODE system
specified by user-specified right-hand-side function `ODE_rhs` \\(f\\)

\begin{equation\*}
y'(t) = f(t, y(t))
\end{equation\*}

In the case where the `rate` vector \\(r\\) is non-zero, this equation becomes:

\begin{equation\*}
y'(t) = f(t, y(t)) + r
\end{equation\*}


## <span class="section-num">2</span> Usage {#usage}

```stan
matrix pmx_solve_[adams || rk45 || bdf](ODE_rhs, int nCmt, time, amt, rate, ii, evid, cmt, addl, ss, theta, [ biovar, tlag, real[,] x_r, int [,] x_i, real rel_tol, real abs_tol, int max_step, real as_rel_tol, real as_abs_tol, int as_max_step ] );
```

Here `[adams || rk45 || bdf]` indicates the
function name can use any of the three suffixes. See below.


## <span class="section-num">3</span> Arguments {#arguments}

-   `ODE_rhs`
    ODE right-hand-side \\(f\\). It should be defined in
    `functions` block and has the following format

<!--listend-->

```stan
vector = f(real t, vector y, real[] param, real[] dat_r, int[] dat_i) {...}
```

Here `t` is time, `y` the unknowns of ODE, `param` the parameters, `dat\_r` the real data, `dat\_i`
the integer data. `param`,
`dat\_r`, and `dat\_i` are from
the entry of `theta`, `x_r`,
and `x_i` corresponding to
`t`, respectively.
\\(f\\) should not include dosing rates in its
definition, as Torsten automatically update \\(f\\)
when the corresponding event indicates infusion dosage.

-   `nCmt`
    The number of compartments, equivalently, the dimension of the ODE system.
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
-   See Tables in Section [Events specification]({{< relref "events" >}}) for the rest of arguments.


## <span class="section-num">4</span> Return value {#return-value}

An `nCmt`-by-`nt` matrix, where `nt` is the size of `time`.


## <span class="section-num">5</span> Note {#note}

-   See Section [ODE  integrator function]({{< relref "ode-integ" >}}) for different types of integrator and general guidance.
-   See Section [ODE  integrator function]({{< relref "ode-integ" >}}) for comments on accuracy and tolerance.
-   The default values of `atol`,
    `rtol`, and `max_step` are
    based on a limited amount of PKPD test problems and should not be considered as
    universally applicable. We strongly recommend user to set these values
    according to physical intuition and numerical tests. See also Section [ODE  integrator function]({{< relref "ode-integ" >}}).
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
