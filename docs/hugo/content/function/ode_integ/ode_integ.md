+++
title = "ODE  integrator Function"
author = ["Yi Zhang"]
date = 2021-06-25T00:00:00-07:00
lastmod = 2021-06-26T09:38:15-07:00
draft = false
weight = 2007
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

Torsten provides its own implementation of ODE solvers that solves

\begin{equation\*}
  y'(t) = f(t, y(t)), \quad y(t\_0) = y\_0
\end{equation\*}

for \\(y\\). These solvers
are customized for Torsten applications and different from those found
in Stan. The general ODE PMX solvers in previous sections are internally powered
by these functions.


## <span class="section-num">2</span> Usage {#usage}

```stan
real[ , ] pmx_integrate_ode_[ adams || bdf || rk45 ](ODE_rhs, real[] y0, real t0, real[] ts, real[] theta, real[] x_r, int[] x_i [ , real rtol, real atol, int max_step ]);
```


## <span class="section-num">3</span> Arguments {#arguments}

\label{sec:ode\_func\_args}

-   `ODE_rhs`
    Function that specifies the right-hand-side \\(f\\).
    It should be defined in
    `functions` block and has the following format

<!--listend-->

```stan
vector = f(real t, vector y, real[] param, real[] dat_r, int[] dat_i) {...}
```

Here `t` is time, `y` the unknowns of ODE, `param` the parameters, `dat\_r` the real data, `dat\_i`
the integer data.

-   `y0`
    Initial condition \\(y\_0\\).
-   `t0`
    Initial time \\(t\_0\\).
-   `ts`
    Output time when solution is seeked.
-   `theta`
    Parameters to be passed to `ODE_rhs` function.
-   `x_r`
    Real data to be passed to `ODE_rhs` function.
-   `x_i`
    Integer data to be passed to `ODE_rhs` function.
-   `rtol`
    Relative tolerance, default to 1.e-6(`rk45`) and 1.e-8(`adams` and `bdf`).
-   `atol`
    Absolute tolerance, default to 1.e-6(`rk45`) and 1.e-8(`adams` and `bdf`).
-   `max_step`
    Maximum number of steps allowed between neighboring time in `ts`,
    default to 100000.


## <span class="section-num">4</span> Return value {#return-value}

An `n`-by-`nd` 2d-array, where `n` is the size of `ts`
and `nd` the dimension of the system.


## <span class="section-num">5</span> Note {#note}

\label{sec:ode\_func\_note}

-   Three numerical integrator are provided:

    -   `pmx_integrate_ode_adams`: Adams-Moulton method,
    -   `pmx_integrate_ode_bdf`: Backward-differentiation formular,
    -   `pmx_integrate_ode_rk45`: Runge-Kutta 4/5 method.

    When not equipped with further understanding of the ODE system, as a
    rule of thumb we suggest user try
    `rk45` integrator first, `bdf`
    integrator when the system is suspected to be stiff, and
    `adams` when a non-stiff system needs to be solved
    with higher accuracy/smaller tolerance.

-   All three integrators support adaptive stepping. To achieve
    that, at step \\(i\\) estimated error \\(e\_i\\) is calculated and
    compared with given tolerance so that

    \begin{equation}
      e\_i < \Vert\text{rtol} \times \tilde{y} + \text{atol}\Vert
    \end{equation}

    Here \\(\tilde{y}\\) is the numerical solution of \\(y\\) at current
    step and \\(\Vert \cdot \Vert\\) indicates certain norm. When the above check fails, the solver attempts
    to reduce step size and retry. The default values of `atol`,
    `rtol`, and `max_step` are
    based on Stan's ODE functions and should not be considered as
    optimal. User should make problem-dependent
    decision on `rtol` and `atol`,
    according to estimated scale of the unknowns, so that the error
    would not affect inference on statistical variance of quantities
    that enter the Stan model. In particular, when an unknown can be neglected
    below certain threshold without affecting the rest of
    the dynamic system, setting
    `atol` greater than that threshold will avoid
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
