+++
title = "ODE  integrator function"
author = ["Yi Zhang"]
date = 2021-06-28T00:00:00-07:00
lastmod = 2021-06-29T14:35:27-07:00
draft = false
weight = 2008
+++

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


#### Return value {#return-value}

An `n`-by-`nd` 2d-array, where `n` is the size of `ts`
and `nd` the dimension of the system.


#### Note {#note}

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
    ([Hindmarsh et al. 2020](#orga5af84e)) and
    1.4 of ([Shampine et al. 2003](#orgb80b28d)) for details.

-   With optional arguments indicated by square bracket, the following calls are allowed:

<!--listend-->

```stan
pmx_integrate_ode_[adams || rk45 || bdf](..., x_i);
pmx_integrate_ode_[adams || rk45 || bdf](..., x_i, rel_tol, abs_tol, max_step);
```


## Bibliography {#bibliography}

<a id="orga5af84e"></a>Hindmarsh, Alan C., Radu Serban, Cody J. Balos, David J. Gardner, Carol S. Woodward, and Daniel R. Reynolds. 2020. _User Documentation for Cvodes V5.4.0_.

<a id="orgb80b28d"></a>Shampine, L. F., I. Gladwell, Larry Shampine, and S. Thompson. 2003. _Solving ODEs with MATLAB_. Cambridge University Press.
