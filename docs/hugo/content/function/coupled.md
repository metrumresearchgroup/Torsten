+++
title = "Coupled ODE Model Function"
author = ["Yi Zhang"]
date = 2021-06-25T00:00:00-07:00
lastmod = 2021-06-29T11:51:47-07:00
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

</div>
<!--endtoc-->



## <span class="section-num">1</span> Description {#description}

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
& two-compartment for PK model, and `rk45` &
`bdf` integrator for nonlinear PD model.


## <span class="section-num">2</span> Usage {#usage}

```stan
matrix pmx_solve_onecpt_[ rk45 || bdf ](reduced_ODE_system, int nOde, time, amt, rate, ii, evid, cmt, addl, ss, theta, biovar, tlag [, real rel_tol, real abs_tol, int max_step, real as_rel_tol, real as_abs_tol, int as_max_step ] );
matrix pmx_solve_twocpt_[ rk45 || bdf ](reduced_ODE_system, int nOde, time, amt, rate, ii, evid, cmt, addl, ss, theta, biovar, tlag [, real rel_tol, real abs_tol, int max_step, real as_rel_tol, real as_abs_tol, int as_max_step ] );
```


## <span class="section-num">3</span> Arguments {#arguments}

-   `reduced\_ODE\_rhs`
    The system  numerically solve (\\(y\_2\\) in the above discussion, also called the
    _reduced system_ and `nOde` the number of equations in
    the \underline{reduced} system. The function that defines a reduced
    system has an almost identical signature to that used for a full
    system, but takes one additional argument: \\(y\_1\\), the PK states,
    i.e. solution to the PK ODEs.

    ```stan
    vector reduced_ODE_rhs(real t, vector y2, vector y1, real[] theta, real[] x_r, int[] x_i)
    ```
-   `nCmt`
    The number of compartments. Equivalently, the dimension of the ODE system.
-   `rel_tol`
    The relative tolerance for numerical integration, default to 1.0E-6.
-   `abs_tol`
    The absolute tolerance for numerical integration, default to 1.0E-6.
-   `max_step`
    The maximum number of steps in numerical integration, default to \\(10^6\\).
-   See Tables in Section [Events specification]({{< relref "events" >}}) for the rest of arguments.


## <span class="section-num">4</span> Return value {#return-value}

An `(nPk + nOde)` &times; `nt` matrix, where `nt` is the size of
`time`, and `nPk` equals to 2 in
`pmx_solve_onecpt_` functions
and 3 in `pmx_solve_twocpt_` functions.
