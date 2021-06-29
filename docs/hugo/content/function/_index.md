+++
title = "Using Torsten"
author = ["Yi Zhang"]
date = 2021-06-25T00:00:00-07:00
lastmod = 2021-06-29T11:51:43-07:00
draft = false
weight = 1006
+++

The reader should have a basic understanding of [how Stan works](https://mc-stan.org/users/documentation/).
In this section we go through the different functions Torsten adds to
Stan. The code for the examples can be found at Torsten's [example models](https://github.com/metrumresearchgroup/Torsten/tree/master/example-models).


### <span class="section-num">0.1</span> Events specification {#events-specification}

Torsten's functions are prefixed with `pmx_`.
For some of their arguments we adopt NM-TRAN format for events
specification(Table [1](#table--tab:event-args)).

<a id="table--tab:event-args"></a>
<div class="table-caption">
  <span class="table-number"><a href="#table--tab:event-args">Table 1</a></span>:
  NM-TRAN compatible event specification arguments. All arrays should have the same length corresponding to the number of events.
</div>

| Argument Name | Definition                  | Stan data type |
|---------------|-----------------------------|----------------|
| `time`        | event time                  | `real[]`       |
| `amt`         | dosing amount               | `real[]`       |
| `rate`        | infusion rate               | `real[]`       |
| `ii`          | interdose interval          | `real[]`       |
| `evid`        | event ID                    | `int[]`        |
| `cmt`         | event compartment           | `int[]`        |
| `addl`        | additionial identical doses | `int[]`        |
| `ss`          | steady-state dosing flag    | `int[]`        |

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
  PMX model parameter overloadings. One can use 1d array <code class="src src-stan"><span style="color: #b58900;">real</span>[]</code> to indicate constants of all events, or 2d array <code class="src src-stan"><span style="color: #b58900;">real</span>[ , ]</code> so that the \(i\)th row of the array describes the model arguments for time interval \((t_{i-1}, t_i)\), and the number of the rows equals to the size of <code>time</code>.
</div>

| Argument Name | Definition               | Stan data type          | Optional           |
|---------------|--------------------------|-------------------------|--------------------|
| `theta`       | model parameters         | `real[]` or `real[ , ]` | N                  |
| `biovar`      | bioavailability fraction | `real[]` or `real[ , ]` | Y (default to 1.0) |
| `tlag`        | lag time                 | `real[]` or `real[ , ]` | Y (default to 0.0) |


### <span class="section-num">0.2</span> One Compartment Model {#one-compartment-model}



#### <span class="section-num">0.2.1</span> Description {#description}

Function `pmx_solve_onecpt` solves a one-compartment PK
model (Figure [1](#org65164c8)). The model obtains plasma concentrations of parent drug \\(c=y\_2/V\_2\\)
by solving for the mass of drug in the central compartment
\\(y\_2\\) from ordinary differential equations(ODEs)

\begin{align}\label{eq:onecpt}
  y\_1' &= -k\_a y\_1, \\\\\\
  y\_2' &= k\_a y\_1 - \left(\frac{CL}{V\_2} + \frac{Q}{V\_2}\right) y\_2.
\end{align}

<a id="org65164c8"></a>

{{< figure src="https://raw.githubusercontent.com/metrumresearchgroup/Torsten/master/docs/graphics/cptModels.png" caption="Figure 1: One and two compartment models with first order absorption implemented in Torsten." width="300" >}}


#### <span class="section-num">0.2.2</span> Usage {#usage}

```stan
matrix = pmx_solve_onecpt(time, amt, rate, ii, evid, cmt, addl, ss, theta [, biovar, tlag ] )
```


#### <span class="section-num">0.2.3</span> Arguments {#arguments}

See Tables in Section [Events specification]({{< relref "events" >}}).


#### <span class="section-num">0.2.4</span> Return value {#return-value}

An `ncmt`-by-`nt` matrix, where `nt` is the number of time steps and `ncmt=2` is the number of compartments.


#### <span class="section-num">0.2.5</span> Note {#note}

-   ODE Parameters `theta` should consist of \\(CL\\), \\(V\_2\\), \\(k\_a\\), in that order.
-   `biovar` and `tlag` are optional, so that the following are allowed:

<!--listend-->

```stan
pmx_solve_onecpt(..., theta);
pmx_solve_onecpt(..., theta, biovar);
pmx_solve_onecpt(..., theta, biovar, tlag);
```

-   Setting \\(k\_a = 0\\) eliminates the first-order absorption.


### <span class="section-num">0.3</span> Two Compartment Model {#two-compartment-model}


#### <span class="section-num">0.3.1</span> Description {#description}

Function `pmx_solve_twocpt` solves a two-compartment PK
model (Figure [{{< relref "one-cpt" >}}]({{< relref "one-cpt" >}})). The model obtains plasma concentrations of parent drug \\(c=y\_2/V\_2\\)
by solving for the mass of drug in the central compartment
\\(y\_2\\) from ordinary differential equations(ODEs)

\begin{align} \label{eq:twocpt}
  y\_1' &= -k\_a y\_1 \\\\\\
  y\_2' &= k\_a y\_1 - \left(\frac{CL}{V\_2} + \frac{Q}{V\_2}\right) y\_2 +  \frac{Q}{V\_3}  y\_3  \\\\\\
  y\_3' &= \frac{Q}{V\_2} y\_2 - \frac{Q}{V\_3} y\_3
\end{align}


#### <span class="section-num">0.3.2</span> Usage {#usage}

```stan
matrix = pmx_solve_twocpt(time, amt, rate, ii, evid, cmt, addl, ss, theta [, biovar, tlag ] )
```


#### <span class="section-num">0.3.3</span> Arguments {#arguments}

See Tables in Section [Events specification]({{< relref "events" >}}).


#### <span class="section-num">0.3.4</span> Return value {#return-value}

An `ncmt`-by-`nt` matrix, where `nt` is the number of time steps and `ncmt=3` is the number of compartments.


#### <span class="section-num">0.3.5</span> Note {#note}

-   ODE Parameters `theta` consists of \\(CL\\), \\(Q\\), \\(V\_2\\), \\(V\_3\\), \\(k\_a\\).
-   `biovar` and `tlag` are optional, so that the following are allowed:

<!--listend-->

```stan
pmx_solve_twocpt(..., theta);
pmx_solve_twocpt(..., theta, biovar);
pmx_solve_twocpt(..., theta, biovar, tlag);
```

-   Setting \\(k\_a = 0\\) eliminates the first-order absorption.


### <span class="section-num">0.4</span> General Linear ODE Model Function {#general-linear-ode-model-function}


#### <span class="section-num">0.4.1</span> Description {#description}

Function `pmx_solve_linode` solves a (piecewise) linear ODEs model with coefficients
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


#### <span class="section-num">0.4.2</span> Usage {#usage}

```stan
matrix = pmx_solve_linode(time, amt, rate, ii, evid, cmt, addl, ss, K, biovar, tlag )
```


#### <span class="section-num">0.4.3</span> Arguments {#arguments}

-   `K`
    System parameters. `K` can be either
    -   a `matrix` for constant parameters in all events, or
    -   an array of matrices `matrix K[nt]` so that the \\(i\\)th entry of the array describes
        the model parameters for time interval \\((t\_{i-1}, t\_i)\\),
        and the number of the rows equals to the number of event time `nt`.
    -   See Tables in Section [Events specification]({{< relref "events" >}}) for the rest of arguments.


#### <span class="section-num">0.4.4</span> Return value {#return-value}

An `n`-by-`nt` matrix, where `nt` is the number of time steps and `n` is the number of rows(columns) of square matrix `K`.


### <span class="section-num">0.5</span> General ODE Model Function {#general-ode-model-function}


#### <span class="section-num">0.5.1</span> Description {#description}

Function `pmx_solve_adams`, `pmx_solve_bdf`, and `pmx_solve_rk45` solve a first-order ODE system
specified by user-specified right-hand-side function `ODE_rhs` \\(f\\)

\begin{equation\*}
y'(t) = f(t, y(t))
\end{equation\*}

In the case where the `rate` vector \\(r\\) is non-zero, this equation becomes:

\begin{equation\*}
y'(t) = f(t, y(t)) + r
\end{equation\*}


#### <span class="section-num">0.5.2</span> Usage {#usage}

```stan
matrix pmx_solve_[adams || rk45 || bdf](ODE_rhs, int nCmt, time, amt, rate, ii, evid, cmt, addl, ss, theta, [ biovar, tlag, real[,] x_r, int [,] x_i, real rel_tol, real abs_tol, int max_step, real as_rel_tol, real as_abs_tol, int as_max_step ] );
```

Here `[adams || rk45 || bdf]` indicates the
function name can use any of the three suffixes. See below.


#### <span class="section-num">0.5.3</span> Arguments {#arguments}

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


#### <span class="section-num">0.5.4</span> Return value {#return-value}

An `nCmt`-by-`nt` matrix, where `nt` is the size of `time`.


#### <span class="section-num">0.5.5</span> Note {#note}

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


### <span class="section-num">0.6</span> Coupled ODE Model Function {#coupled-ode-model-function}



#### <span class="section-num">0.6.1</span> Description {#description}

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


#### <span class="section-num">0.6.2</span> Usage {#usage}

```stan
matrix pmx_solve_onecpt_[ rk45 || bdf ](reduced_ODE_system, int nOde, time, amt, rate, ii, evid, cmt, addl, ss, theta, biovar, tlag [, real rel_tol, real abs_tol, int max_step, real as_rel_tol, real as_abs_tol, int as_max_step ] );
matrix pmx_solve_twocpt_[ rk45 || bdf ](reduced_ODE_system, int nOde, time, amt, rate, ii, evid, cmt, addl, ss, theta, biovar, tlag [, real rel_tol, real abs_tol, int max_step, real as_rel_tol, real as_abs_tol, int as_max_step ] );
```


#### <span class="section-num">0.6.3</span> Arguments {#arguments}

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


#### <span class="section-num">0.6.4</span> Return value {#return-value}

An `(nPk + nOde)` &times; `nt` matrix, where `nt` is the size of
`time`, and `nPk` equals to 2 in
`pmx_solve_onecpt_` functions
and 3 in `pmx_solve_twocpt_` functions.


### <span class="section-num">0.7</span> General ODE-based Population Model Function {#general-ode-based-population-model-function}


#### <span class="section-num">0.7.1</span> Description {#description}

All the preivous functions solves for a single sunject. Torsten also
provides population modeling counterparts for ODE solutions. The
functions solve for a population that share an ODE model but with
subject-level parameters and event specifications and have similar
signatures to single-subject functions, except that now events
arguments `time`, `amt`, `rate`, `ii`,
`evid`, `cmt`,
`addl`, `ss` specifies the entire
population, one subject after another.


#### <span class="section-num">0.7.2</span> Usage {#usage}

```stan
matrix pmx_solve_group_[adams || rk45 || bdf](ODE_rhs, int nCmt, int[] len, time, amt, rate, ii, evid, cmt, addl, ss, theta, [ biovar, tlag, real[,] x_r, int [,] x_i, real rel_tol, real abs_tol, int max_step, real as_rel_tol, real as_abs_tol, int as_max_step ] );
```

Here `[adams || rk45 || bdf]` indicates the
function name can be of any of the three suffixes. See Section [ODE  integrator function]({{< relref "ode-integ" >}}).


#### <span class="section-num">0.7.3</span> Arguments {#arguments}

-   `ODE_rhs`
    Same as in Section [ODE  integrator function]({{< relref "ode-integ" >}}).
-   `time`, `amt`, `rate`, `ii`, `evid`, `cmt`, `addl`, `ss`
    2d-array arguments that describe data record for the
    entire population (see also Tables in Section [Events specification]({{< relref "events" >}})). They must have same size in the first
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


#### <span class="section-num">0.7.4</span> Return value {#return-value}

An `nCmt`-by-`nt` matrix, where `nt` is the total size of
events \\(\sum\_{i}n\_i\\).


#### <span class="section-num">0.7.5</span> Note {#sec:ode_pop_note}

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


### <span class="section-num">0.8</span> ODE  integrator function {#ode-integrator-function}


#### <span class="section-num">0.8.1</span> Description {#description}

Torsten provides its own implementation of ODE solvers that solves

\begin{equation\*}
  y'(t) = f(t, y(t)), \quad y(t\_0) = y\_0
\end{equation\*}

for \\(y\\). These solvers
are customized for Torsten applications and different from those found
in Stan. The general ODE PMX solvers in previous sections are internally powered
by these functions.


#### <span class="section-num">0.8.2</span> Usage {#usage}

```stan
real[ , ] pmx_integrate_ode_[ adams || bdf || rk45 ](ODE_rhs, real[] y0, real t0, real[] ts, real[] theta, real[] x_r, int[] x_i [ , real rtol, real atol, int max_step ]);
```


#### <span class="section-num">0.8.3</span> Arguments {#arguments}

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


#### <span class="section-num">0.8.4</span> Return value {#return-value}

An `n`-by-`nd` 2d-array, where `n` is the size of `ts`
and `nd` the dimension of the system.


#### <span class="section-num">0.8.5</span> Note {#note}

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
    ([Hindmarsh et al. 2020](#org8935a5e)) and
    1.4 of ([Shampine et al. 2003](#orgd3b17e4)) for details.

-   With optional arguments indicated by square bracket, the following calls are allowed:

<!--listend-->

```stan
pmx_integrate_ode_[adams || rk45 || bdf](..., x_i);
pmx_integrate_ode_[adams || rk45 || bdf](..., x_i, rel_tol, abs_tol, max_step);
```


### <span class="section-num">0.9</span> ODE group  integrator Function {#ode-group-integrator-function}



#### <span class="section-num">0.9.1</span> Description {#description}

All the preivous functions solves for a single ODE system. Torsten also
provides group modeling counterparts for ODE integrators. The
functions solve for a group of ODE systems that share an ODE RHS but with
different parameters. They have similar
signatures to single-ODE integration functions.


#### <span class="section-num">0.9.2</span> Usage {#usage}

```stan
matrix pmx_integrate_ode_group_[adams || rk45 || bdf](ODE_system, real[ , ] y0, real t0, int[] len, real[] ts, real[ , ] theta, real[ , ] x_r, int[ , ] x_i, [ real rtol, real atol, int max_step ] );
```

Here `[adams || rk45 || bdf]` indicates the
function name can be of any of the three suffixes. See Section [ODE  integrator function]({{< relref "ode-integ" >}}).

-   `ODE_rhs`
    Function that specifies the right-hand-side \\(f\\). See Section [ODE  integrator function]({{< relref "ode-integ" >}}).
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


#### <span class="section-num">0.9.3</span> Return value {#return-value}

An `n`-by-`nd` matrix, where `n` is the size of `ts`
and `nd` the dimension of the system.


#### <span class="section-num">0.9.4</span> Note {#note}

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


### <span class="section-num">0.10</span> Univariate integral {#univariate-integral}

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


### <span class="section-num">0.11</span> Piecewise linear interpolation {#piecewise-linear-interpolation}

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


## <span class="section-num">1</span> Bibliography {#bibliography}

<a id="org8935a5e"></a>Hindmarsh, Alan C., Radu Serban, Cody J. Balos, David J. Gardner, Carol S. Woodward, and Daniel R. Reynolds. 2020. _User Documentation for Cvodes V5.4.0_.

<a id="orgd3b17e4"></a>Shampine, L. F., I. Gladwell, Larry Shampine, and S. Thompson. 2003. _Solving ODEs with MATLAB_. Cambridge University Press.
