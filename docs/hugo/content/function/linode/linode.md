+++
title = "General Linear ODE Model Function"
author = ["Yi Zhang"]
date = 2021-06-25T00:00:00-07:00
lastmod = 2021-06-26T09:38:13-07:00
draft = false
weight = 2003
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


## <span class="section-num">2</span> Usage {#usage}

```stan
matrix = pmx_solve_linode(time, amt, rate, ii, evid, cmt, addl, ss, K, biovar, tlag )
```


## <span class="section-num">3</span> Arguments {#arguments}

-   `K`
    System parameters. `K` can be either
    -   a `matrix` for constant parameters in all events, or
    -   an array of matrices `matrix K[nt]` so that the \\(i\\)th entry of the array describes
        the model parameters for time interval \\((t\_{i-1}, t\_i)\\),
        and the number of the rows equals to the number of event time `nt`.
-   See Table [tab:event_args](#tab:event_args) and Table [tab:event_params](#tab:event_params) for the rest of arguments.


## <span class="section-num">4</span> Return value {#return-value}

An `n`-by-`nt` matrix, where `nt` is the number of time steps and `n` is the number of rows(columns) of square matrix `K`.
