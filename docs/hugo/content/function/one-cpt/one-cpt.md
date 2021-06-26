+++
title = "One Compartment Model"
author = ["Yi Zhang"]
date = 2021-06-25T00:00:00-07:00
lastmod = 2021-06-26T09:38:12-07:00
draft = false
weight = 2001
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

Function `pmx_solve_onecpt` solves a one-compartment PK
model (Figure [1](#org92aa0a2)). The model obtains plasma concentrations of parent drug \\(c=y\_2/V\_2\\)
by solving for the mass of drug in the central compartment
\\(y\_2\\) from ordinary differential equations(ODEs)

\begin{align}\label{eq:onecpt}
  y\_1' &= -k\_a y\_1, \\\\\\
  y\_2' &= k\_a y\_1 - \left(\frac{CL}{V\_2} + \frac{Q}{V\_2}\right) y\_2.
\end{align}

<a id="org92aa0a2"></a>

{{< figure src="/ox-hugo/cptModels.png" caption="Figure 1: One and two compartment models with first order absorption implemented in Torsten." >}}


## <span class="section-num">2</span> Usage {#usage}

```stan
matrix = pmx_solve_onecpt(time, amt, rate, ii, evid, cmt, addl, ss, theta [, biovar, tlag ] )
```


## <span class="section-num">3</span> Arguments {#arguments}

See Table [tab:event_args](#tab:event_args) and Table [tab:event_params](#tab:event_params).


## <span class="section-num">4</span> Return value {#return-value}

An `ncmt`-by-`nt` matrix, where `nt` is the number of time steps and `ncmt=2` is the number of compartments.


## <span class="section-num">5</span> Note {#note}

-   ODE Parameters `theta` should consist of \\(CL\\), \\(V\_2\\), \\(k\_a\\), in that order.
-   `biovar` and `tlag` are optional, so that the following are allowed:

<!--listend-->

```stan
pmx_solve_onecpt(..., theta);
pmx_solve_onecpt(..., theta, biovar);
pmx_solve_onecpt(..., theta, biovar, tlag);
```

-   Setting \\(k\_a = 0\\) eliminates the first-order absorption.
