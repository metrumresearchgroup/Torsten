+++
title = "Two Compartment Model"
author = ["Yi Zhang"]
date = 2021-06-25T00:00:00-07:00
lastmod = 2021-06-28T19:35:14-07:00
draft = false
weight = 2002
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

<a id="org1030393"></a>


## <span class="section-num">1</span> Description {#description}

Function `pmx_solve_twocpt` solves a two-compartment PK
model (Figure [{{< relref "one-cpt" >}}]({{< relref "one-cpt" >}})). The model obtains plasma concentrations of parent drug \\(c=y\_2/V\_2\\)
by solving for the mass of drug in the central compartment
\\(y\_2\\) from ordinary differential equations(ODEs)

\begin{align} \label{eq:twocpt}
  y\_1' &= -k\_a y\_1 \\\\\\
  y\_2' &= k\_a y\_1 - \left(\frac{CL}{V\_2} + \frac{Q}{V\_2}\right) y\_2 +  \frac{Q}{V\_3}  y\_3  \\\\\\
  y\_3' &= \frac{Q}{V\_2} y\_2 - \frac{Q}{V\_3} y\_3
\end{align}


## <span class="section-num">2</span> Usage {#usage}

```stan
matrix = pmx_solve_twocpt(time, amt, rate, ii, evid, cmt, addl, ss, theta [, biovar, tlag ] )
```


## <span class="section-num">3</span> Arguments {#arguments}

See Table [tab:event_args](#tab:event_args) and Table [tab:event_params](#tab:event_params).


## <span class="section-num">4</span> Return value {#return-value}

An `ncmt`-by-`nt` matrix, where `nt` is the number of time steps and `ncmt=3` is the number of compartments.


## <span class="section-num">5</span> Note {#note}

-   ODE Parameters `theta` consists of \\(CL\\), \\(Q\\), \\(V\_2\\), \\(V\_3\\), \\(k\_a\\).
-   `biovar` and `tlag` are optional, so that the following are allowed:

<!--listend-->

```stan
pmx_solve_twocpt(..., theta);
pmx_solve_twocpt(..., theta, biovar);
pmx_solve_twocpt(..., theta, biovar, tlag);
```

-   Setting \\(k\_a = 0\\) eliminates the first-order absorption.
