+++
title = "Introduction"
author = ["Yi Zhang"]
date = 2021-06-25T00:00:00-07:00
lastmod = 2021-06-29T14:35:18-07:00
draft = false
weight = 1003
+++

[Stan](https://mc-stan.org/) is an open source probabilistic programing language designed
primarily to do Bayesian data analysis
([Carpenter et al. 2017](#org4e1e0b2)). It provides an expressive syntax for statistic
modeling and contains an efficient variant of No U-Turn
Sampler(NUTS), an adaptative Hamiltonian Monte Carlo
algorithm that was proven more efficient than commonly used Monte Carlo Markov Chains
(MCMC) samplers for complex high dimensional problems ([Hoffman and Gelman 2011](#orgdb7f16c); [Betancourt 2018](#orgaa94fbb)).

{{< figure src="https://metrumrg.com/wp-content/uploads/2019/07/torsten-white-stan-cropped.png" width="300" >}}

Torsten is a collection of Stan functions to facilitate analysis of
pharmacometric data. Given an event schedule and an ODE system, it calculates amounts
in each compartment. The current version includes&nbsp;[^fn:1]:

-   Specific linear compartment models:
    -   One compartment model with first order absorption.
    -   Two compartment model with elimination from and first order absorption into central compartment
-   General linear compartment model described by a system of first-order <span class="underline">linear</span> Ordinary Differential Equations (ODEs).
-   General compartment model described by a system of first order ODEs.
-   Coupled model with PK forcing function described by a linear one or two compartment model and PD components solved by numerical ODE integration.

The models and data format are based on
NONMEM \textregistered{}&nbsp;[^fn:2]/NMTRAN/PREDPP
conventions including:

-   recursive calculation of model predictions, which permits piecewise constant covariate values,
-   bolus or constant rate inputs into any compartment,
-   single dose and multiple dose events,
-   steady state dosing events,
-   NMTRAN-compartible data items such as TIME, EVID, CMT, AMT, RATE, ADDL, II, and SS.

All real variable arguments in Torsten functions can be passed as Stan `parameters`.


#### Implementation summary {#implementation-summary}

-   Current Torsten v0.89rc is based on Stan v2.27.0.
-   All functions are programmed in C++ and are compatible
    with the Stan math automatic differentiation library ([Carpenter et al. 2015](#orgcf8b534))
-   One and two compartment models are based on analytical solutions of governing ODEs.
-   General linear compartment models are based on semi-analytical solutions using the built-in matrix exponential function
-   General compartment models are solved numerically using built-in ODE integrators in Stan. The tuning parameters of the solver are adjustable. The steady state solution is calculated using a numerical algebraic solver.
-   Coupled model that has PK forcing function solved analytically and PD ODE components solved numerically.


#### Development plans {#development-plans}

Our current plans for future development of Torsten include the
following:

-   Build a system to easily share packages of Stan functions
    (written in C++ or in the Stan language)
-   Optimize Matrix exponential functions
    -   Function for the action of Matrix Exponential on a vector
    -   Hand-coded gradients
    -   Special algorithm for matrices with special properties
-   Develop new method for large-scale hierarchical models with costly
    ODE solving.
-   Fix issue that arises when computing the adjoint of the lag time
    parameter (in a dosing compartment) evaluated at \\(t\_{\text{lag}} = 0\\).
-   Extend formal tests
    -   More unit tests and better CD/CI support.
    -   Comparison with simulations from the R package
        `mrgsolve` and the software NONMEM\textregistered{}
    -   Recruit non-developer users to conduct beta testing


## Bibliography {#bibliography}

<a id="orgaa94fbb"></a>Betancourt, Michael. 2018. “A Conceptual Introduction to Hamiltonian Monte Carlo.” <https://arxiv.org/abs/1701.02434>.

<a id="org4e1e0b2"></a>Carpenter, Bob, Andrew Gelman, Matthew D. Hoffman, Daniel Lee, Ben Goodrich, Michael Betancourt, Marcus Brubaker, Jiqiang Guo, Peter Li, and Allen Riddell. 2017. “Stan: A Probabilistic Programming Language.” _Journal of Statistical Software_ 76.

<a id="orgcf8b534"></a>Carpenter, Bob, Matthew D. Hoffman, Marcus Brubaker, Daniel Lee, Peter Li, and Michael Betancourt. 2015. “The Stan Math Library: Reverse-Mode Automatic Differentiation in C++.” _arXiv:1509.07164 [Cs]_, September. <http://arxiv.org/abs/1509.07164>.

<a id="orgdb7f16c"></a>Hoffman, Matthew D., and Andrew Gelman. 2011. “The No-U-Turn Sampler: Adaptively Setting Path Lengths in Hamiltonian Monte Carlo.” _arXiv:1111.4246 [Cs, Stat]_, November. <http://arxiv.org/abs/1111.4246>.

[^fn:1]: **WARNING:** The current version of Torsten is a _prototype_. It is being released for review and comment, and to support limited research applications. It has not been rigorously tested and should not be used for critical applications without further testing or cross-checking by comparison with other methods. We encourage interested users to try Torsten out and are happy to assist. Please report issues, bugs, and feature requests on [our GitHub page](https://github.com/metrumresearchgroup/stan).
[^fn:2]: NONMEM\textregistered{} is licensed and distributed by ICON Development Solutions.
