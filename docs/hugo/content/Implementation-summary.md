+++
title = "Implementation summary"
author = ["Yi Zhang"]
date = 2021-06-25T00:00:00-07:00
lastmod = 2021-06-28T14:19:25-07:00
draft = false
weight = 2002
+++

-   Current Torsten v0.89rc is based on Stan v2.27.0.
-   All functions are programmed in C++ and are compatible
    with the Stan math automatic differentiation library ([Carpenter et al. 2015](#org21d88de))
-   One and two compartment models are based on analytical solutions of governing ODEs.
-   General linear compartment models are based on semi-analytical solutions using the built-in matrix exponential function
-   General compartment models are solved numerically using built-in ODE integrators in Stan. The tuning parameters of the solver are adjustable. The steady state solution is calculated using a numerical algebraic solver.
-   A mix compartment model's PK forcing function is solved analytically, and its forced ODE system is solved numerically.


## <span class="section-num">1</span> Bibliography {#bibliography}

<a id="org21d88de"></a>Carpenter, Bob, Matthew D. Hoffman, Marcus Brubaker, Daniel Lee, Peter Li, and Michael Betancourt. 2015. “The Stan Math Library: Reverse-Mode Automatic Differentiation in C++.” _arXiv:1509.07164 [Cs]_, September. <http://arxiv.org/abs/1509.07164>.
