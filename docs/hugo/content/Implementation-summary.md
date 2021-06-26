+++
title = "Implementation summary"
author = ["Yi Zhang"]
date = 2021-06-25T00:00:00-07:00
lastmod = 2021-06-26T09:38:07-07:00
draft = false
weight = 2002
+++

-   Current Torsten v0.89rc is based on Stan v2.27.0.
-   All functions are programmed in C++ and are compatible
    with the Stan math automatic differentiation library <sup id="5edcf88ae76ae7132583b12fffe35055"><a href="#carpenter15_stan_math_librar" title="Carpenter, Hoffman, , Brubaker, Lee, Li, Peter \&amp; Betancourt, The {Stan} {Math} {Library}:  {Reverse}-{Mode} {Automatic}  {Differentiation} in {C}++, {arXiv:1509.07164 [cs]}, v(), (2015).">carpenter15_stan_math_librar</a></sup>
-   One and two compartment models are based on analytical solutions of governing ODEs.
-   General linear compartment models are based on semi-analytical solutions using the built-in matrix exponential function
-   General compartment models are solved numerically using built-in ODE integrators in Stan. The tuning parameters of the solver are adjustable. The steady state solution is calculated using a numerical algebraic solver.
-   A mix compartment model's PK forcing function is solved analytically, and its forced ODE system is solved numerically.
