+++
title = "Development plans"
author = ["Yi Zhang"]
date = 2021-06-25T00:00:00-07:00
lastmod = 2021-06-29T14:35:19-07:00
draft = false
weight = 3002
+++

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
