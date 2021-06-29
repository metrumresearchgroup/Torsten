+++
title = "Piecewise linear interpolation"
author = ["Yi Zhang"]
date = 2021-06-25T00:00:00-07:00
lastmod = 2021-06-29T11:51:50-07:00
draft = false
weight = 2011
+++

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
