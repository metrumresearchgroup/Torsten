#ifndef STAN_MATH_TORSTEN_PKMODEL_UNPROMOTE_HPP
#define STAN_MATH_TORSTEN_PKMODEL_UNPROMOTE_HPP

#include <iostream>
#include <stan/math/rev/core.hpp>
#include <stan/math/fwd/core.hpp>

// DEV: Not a big fan of having to use this function. 

/**
 * Functions that converts an autodiff variable into a double. 
 * The variable will either be a stan::math::var
 * or a double (in which case it will not be modified).
 */
inline double unpromote(stan::math::var x) { return x.val(); }
inline double unpromote(double x) { return x; }

#endif
