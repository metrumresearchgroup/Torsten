#ifndef STAN_MATH_TORSTEN_DSOLVE_ODE_FUNC_HPP
#define STAN_MATH_TORSTEN_DSOLVE_ODE_FUNC_HPP

#include <stan/math/torsten/cvodes_sens_method.hpp>

namespace torsten {
  namespace dsolve {
    /*
     * Retrieve ODE template RHS functor type
     */
    template <typename Ode>
    struct ode_func;

    /*
     * Retrieve ODE template RHS functor type
     */
    template<template<typename> class Ode, typename F>
    struct ode_func<Ode<F>> {
      using type = F;
    };

    /*
     * Retrieve ODE template RHS functor type
     */
    template<template<typename, typename...> class Ode, typename F, typename... Ts>
    struct ode_func<Ode<F, Ts...>> {
      using type = F;
    };
  }
}

#endif
