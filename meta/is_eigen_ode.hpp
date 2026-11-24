#ifndef STAN_MATH_TORSTEN_META_IS_EIGEN_ODE_HPP
#define STAN_MATH_TORSTEN_META_IS_EIGEN_ODE_HPP

#include <stan/math/prim/fun/Eigen.hpp>
#include <stdexcept>
#include <vector>

namespace torsten {
  
  /** 
   * Whether a Stan user-defined function can be used as RHS of ode solver
   * The solver is supposed to accept "new" Stan ode integrator
   * signature that supports variadic params
   * 
   */
  template <typename F, typename... Args>
  struct is_eigen_ode {
    using yes = double;
    using no = bool;

    template <typename Functor>
    static double test(decltype(std::declval<Functor&>()(0.0, Eigen::VectorXd(), nullptr,
                                                         std::declval<Args>()...))*);

    template <typename Functor>
    static no test(...);

    static constexpr bool value = sizeof(test<F>(nullptr)) == sizeof(yes);
  };

  template <typename F, typename... Args>
  struct eigen_ode {
    using type = std::enable_if_t<torsten::is_eigen_ode<F, Args...>::value>;
  };
}

#endif
