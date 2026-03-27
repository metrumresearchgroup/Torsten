#ifndef STAN_MATH_TORSTEN_IS_STD_ODE_HPP
#define STAN_MATH_TORSTEN_IS_STD_ODE_HPP

#include <stdexcept>
#include <vector>

namespace torsten {
  /** 
   * Whether a Stan user-defined function can be used as RHS of ode solver
   * The solver is supposed to accept "old" Stan ode integrator signature.
   * 
   */
  template <typename F>
  struct is_std_ode {
    using yes = double;
    using no = bool;

    template <typename Functor>
    static double test(decltype(std::declval<Functor&>()(0.0, std::vector<double>(),
                                                         std::vector<double>(),
                                                         std::vector<double>(),
                                                         std::vector<int>(), nullptr))*);

    template <typename Functor>
    static no test(...);

    static constexpr bool value = sizeof(test<F>(nullptr)) == sizeof(yes);
  };

  template <typename F>
  struct std_ode {
    using type = std::enable_if_t<torsten::is_std_ode<F>::value>;
  };
}

#endif
