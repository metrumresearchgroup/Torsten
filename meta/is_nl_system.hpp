#ifndef STAN_MATH_TORSTEN_META_IS_NL_SYSTEM_HPP
#define STAN_MATH_TORSTEN_META_IS_NL_SYSTEM_HPP

#include <stan/math/prim/fun/Eigen.hpp>
#include <stan/math/prim/meta/return_type.hpp>
#include <stdexcept>
#include <vector>

namespace torsten {
  
  /** 
   * Whether a Stan user-defined function can be used as nonlinear
   * system of which is the zero is to be seeked
   * The solver is supposed to accept variadic params.
   * 
   */
  template <typename F, typename... Args>
  struct is_nl_system {
    using yes = double;
    using no = bool;

    template <typename Functor>
    static double test(decltype(std::declval<Functor&>()(Eigen::VectorXd(), nullptr,
                                                         std::declval<Args>()...))*);

    template <typename Functor>
    static no test(...);

    static constexpr bool value = sizeof(test<F>(nullptr)) == sizeof(yes);
  };


  /** 
   * Adaptor for existing Stan nl system signature.
   * 
   */
  template<typename F>
  struct nl_system_adaptor {
    F const& f_;

    nl_system_adaptor(F const& f) : f_(f) {}

    template <typename T0, typename T1>
    Eigen::Matrix<stan::return_type_t<T0, T1>, -1, 1>
    operator()(const Eigen::Matrix<T0, Eigen::Dynamic, 1>& x,
               std::ostream* pstream__,
               const Eigen::Matrix<T1, Eigen::Dynamic, 1>& y,
               const std::vector<double>& dat, const std::vector<int>& dat_int) const {
      return f_(x, y, dat, dat_int, pstream__);
    }
  };
}

#endif
