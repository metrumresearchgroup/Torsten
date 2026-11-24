#ifndef STAN_MATH_TORSTEN_DEF_HPP
#define STAN_MATH_TORSTEN_DEF_HPP

#include <Eigen/Dense>
#include <stan/math/prim/meta/return_type.hpp>
#include <stan/math/prim/meta/is_var.hpp>
#include <stan/math/rev/meta/is_var.hpp>

namespace torsten {

  template<typename T>
  using PKRec = Eigen::Matrix<T, -1, 1>;

  template<typename T>
  using PMXLin = Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>;

  // template <class T>
  // using rate_t = typename std::decay<decltype(std::declval<T>().rate())>::type::value_type; // NOLINT

  // template <class T>
  // using init_t = typename std::decay<decltype(std::declval<T>().y0())>::type::value_type; // NOLINT

  template <class T>
  using par_t = typename std::decay<decltype(std::declval<T>().par())>::type::value_type; // NOLINT

  // template <class T>
  // using aug_par_t = typename stan::return_type<rate_t<T>, par_t<T>, init_t<T>>::type; // NOLINT

  template <class T>
  using f_t = typename std::decay<decltype(std::declval<T>().f())>::type; // NOLINT
}

#endif
