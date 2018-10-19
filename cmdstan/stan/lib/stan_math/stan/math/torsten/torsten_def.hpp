#ifndef TORSTEN_DEF_HPP
#define TORSTEN_DEF_HPP

#include <stan/math/rev/core.hpp>
#include <stan/math/torsten/is_detected.hpp>

namespace refactor {
    template<typename T>
    using PKRec = Eigen::Matrix<T, 1, Eigen::Dynamic>;

    template<typename T>
    using PKLin = Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>;

}

namespace torsten {

  template <class T>
  using rate_t = typename std::decay<decltype(std::declval<T>().rate())>::type::value_type; // NOLINT

  template <class T>
  using init_t = typename std::decay<decltype(std::declval<T>().y0())>::type::value_type; // NOLINT

  template <class T>
  using par_t = typename std::decay<decltype(std::declval<T>().par())>::type::value_type; // NOLINT

  template <class T>
  using detime_t = typename std::decay<decltype(std::declval<T>().t0())>::type; // NOLINT

  template <class T>
  using scalar_t = typename stan::return_type<detime_t<T>, rate_t<T>, par_t<T>, init_t<T>>::type; // NOLINT

  template <class T>
  using aug_par_t = typename stan::return_type<rate_t<T>, par_t<T>, init_t<T>>::type; // NOLINT

  template <class T>
  using f_t = typename std::decay<decltype(std::declval<T>().f())>::type; // NOLINT

  template <class T>
  using has_var_rate = stan::is_var<rate_t<T> >;

  template <class T>
  using has_var_init = stan::is_var<init_t<T> >;

  template <class T>
  using has_var_par = stan::is_var<par_t<T> >;

}

#endif
