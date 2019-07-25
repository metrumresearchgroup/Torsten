#ifndef STAN_MATH_TORSTEN_RETURN_TYPE_HPP
#define STAN_MATH_TORSTEN_RETURN_TYPE_HPP

#include <stan/math/rev/core/var.hpp>
#include <stan/math/prim/scal/meta/return_type.hpp>

namespace torsten {
  template<typename... Ts>
  struct return_t;

  template<typename T>
  struct return_t<T>
  {
    using type = T;
  };

  template<typename T1, typename... Ts>
  struct return_t<T1, Ts...>
  {
    using type = typename stan::return_type<T1, typename torsten::return_t<Ts...>::type>::type; //NOLINT
  };
}

#endif
