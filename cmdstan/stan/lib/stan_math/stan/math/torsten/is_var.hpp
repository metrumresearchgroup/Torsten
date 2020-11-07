#ifndef STAN_MATH_TORSTEN_IS_VAR_HPP
#define STAN_MATH_TORSTEN_IS_VAR_HPP

#include <stan/math/prim/meta/is_var.hpp>
#include <stan/math/rev/core.hpp>
#include <stan/math/prim/meta/return_type.hpp>

namespace torsten {
  
  /*
   * sugar template to simplify @c is_var check
   */
  template<typename... Ts>
  struct is_var {
    enum { value = stan::is_var<typename stan::return_type_t<Ts...>>::value };
  };
}

#endif
