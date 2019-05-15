#ifndef STAN_MATH_TORSTEN_IS_VAR_HPP
#define STAN_MATH_TORSTEN_IS_VAR_HPP

#include <stan/math/prim/scal/meta/is_var.hpp>
#include <stan/math/rev/scal/meta/is_var.hpp>
#include <stan/math/torsten/return_type.hpp>

namespace torsten {
  
  /*
   * sugar template to simplify @c is_var check
   */
  template<typename... Ts>
  struct is_var {
    enum { value = stan::is_var<typename torsten::return_t<Ts...>::type>::value };
    // static const bool value = stan::is_var<typename torsten::return_t<Ts...>::type>::value;
  };
}

#endif
