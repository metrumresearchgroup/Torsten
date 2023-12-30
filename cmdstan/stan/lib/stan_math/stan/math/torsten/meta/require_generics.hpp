#ifndef STAN_MATH_TORSTEN_META_REQUIRE_GENERICS_HPP
#define STAN_MATH_TORSTEN_META_REQUIRE_GENERICS_HPP

#include <stan/math/prim/meta/require_helpers.hpp>
#include <stan/math/prim/meta/is_vector.hpp>
#include <stan/math/torsten/meta/is_std_ode.hpp>
#include <stan/math/torsten/meta/is_eigen_ode.hpp>

namespace torsten {

  using stan::require_t;
  using stan::require_not_t;
  using stan::require_all_t;
  using stan::require_any_t;
  using stan::require_any_not_t;
  using stan::require_all_not_t;

  STAN_ADD_REQUIRE_UNARY(std_ode, torsten::is_std_ode, torsten_meta)
  STAN_ADD_REQUIRE_UNARY(eigen_ode, torsten::is_eigen_ode, torsten_meta)
  STAN_ADD_REQUIRE_UNARY(std_vector, stan::is_std_vector, torsten_meta)
}

#endif
