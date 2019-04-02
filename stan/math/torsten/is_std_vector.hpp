#ifndef STAN_MATH_TORSTEN_IS_STD_VECTOR_HPP
#define STAN_MATH_TORSTEN_IS_STD_VECTOR_HPP

#include <type_traits>

namespace torsten {
  template<typename>
  struct is_std_vector : std::false_type {};

  template<typename T, typename... Ts>
  struct is_std_vector<std::vector<T, Ts...>> : std::true_type {};
}

#endif
