#ifndef STAN_MATH_TORSTEN_IS_STD_VECTOR_HPP
#define STAN_MATH_TORSTEN_IS_STD_VECTOR_HPP

#include <type_traits>
#include <vector>

namespace torsten {
  template<typename T1, typename... Tn>
  struct is_std_vector : std::false_type {};

  template<typename T, typename... Ts>
  struct is_std_vector<std::vector<T, Ts...>> : std::true_type {};

  template<typename T, typename... Ts, typename... Tn>
  struct is_std_vector<std::vector<T, Ts...>, Tn...> : torsten::is_std_vector<Tn...> {};

  template<typename T>
  struct value_type {
    using type = T;
  };

  template<typename T>
  struct value_type<std::vector<T> > {
    using type = T;
  };

}

#endif
