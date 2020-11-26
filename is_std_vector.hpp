#ifndef STAN_MATH_TORSTEN_IS_STD_VECTOR_HPP
#define STAN_MATH_TORSTEN_IS_STD_VECTOR_HPP

#include <type_traits>
#include <vector>

namespace torsten {
  /**
   * type trait for a parameter pack with types that are not all
   * <code>std::vector</code>
   * 
   */
  template<typename T1, typename... Tn>
  struct is_std_vector : std::false_type {};

  template<typename T, typename... Ts>
  struct is_std_vector<std::vector<T, Ts...>> : std::true_type {};

  template<typename T, typename... Ts, typename... Tn>
  struct is_std_vector<std::vector<T, Ts...>, Tn...> : torsten::is_std_vector<Tn...> {};

  /**
   * type trait for a parameter pack with types none of which is
   * <code>std::vector</code>
   * 
   */
  template<typename... Ts>
  struct none_std_vector;

  template<typename T>
  struct none_std_vector<T> : std::true_type {};

  template<typename T, typename... Ts>
  struct none_std_vector<std::vector<T, Ts...>> : std::false_type {};

  template<typename T, typename... Tn>
  struct none_std_vector<T, Tn...> :
    std::conditional<none_std_vector<T>::value, none_std_vector<Tn...>, std::false_type>::type {};


  // value type of std vector
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
