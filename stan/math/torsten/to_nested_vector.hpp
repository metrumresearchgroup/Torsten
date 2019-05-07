#ifndef STAN_MATH_TORSTEN_TO_NESTED_VECTOR_HPP
#define STAN_MATH_TORSTEN_TO_NESTED_VECTOR_HPP

#include <stan/math/torsten/is_std_vector.hpp>

namespace torsten {
  template<typename T,
           typename std::enable_if_t<!torsten::is_std_vector<T>::value >* = nullptr> //NOLINT
  const std::vector<std::vector<T>> to_nested_vector(const std::vector<T>& v)
  {
    return {v};
  }

  template<typename T>
  const std::vector<std::vector<T>>& to_nested_vector(const std::vector<std::vector<T>>& v)
  {
    return v;
  }
}

#endif
