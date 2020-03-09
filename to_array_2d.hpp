#ifndef STAN_MATH_TORSTEN_TO_ARRAY_2D_HPP
#define STAN_MATH_TORSTEN_TO_ARRAY_2D_HPP

#include <stan/math/torsten/is_std_vector.hpp>

namespace torsten {
  template<typename T,
           typename std::enable_if_t<!torsten::is_std_vector<T>::value >* = nullptr> //NOLINT
  const std::vector<std::vector<T>> to_array_2d(const std::vector<T>& v)
  {
    return {v};
  }

  template<typename T>
  const std::vector<std::vector<T>>  to_array_2d(const T* p, int m, int n)
  {
    std::vector<std::vector<T>> res(m, std::vector<T>(n));
    for (int i = 0, k = 0; i < m; ++i) {
      for (int j = 0; j < n; ++j, ++k) {
        res[i][j] = *(p + k);
      }
    }
    return res;
  }

  template<typename T>
  const std::vector<std::vector<T>>& to_array_2d(const std::vector<std::vector<T>>& v)
  {
    return v;
  }
}

#endif
