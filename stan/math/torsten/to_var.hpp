#ifndef STAN_MATH_TORSTEN_TO_VAR
#define STAN_MATH_TORSTEN_TO_VAR

#include <stan/math/rev/mat/fun/typedefs.hpp>
#include <stan/math/prim/arr/fun/value_of.hpp>
#include <stan/math/rev/core/var.hpp>
#include <iostream>
#include <sstream>
#include <vector>
#include <string>

namespace torsten {
  template<typename T>
  std::vector<std::vector<stan::math::var> > to_var(const std::vector<std::vector<T> >& d)
  {
    std::vector<std::vector<stan::math::var> > res(d.size());
    for (size_t i = 0; i < d.size(); ++i) {
      res[i].resize(d[i].size());
      for (size_t j = 0; j < d[i].size(); ++j) {
        res[i][j] = stan::math::value_of(d[i][j]);
      }
    }
    return res;
  }
    
  template<typename T>
  std::vector<stan::math::var> to_var(const std::vector<T>& d)
  {
    std::vector<stan::math::var> res(d.size());
    for (size_t i = 0; i < d.size(); ++i) {
      res[i] = stan::math::value_of(d[i]);
    }
    return res;
  }

  template<typename T>
  std::vector<stan::math::matrix_v> to_var(const std::vector<Eigen::Matrix<T, -1, -1> >& d)
  {
    std::vector<stan::math::matrix_v> res(d.size());
    for (size_t i = 0; i < d.size(); ++i) {
      res[i] = stan::math::value_of(d[i]);
    }
    return res;
  }
}

#endif
