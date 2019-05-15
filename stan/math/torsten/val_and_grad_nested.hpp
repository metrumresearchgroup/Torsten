#ifndef STAN_MATH_TORSTEN_VAL_AND_GRAD_HPP
#define STAN_MATH_TORSTEN_VAL_AND_GRAD_HPP

#include <stan/math/rev/core.hpp>
#include <stan/math/rev/mat/fun/typedefs.hpp>

namespace torsten {

  inline Eigen::VectorXd
  val_and_grad_nested(stan::math::vector_v& yv, const std::vector<stan::math::var>& theta) {
    const int n = yv.size();
    const int m = theta.size();
    Eigen::VectorXd yd( n + n * m );
    for (int i = 0; i < n; ++i) {
      yd[i] = yv[i].val();
      stan::math::set_zero_all_adjoints_nested();
      yv[i].grad();
      for (int j = 0; j < m; ++j) {
        yd[n + j * n + i] = theta[j].adj();
      }
    }
    return yd;
  }

}
#endif
