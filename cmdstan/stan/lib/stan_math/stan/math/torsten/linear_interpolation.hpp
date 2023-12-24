#ifndef STAN_MATH_TORSTEN_LINEAR_INTERPOLATION_HPP
#define STAN_MATH_TORSTEN_LINEAR_INTERPOLATION_HPP

#include <stan/math/rev.hpp>
#include <stan/math/prim.hpp>
#include <stan/math/prim/err/check_nonzero_size.hpp>
#include <stan/math/torsten/PKModel/SearchReal.hpp>
#include <stan/math/prim/meta/return_type.hpp>
#include <vector>
#include <algorithm>

namespace torsten {

  /**
   * Return the values of a piecewise linear function at specifed values of the 
   * function argument. The function is specified in terms of a set of x y pairs.
   *
   * @tparam T0 Type of value to be interpolated
   * @tparam T1 Type of elements contained in x array.
   * @tparam T2 Type of elements contained in y array.
   * @param xout Scalar or array of function argument values.
   * @param x Array of x values. Must be in increasing order.
   * @param y Array of y values. Must have same length as x.
   * @return Scalar or array of function values.
   */

  template <typename T0, typename T1, typename T2>
  typename stan::return_type_t<T0, T1, T2>
  inline linear_interpolation(const T0& xout,
                              const std::vector<T1>& x,
                              const std::vector<T2>& y) {
    stan::math::check_finite("linear_interpolation", "xout", xout);
    stan::math::check_finite("linear_interpolation", "x", x);
    stan::math::check_finite("linear_interpolation", "y", y);
    stan::math::check_nonzero_size("linear_interpolation", "x", x);
    stan::math::check_nonzero_size("linear_interpolation", "y", y);
    stan::math::check_ordered("linear_interpolation", "x", x);
    stan::math::check_matching_sizes("linear_interpolation", "x", x, "y", y);

    typename stan::return_type_t<T0, T1, T2> yout;
    if (xout < x.front()) {
      yout = y.front();
    } else if (xout > x.back()) {
      yout = y.back();
    } else {
      int i = std::lower_bound(x.begin(), x.end(), xout) - x.begin() - 1;
      yout = y.at(i) + (y.at(i+1) - y.at(i)) / (x.at(i+1) - x.at(i)) * (xout - x.at(i));
    }
    return yout;
  }

  template <typename T0, typename T1, typename T2>
  std::vector<typename stan::return_type_t<T0, T1, T2>>
  inline linear_interpolation(const std::vector<T0>& xout,
                              const std::vector<T1>& x,
                              const std::vector<T2>& y) {
    stan::math::check_nonzero_size("linear_interpolation", "xout", xout);
    stan::math::check_finite("linear_interpolation", "xout", xout);

    std::vector<typename stan::return_type_t<T0, T1, T2>> yout(xout.size());
    std::transform(xout.begin(), xout.end(), yout.begin(),
                   [&x, &y](const T0& xi) { return linear_interpolation(xi, x, y); });
    return yout;
  }
}
#endif
