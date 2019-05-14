#ifndef STAN_MATH_TORSTEN_LINEAR_INTERPOLATION_HPP
#define STAN_MATH_TORSTEN_LINEAR_INTERPOLATION_HPP

#include <stan/math/rev/mat.hpp>
#include <stan/math/prim/arr.hpp>
#include <stan/math/prim/arr/err/check_nonzero_size.hpp>
#include <stan/math/torsten/PKModel/SearchReal.hpp>
#include <vector>
#include <algorithm>

namespace torsten {

  /**
   * Return the values of a piecewise linear function at specifed values of the 
   * function argument. The function is specified in terms of a set of x y pairs.
   *
   * @tparam T0 Type of elements contained in xout array.
   * @tparam T1 Type of elements contained in y array.
   * @param xout Scalar or array of function argument values.
   * @param x Array of x values. Must be in increasing order.
   * @param y Array of y values. Must have same length as x.
   * @return Scalar or array of function values.
   */

  template <typename T0, typename T1, typename T2>
  typename boost::math::tools::promote_args <T0, T1, T2>::type
  linear_interpolation(const T0& xout,
                       const std::vector<T1>& x,
                       const std::vector<T2>& y) {
    typedef typename boost::math::tools::promote_args <T0, T1, T2>::type
      scalar;
    using std::vector;
    int nx = x.size();
    scalar yout;

    stan::math::check_finite("linear_interpolation", "xout", xout);
    stan::math::check_finite("linear_interpolation", "x", x);
    stan::math::check_finite("linear_interpolation", "y", y);
    stan::math::check_nonzero_size("linear_interpolation", "x", x);
    stan::math::check_nonzero_size("linear_interpolation", "y", y);
    stan::math::check_ordered("linear_interpolation", "x", x);
    stan::math::check_matching_sizes("linear_interpolation", "x", x, "y", y);

    if (xout <= x.front()) {
      yout = y.front();
    } else if (xout >= x.back()) {
      yout = y.back();
    } else {
      int i = SearchReal(x, nx, xout) - 1;
      yout = y[i] + (y[i+1] - y[i]) / (x[i+1] - x[i]) * (xout - x[i]);
    }

    return yout;
  }

  template <typename T0, typename T1, typename T2>
  std::vector <typename boost::math::tools::promote_args <T0, T1, T2>::type>
  linear_interpolation(const std::vector<T0>& xout,
                       const std::vector<T1>& x,
                       const std::vector<T2>& y) {
    typedef typename boost::math::tools::promote_args <T0, T1, T2>::type
      scalar;
    using std::vector;

    int nxout = xout.size();
    vector<scalar> yout(nxout);

    stan::math::check_nonzero_size("linear_interpolation", "xout", xout);
    stan::math::check_finite("linear_interpolation", "xout", xout);

    for (int i = 0; i < nxout; i++) {
      yout[i] = linear_interpolation(xout[i], x, y);
    }
    return yout;
  }

}
#endif
