#ifndef STAN_MATH_TORSTEN_FINITE_DIFF_GRADIENT_HPP
#define STAN_MATH_TORSTEN_FINITE_DIFF_GRADIENT_HPP

#include <stan/math/prim/fun/Eigen.hpp>

namespace torsten {

/**
 * Calculate the value and the gradient of the specified function
 * at the specified argument using finite difference.
 *
 * <p>The functor must implement
 *
 * <code>
 * double
 * operator()(const
 * Eigen::Matrix<double, Eigen::Dynamic, 1>&)
 * </code>
 *
 * Error should be on order of epsilon ^ 6.
 * The reference for this algorithm is:
 *
 * De Levie: An improved numerical approximation
 * for the first derivative, page 3
 *
 * This function involves 6 calls to f.
 *
 * @tparam F Type of function
 * @param[in] f Function that returns a scalar
 * @param[in] x Argument to function
 * @param[out] fx Function applied to argument
 * @param[out] grad_fx Gradient of function at argument
 * @param[in] epsilon perturbation size
 */
template <typename F>
auto finite_diff_gradient(const F& f, const std::vector<double>& x,
                          double epsilon = 1e-03) {
  std::vector<double> x_temp(x);

  int d = x.size();
  using scalar_t = typename std::decay<decltype(f(x))>::type;
  std::vector<scalar_t> grad_fx(d);

  for (int i = 0; i < d; ++i) {
    scalar_t delta_f = 0.0;

    x_temp[i] = x[i] + 3.0 * epsilon;
    delta_f = f(x_temp);

    x_temp[i] = x[i] + 2.0 * epsilon;
    delta_f -= 9.0 * f(x_temp);

    x_temp[i] = x[i] + epsilon;
    delta_f += 45.0 * f(x_temp);

    x_temp[i] = x[i] + -3.0 * epsilon;
    delta_f -= f(x_temp);

    x_temp[i] = x[i] + -2.0 * epsilon;
    delta_f += 9.0 * f(x_temp);

    x_temp[i] = x[i] + -epsilon;
    delta_f -= 45.0 * f(x_temp);

    delta_f /= 60 * epsilon;

    x_temp[i] = x[i];
    grad_fx[i] = delta_f;
  }
  return grad_fx;
}

  template<typename T>
  inline const T& ith_of(const std::vector<T, T>&& d, int i) {
    return d.at(i);
  }

  template<typename T, int R, int C>
  inline const T& ith_of(const Eigen::Matrix<T, R, C>&& d, int i) {
    return d(i);
  }

  /*
   * @param[in] f Function that returns a vector
   */
template <typename F>
auto finite_diff_gradient(const F& f, const std::vector<double>& x, int k,
                          double epsilon = 1e-03) {
  auto f1 = [&f, &k](const std::vector<double>& x) {
    return ith_of(f(x), k);
  };

  return finite_diff_gradient(f1, x, epsilon);
}

}  // namespace torsten
#endif
