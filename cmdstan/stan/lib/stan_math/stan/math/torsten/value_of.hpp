#ifndef STAN_MATH_TORSTEN_FUN_VALUE_OF_HPP
#define STAN_MATH_TORSTEN_FUN_VALUE_OF_HPP

#include <stan/math/prim/meta.hpp>
#include <stan/math/prim/fun/Eigen.hpp>
#include <stan/math/prim/fun/value_of.hpp>
#include <stan/math/rev/fun/value_of.hpp>
#include <cstddef>
#include <vector>

namespace torsten {

/**
 * Inputs that are arithmetic types or containers of airthmetric types
 * are returned from value_of unchanged
 *
 * @tparam T Input type
 * @param[in] x Input argument
 * @return Forwarded input argument
 **/
template <typename T, stan::require_st_arithmetic<T>* = nullptr>
inline T value_of(T&& x) {
  return stan::math::value_of(x);
}

/**
 * Return the value of the specified variable.
 *
 * <p>This function is used internally by autodiff functions along
 * with <code>value_of(T x)</code> to extract the
 * <code>double</code> value of either a scalar or an autodiff
 * variable.  This function will be called when the argument is a
 * <code>var</code> even if the function is not
 * referred to by namespace because of argument-dependent lookup.
 *
 * @param v Variable.
 * @return Value of variable.
 */
template <typename T>
inline auto& value_of(const stan::math::var_value<T>& v) {
  return stan::math::value_of(v);
}

/**
 * For std::vectors of non-arithmetic types, return a std::vector composed
 * of value_of applied to each element.
 *
 * @tparam T Input element type
 * @param[in] x Input std::vector
 * @return std::vector of values
 **/
template <typename T, stan::require_not_st_arithmetic<T>* = nullptr>
inline auto value_of(const std::vector<T>& x) {
  return stan::math::value_of(x);
}

/**
 * For Eigen matrices and expressions of non-arithmetic types, return an
 *expression that represents the Eigen::Matrix resulting from applying value_of
 *elementwise
 *
 * @tparam EigMat type of the matrix
 *
 * @param[in] M Matrix to be converted
 * @return Matrix of values
 **/
  template <typename T>
  inline Eigen::VectorXd value_of(const Eigen::Matrix<T, -1, 1>& M) {
    Eigen::VectorXd res(M.size());
    for (auto i = 0; i < M.size(); ++i) {
      res[i] = stan::math::value_of(M[i]);
    }
    return res;
  }

  template <typename T, int R, int C>
  inline Eigen::Matrix<double, R, C> value_of(const Eigen::Matrix<T, R, C>& M) {
    Eigen::Matrix<double, R, C> res(R, C);
    for (auto i = 0; i < M.size(); ++i) {
      res(i) = stan::math::value_of(M(i));
    }
    return res;
  }
}

#endif
