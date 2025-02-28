#ifndef STAN_MATH_FWD_CORE_OPERATOR_LOGICAL_OR_HPP
#define STAN_MATH_FWD_CORE_OPERATOR_LOGICAL_OR_HPP

#include <stan/math/fwd/core/fvar.hpp>

namespace stan {
namespace math {

/**
 * Return the logical disjunction of the values of the two
 * arguments as defined by <code>||</code>.
 *
 * @tparam value and tangent type for variables
 * @param[in] x first argument
 * @param[in] y second argument
 * @return disjunction of the arguments' values
 */
template <typename T>
inline bool operator||(const fvar<T>& x, const fvar<T>& y) {
  return x.val_ || y.val_;
}

/**
 * Return the logical disjunction of the values of the two
 * arguments as defined by <code>||</code>.
 *
 * @tparam value and tangent type for variables
 * @param[in] x first argument
 * @param[in] y second argument
 * @return disjunction of first argument's value and second
 * argument
 */
template <typename T>
inline bool operator||(const fvar<T>& x, double y) {
  return x.val_ || y;
}

/**
 * Return the logical disjunction of the values of the two
 * arguments as defined by <code>||</code>.
 *
 * @tparam value and tangent type for variables
 * @param[in] x first argument
 * @param[in] y second argument
 * @return disjunction of first argument and the second
 * argument's value
 */
template <typename T>
inline bool operator||(double x, const fvar<T>& y) {
  return x || y.val_;
}

}  // namespace math
}  // namespace stan
#endif
