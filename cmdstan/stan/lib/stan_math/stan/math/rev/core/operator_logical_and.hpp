#ifndef STAN_MATH_REV_CORE_OPERATOR_LOGICAL_AND_HPP
#define STAN_MATH_REV_CORE_OPERATOR_LOGICAL_AND_HPP

#include <stan/math/rev/core/var.hpp>
#include <stan/math/prim/meta.hpp>

namespace stan {
namespace math {

/**
 * Return the logical conjunction of the values of the two
 * arguments as defined by <code>&amp;&amp;</code>.
 *
 * @param[in] x first argument
 * @param[in] y second argument
 * @return conjunction of the arguments' values
 */
inline bool operator&&(const var& x, const var& y) {
  return x.val() && y.val();
}

/**
 * Return the logical conjunction of the values of the two
 * arguments as defined by <code>&amp;&amp;</code>.
 *
 * @tparam Var value type of a var
 * @tparam Arith An arithmetic type
 * @param[in] x first argument
 * @param[in] y second argument
 * @return conjunction of first argument's value and second
 * argument
 */
template <typename Arith, require_arithmetic_t<Arith>* = nullptr>
inline bool operator&&(const var& x, Arith y) {
  return x.val() && y;
}

/**
 * Return the logical conjunction of the values of the two
 * arguments as defined by <code>&amp;&amp;</code>.
 *
 * @tparam Var value type of a var
 * @tparam Arith An arithmetic type
 * @param[in] x first argument
 * @param[in] y second argument
 * @return conjunction of first argument and second argument's
 * value
 */
template <typename Arith, require_arithmetic_t<Arith>* = nullptr>
inline bool operator&&(Arith x, const var& y) {
  return x && y.val();
}

}  // namespace math
}  // namespace stan
#endif
