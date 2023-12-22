#ifndef STAN_MATH_REV_CORE_OPERATOR_LESS_THAN_HPP
#define STAN_MATH_REV_CORE_OPERATOR_LESS_THAN_HPP

#include <stan/math/rev/core/var.hpp>
#include <stan/math/prim/meta.hpp>

namespace stan {
namespace math {

/**
 * Less than operator comparing variables' values (C++).
 *
 \f[
   \mbox{operator\textless}(x, y) =
   \begin{cases}
     0 & \mbox{if } x \geq y \\
     1 & \mbox{if } x < y \\[6pt]
     0 & \mbox{if } x = \textrm{NaN or } y = \textrm{NaN}
   \end{cases}
 \f]
 * @param a First variable.
 * @param b Second variable.
 * @return True if first variable's value is less than second's.
 */
inline bool operator<(const var& a, const var& b) { return a.val() < b.val(); }

/**
 * Less than operator comparing variable's value and a double
 * (C++).
 *
 * @tparam Arith An arithmetic type
 * @param a First variable.
 * @param b Second value.
 * @return True if first variable's value is less than second value.
 */
template <typename Arith, require_arithmetic_t<Arith>* = nullptr>
inline bool operator<(const var& a, Arith b) {
  return a.val() < b;
}

/**
 * Less than operator comparing a double and variable's value
 * (C++).
 *
 * @tparam Arith An arithmetic type
 * @param a First value.
 * @param b Second variable.
 * @return True if first value is less than second variable's value.
 */
template <typename Arith, require_arithmetic_t<Arith>* = nullptr>
inline bool operator<(Arith a, const var& b) {
  return a < b.val();
}

}  // namespace math
}  // namespace stan
#endif
