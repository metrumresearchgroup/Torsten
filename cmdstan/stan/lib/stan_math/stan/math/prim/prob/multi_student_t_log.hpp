#ifndef STAN_MATH_PRIM_PROB_MULTI_STUDENT_T_LOG_HPP
#define STAN_MATH_PRIM_PROB_MULTI_STUDENT_T_LOG_HPP

#include <stan/math/prim/meta.hpp>
#include <stan/math/prim/prob/multi_student_t_lpdf.hpp>

namespace stan {
namespace math {

/** \ingroup multivar_dists
 * Return the log of the multivariate Student t distribution
 * at the specified arguments.
 *
 * @deprecated use <code>multi_student_t_lpdf</code>
 *
 * @tparam propto Carry out calculations up to a proportion
 */
template <bool propto, typename T_y, typename T_dof, typename T_loc,
          typename T_scale>
return_type_t<T_y, T_dof, T_loc, T_scale> multi_student_t_log(
    const T_y& y, const T_dof& nu, const T_loc& mu, const T_scale& Sigma) {
  return multi_student_t_lpdf<propto, T_y, T_dof, T_loc, T_scale>(y, nu, mu,
                                                                  Sigma);
}

/** \ingroup multivar_dists
 * @deprecated use <code>multi_student_t_lpdf</code>
 */
template <typename T_y, typename T_dof, typename T_loc, typename T_scale>
inline return_type_t<T_y, T_dof, T_loc, T_scale> multi_student_t_log(
    const T_y& y, const T_dof& nu, const T_loc& mu, const T_scale& Sigma) {
  return multi_student_t_lpdf<T_y, T_dof, T_loc, T_scale>(y, nu, mu, Sigma);
}

}  // namespace math
}  // namespace stan
#endif
