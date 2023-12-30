#ifndef STAN_MATH_TORSTEN_ONECTP_HPP
#define STAN_MATH_TORSTEN_ONECTP_HPP

#include <Eigen/Dense>
#include <boost/math/tools/promotion.hpp>
#include <stan/math/torsten/to_array_2d.hpp>
#include <stan/math/torsten/pmx_solve_ode.hpp>
#include <stan/math/torsten/pmx_solve_cpt.hpp>
#include <stan/math/prim/err/check_positive_finite.hpp>
#include <stan/math/torsten/ev_solver.hpp>
#include <stan/math/torsten/ev_manager.hpp>
#include <stan/math/torsten/pmx_onecpt_model.hpp>
#include <string>
#include <vector>

namespace torsten {

/**
 * Computes the predicted amounts in each compartment at each event
 * for a one-compartment model with analytical solution.
 *
 * @tparam Ts types of parameters, see <code>pmx_solve_cpt</code> for
 *         details. 
 * @return a matrix with predicted amount in each compartment 
 *         at each event. 
 *
 */
  template <typename... Ts>
  auto pmx_solve_onecpt(Ts... args) {
    return PMXSolveCPT<PMXOneCptModel>::solve(args...);
  }

  // old version
  template <typename T0, typename T1, typename T2, typename T3, typename T4,
            typename T5, typename T6>
  stan::matrix_return_t<T0, T1, T2, T3, T4, T5, T6>
  PKModelOneCpt(const std::vector<T0>& time,
                const std::vector<T1>& amt,
                const std::vector<T2>& rate,
                const std::vector<T3>& ii,
                const std::vector<int>& evid,
                const std::vector<int>& cmt,
                const std::vector<int>& addl,
                const std::vector<int>& ss,
                const std::vector<std::vector<T4> >& pMatrix,
                const std::vector<std::vector<T5> >& biovar,
                const std::vector<std::vector<T6> >& tlag) {
    auto x = pmx_solve_onecpt(time, amt, rate, ii, evid, cmt, addl, ss,
                              pMatrix, biovar, tlag);
    return x.transpose();
  }

  template <typename T0, typename T1, typename T2, typename T3,
            typename T_par, typename T_biovar, typename T_tlag,
            typename = require_any_not_std_vector_t<T_par, T_biovar, T_tlag> >
  stan::matrix_return_t<T0, T1, T2, T3, T_par, T_biovar, T_tlag>
  PKModelOneCpt(const std::vector<T0>& time,
                   const std::vector<T1>& amt,
                   const std::vector<T2>& rate,
                   const std::vector<T3>& ii,
                   const std::vector<int>& evid,
                   const std::vector<int>& cmt,
                   const std::vector<int>& addl,
                   const std::vector<int>& ss,
                   const std::vector<T_par>& pMatrix,
                   const std::vector<T_biovar>& biovar,
                   const std::vector<T_tlag>& tlag) {
    auto x = pmx_solve_onecpt(time, amt, rate, ii, evid, cmt, addl, ss,
                            pMatrix, biovar, tlag);
    return x.transpose();
  }

}

#endif
