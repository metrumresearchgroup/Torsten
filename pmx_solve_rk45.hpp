#ifndef STAN_MATH_TORSTEN_SOLVE_RK45_HPP
#define STAN_MATH_TORSTEN_SOLVE_RK45_HPP

#include <Eigen/Dense>
#include <stan/math/torsten/to_array_2d.hpp>
#include <stan/math/torsten/ev_manager.hpp>
#include <stan/math/torsten/pmx_population_check.hpp>
#include <stan/math/torsten/ev_solver.hpp>
#include <stan/math/torsten/pmx_ode_model.hpp>
#include <stan/math/torsten/pmx_check.hpp>
#include <stan/math/torsten/pmx_solve_ode.hpp>
#include <stan/math/torsten/pmx_solve_group_ode.hpp>
#include <vector>

namespace torsten {

/**
 * Computes the predicted amounts in each compartment at each event
 * for a general compartment model, defined by a system of ordinary
 * differential equations. Uses the torsten::pmx_integrate_ode_rk45
 * function. 
 *
 * @tparam Ts types of parameters, see <code>pmx_solve_ode</code> for
 *         details. The last type of <code>Ts...</code> is
 *         <code>ostream*</code> as required by stan
 * @tparam F type of ODE system function.
 * @param[in] f functor for base ordinary differential equation that defines 
 *            compartment model.
 * @param[in] nCmt number of compartments in model
 * @return a matrix with predicted amount in each compartment 
 *         at each event. 
 *
 */
  template <typename F, typename... Ts>
  auto pmx_solve_rk45(const F& f, const int nCmt, Ts... args) {
    using scheme_t = torsten::dsolve::odeint_scheme_rk45;
    return PMXSolveODE<dsolve::PMXOdeIntegrator<dsolve::PMXVariadicOdeSystem, dsolve::PMXOdeintIntegrator<scheme_t>>>::solve(f, nCmt, args...);
}

  /*
   * For backward compatibility we keep old version of
   * return type using transpose. This is less efficient and
   * will be decomissioned in formal release.
   */
  template <typename T0, typename T1, typename T2, typename T3,
            typename T_par, typename T_biovar, typename T_tlag,
            typename F>
  stan::matrix_return_t<T0, T1, T2, T3, T_par, T_biovar, T_tlag>
  generalOdeModel_rk45(const F& f,
                       const int nCmt,
                       const std::vector<T0>& time,
                       const std::vector<T1>& amt,
                       const std::vector<T2>& rate,
                       const std::vector<T3>& ii,
                       const std::vector<int>& evid,
                       const std::vector<int>& cmt,
                       const std::vector<int>& addl,
                       const std::vector<int>& ss,
                       const std::vector<T_par>& pMatrix,
                       const std::vector<T_biovar>& biovar,
                       const std::vector<T_tlag>& tlag,
                       double rel_tol = 1e-6,
                       double abs_tol = 1e-6,
                       long int max_num_steps = 1e6,
                       std::ostream* msgs = 0) {
    auto x = pmx_solve_rk45(f, nCmt,
                            time, amt, rate, ii, evid, cmt, addl, ss,
                            pMatrix, biovar, tlag,
                            rel_tol, abs_tol, max_num_steps,
                            msgs);
    return x.transpose();
  }

  /*
   * For population models, more often we use ragged arrays
   * to describe the entire population, so in addition we need the arrays of
   * the length of each individual's data. The size of that
   * vector is the size of the population.
   */
  template <typename F, typename... Ts>
  auto pmx_solve_group_rk45(const F& f, const int nCmt,
                            const std::vector<int>& len, Ts... args) {
    using scheme_t = torsten::dsolve::odeint_scheme_rk45;
    return PMXSolveGroupODE<dsolve::PMXOdeIntegrator<dsolve::PMXVariadicOdeSystem, dsolve::PMXOdeintIntegrator<scheme_t>>>::solve(f, nCmt, len, args...);
  }
}

#endif
