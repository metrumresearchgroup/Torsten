#ifndef STAN_MATH_TORSTEN_DSOLVE_INTEGRATE_DAE_HPP
#define STAN_MATH_TORSTEN_DSOLVE_INTEGRATE_DAE_HPP

#include <stan/math/torsten/dsolve/pmx_idas_fwd_system.hpp>
#include <stan/math/torsten/dsolve/pmx_idas_integrator.hpp>
#include <ostream>
#include <vector>

namespace torsten {
namespace dsolve {
/**
 * Return the solutions for a semi-explicit DAE system with residual
 * specified by functor F,
 * given the specified consistent initial state yy0 and yp0.
 *
 * @tparam DAE type of DAE system
 * @tparam Tpar scalar type of parameter theta
 * @param[in] f functor for the base ordinary differential equation
 * @param[in] yy0 initial state
 * @param[in] yp0 initial derivative state
 * @param[in] t0 initial time
 * @param[in] ts times of the desired solutions, in strictly
 * increasing order, all greater than the initial time
 * @param[in] theta parameters
 * @param[in] x_r real data
 * @param[in] x_i int data
 * @param[in] rtol relative tolerance passed to IDAS, requred <10^-3
 * @param[in] atol absolute tolerance passed to IDAS, problem-dependent
 * @param[in] max_num_steps maximal number of admissable steps
 * between time-points
 * @param[in] msgs message
 * @return a vector of states, each state being a vector of the
 * same size as the state variable, corresponding to a time in ts.
 */
template <typename F, typename Tpar>
std::vector<std::vector<Tpar> > pmx_integrate_dae(
    const F& f, const std::vector<double>& yy0, const std::vector<double>& yp0,
    double t0, const std::vector<double>& ts, const std::vector<Tpar>& theta,
    const std::vector<double>& x_r, const std::vector<int>& x_i,
    const double rtol, const double atol,
    const int64_t max_num_steps = PMXIdasIntegrator::IDAS_MAX_STEPS,
    std::ostream* msgs = nullptr) {
  /* it doesn't matter here what values @c eq_id has, as we
     don't allow yy0 or yp0 to be parameters */
  const std::vector<int> dummy_eq_id(yy0.size(), 0);

  torsten::dsolve::PMXIdasIntegrator solver(rtol, atol, max_num_steps);
  torsten::dsolve::PMXIdasFwdSystem<F, double, double, Tpar> dae{
      f, dummy_eq_id, yy0, yp0, theta, x_r, x_i, msgs};

  dae.check_ic_consistency(t0, atol);

  return solver.integrate(dae, t0, ts);
}
}  // namespace dsolve
}  // namespace torsten

#endif
