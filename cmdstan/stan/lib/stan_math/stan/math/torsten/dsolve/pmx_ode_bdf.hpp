#ifndef STAN_MATH_TORSTEN_ODE_BDF_HPP
#define STAN_MATH_TORSTEN_ODE_BDF_HPP

#include <stan/math/torsten/dsolve/pmx_ode_integrator.hpp>
#include <stan/math/torsten/dsolve/pmx_cvodes_integrator.hpp>
#include <ostream>
#include <vector>

namespace torsten {
  /**
   * solve an ODE given its RHS with CVODES' BDF solver.
   *
   * @tparam F Functor type for RHS of ODE
   * @tparam Tt type of time
   * @tparam T_initial type of initial condition @c y0
   * @tparam T_param type of parameter @c theta
   * @param f RHS functor of ODE system
   * @param y0 initial condition
   * @param t0 initial time
   * @param ts time steps
   * @param theta parameters for ODE
   * @param x_r data used in ODE
   * @param x_i integer data used in ODE
   * @param msgs output stream
   * @param rtol relative tolerance
   * @param atol absolute tolerance
   * @param max_num_step maximum number of integration steps allowed.
   * @return a vector of vectors for results in each time step.
   */
  template <typename F, typename Tt, typename T_initial, typename... T_param>
  std::vector<Eigen::Matrix<typename stan::return_type_t<Tt, T_initial, T_param...>,-1,1> >
  pmx_ode_bdf_ctrl(const F& f,
                   const Eigen::Matrix<T_initial, -1, 1>& y0,
                   double t0,
                   const std::vector<Tt>& ts,
                   std::ostream* msgs,
                   double rtol, double atol, int max_num_step,
                   const T_param&... args) {
    using dsolve::PMXOdeIntegrator;
    using dsolve::PMXCvodesIntegrator;
    using dsolve::PMXVariadicOdeSystem;
    PMXOdeIntegrator<PMXVariadicOdeSystem, PMXCvodesIntegrator<CV_BDF, CV_STAGGERED>> solver(rtol, atol, max_num_step, msgs);
    return solver(f, y0, t0, ts, args...);
  }

  /**
   * overload with default ode controls
   */
  template <typename F, typename Tt, typename T_initial, typename... T_param>
  std::vector<Eigen::Matrix<typename stan::return_type_t<Tt, T_initial, T_param...>,-1,1> >
  pmx_ode_bdf(const F& f,
              const Eigen::Matrix<T_initial, -1, 1>& y0,
              double t0,
              const std::vector<Tt>& ts,
              std::ostream* msgs,
              const T_param&... args) {
    return pmx_ode_bdf_ctrl(f, y0, t0, ts, msgs, 1.e-10, 1.e-10, 1e8, args...);
  }
}
#endif
