#ifndef STAN_MATH_TORSTEN_DSOLVE_INTEGRATE_ODE_BDF_HPP
#define STAN_MATH_TORSTEN_DSOLVE_INTEGRATE_ODE_BDF_HPP

#include <stan/math/torsten/dsolve/pmx_cvodes_integrator.hpp>
#include <stan/math/torsten/dsolve/ode_check.hpp>
#include <stan/math/torsten/mpi.hpp>
#include <ostream>
#include <vector>

namespace torsten {

  /*
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
  template <typename F, typename Tt, typename T_initial, typename T_param>
  std::vector<std::vector<typename stan::return_type<Tt, T_initial, T_param>::type> >
  pmx_integrate_ode_bdf(const F& f,
                         const std::vector<T_initial>& y0,
                         double t0,
                         const std::vector<Tt>& ts,
                         const std::vector<T_param>& theta,
                         const std::vector<double>& x_r,
                         const std::vector<int>& x_i,
                         std::ostream* msgs = nullptr,
                         double rtol = 1e-10,
                         double atol = 1e-10,
                         long int max_num_step = 1e8) {
    static const char* caller = "pmx_integrate_ode_bdf";
    dsolve::ode_check(y0, t0, ts, theta, x_r, x_i, caller);

    using Ode = dsolve::PMXCvodesFwdSystem<F, Tt, T_initial, T_param, CV_BDF, AD>;
    const int m = theta.size();
    const int n = y0.size();

    static dsolve::PMXOdeService<typename Ode::Ode> serv(n, m);

    Ode ode{serv, f, t0, ts, y0, theta, x_r, x_i, msgs};
    dsolve::PMXCvodesIntegrator solver(rtol, atol, max_num_step);
    return solver.integrate(ode);
}
}
#endif
