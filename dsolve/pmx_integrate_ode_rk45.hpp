#ifndef STAN_MATH_TORSTEN_DSOLVE_INTEGRATE_ODE_RK45_HPP
#define STAN_MATH_TORSTEN_DSOLVE_INTEGRATE_ODE_RK45_HPP

#include <stan/math/torsten/dsolve/pmx_odeint_system.hpp>
#include <stan/math/torsten/dsolve/pmx_odeint_integrator.hpp>
#include <stan/math/torsten/dsolve/ode_check.hpp>
#include <ostream>
#include <vector>

namespace torsten {
  /*
   * solve an ODE given its RHS with Boost Odeint's Rk45 solver.
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
  inline std::vector<std::vector<typename stan::return_type<Tt, T_initial, T_param>::type> >
  pmx_integrate_ode_rk45(const F& f,
                         const std::vector<T_initial>& y0,
                         double t0,
                         const std::vector<Tt>& ts,
                         const std::vector<T_param>& theta,
                         const std::vector<double>& x_r,
                         const std::vector<int>& x_i,
                         double rtol,
                         double atol,
                         long int max_num_steps,
                         std::ostream* msgs = nullptr) {
    static const char* caller = "pmx_integrate_ode_rk45";
    dsolve::ode_check(y0, t0, ts, theta, x_r, x_i, caller);

    using Ode = dsolve::PMXOdeintSystem<F, Tt, T_initial, T_param>;
    const int m = theta.size();
    const int n = y0.size();

    dsolve::PMXOdeService<Ode, dsolve::Odeint> serv(n, m);

    Ode ode{serv, f, t0, ts, y0, theta, x_r, x_i, msgs};
    using scheme_t = boost::numeric::odeint::runge_kutta_dopri5<std::vector<double>, double, std::vector<double>, double>;
    dsolve::PMXOdeintIntegrator<scheme_t> solver(rtol, atol, max_num_steps);

    return solver.integrate(ode);
}

  /*
   * overload with default ode controls
   */
  template <typename F, typename Tt, typename T_initial, typename T_param>
  inline std::vector<std::vector<typename stan::return_type<Tt, T_initial, T_param>::type> >
  pmx_integrate_ode_rk45(const F& f,
                         const std::vector<T_initial>& y0,
                         double t0,
                         const std::vector<Tt>& ts,
                         const std::vector<T_param>& theta,
                         const std::vector<double>& x_r,
                         const std::vector<int>& x_i,
                         std::ostream* msgs = nullptr) {
    return pmx_integrate_ode_rk45(f, y0, t0, ts, theta, x_r, x_i,
                                  1.e-6, 1.e-6, 1e6,
                                  msgs);
  }
}
#endif
