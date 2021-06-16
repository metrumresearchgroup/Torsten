#ifndef STAN_MATH_TORSTEN_DSOLVE_INTEGRATE_ODE_GROUP_RK45_HPP
#define STAN_MATH_TORSTEN_DSOLVE_INTEGRATE_ODE_GROUP_RK45_HPP

#include <stan/math/torsten/mpi/pmx_population_integrator.hpp>
#include <stan/math/torsten/dsolve/pmx_integrate_ode_rk45.hpp>

namespace torsten {
  /**
   * Solve population ODE model by delegating the population
   * ODE integration task to multiple processors through
   * MPI, then gather the results, before generating @c var arrays.
   * Each entry has an additional level of nested vector to
   * identifiy the individual among a population of ODE parameters/data.
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
   * @return res nested vector that contains results for
   * (individual i, time j, equation k)
   **/
  template <typename F, typename Tt, typename T_initial, typename T_param>
  Eigen::Matrix<typename stan::return_type_t<Tt, T_initial, T_param>,
                Eigen::Dynamic, Eigen::Dynamic>
  pmx_integrate_ode_group_rk45(const F& f,
                               const std::vector<std::vector<T_initial> >& y0,
                               double t0,
                               const std::vector<int>& len,
                               const std::vector<Tt>& ts,
                               const std::vector<std::vector<T_param> >& theta,
                               const std::vector<std::vector<double> >& x_r,
                               const std::vector<std::vector<int> >& x_i,
                               double rtol,
                               double atol,
                               long int max_num_step,
                               std::ostream* msgs = nullptr) {
    static const char* caller("pmx_integrate_ode_rk45");
    dsolve::ode_group_check(y0, t0, len, ts, theta, x_r, x_i, caller);

    using scheme_t = boost::numeric::odeint::runge_kutta_dopri5<std::vector<double>, double, std::vector<double>, double>;
    dsolve::PMXOdeintIntegrator<scheme_t> integrator(rtol, atol, max_num_step);
    torsten::mpi::PMXPopulationIntegrator<F, dsolve::PMXOdeintIntegrator<scheme_t>,
                                          dsolve::PMXOdeSystem> solver(integrator);

    return solver(f, y0, t0, len, ts, theta, x_r, x_i, msgs);
  }

  /*
   * overload with default ode integrator controls
   */
  template <typename F, typename Tt, typename T_initial, typename T_param>
  Eigen::Matrix<typename stan::return_type_t<Tt, T_initial, T_param>,
                Eigen::Dynamic, Eigen::Dynamic>
  pmx_integrate_ode_group_rk45(const F& f,
                               const std::vector<std::vector<T_initial> >& y0,
                               double t0,
                               const std::vector<int>& len,
                               const std::vector<Tt>& ts,
                               const std::vector<std::vector<T_param> >& theta,
                               const std::vector<std::vector<double> >& x_r,
                               const std::vector<std::vector<int> >& x_i,
                               std::ostream* msgs = nullptr) {
    return pmx_integrate_ode_group_rk45(f, y0, t0, len, ts, theta, x_r, x_i, 
                                        1.e-6, 1.e-6, 1e6,
                                        msgs);
  }


}
#endif
