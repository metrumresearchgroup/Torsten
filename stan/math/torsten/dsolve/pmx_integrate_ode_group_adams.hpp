#ifndef STAN_MATH_TORSTEN_DSOLVE_INTEGRATE_ODE_GROUP_ADAMS_HPP
#define STAN_MATH_TORSTEN_DSOLVE_INTEGRATE_ODE_GROUP_ADAMS_HPP

#include <stan/math/torsten/mpi/pmx_population_integrator.hpp>
#include <stan/math/torsten/dsolve/pmx_integrate_ode_adams.hpp>
#include <stan/math/torsten/dsolve/ode_check.hpp>

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
  Eigen::Matrix<typename torsten::return_t<Tt, T_initial, T_param>::type,
                Eigen::Dynamic, Eigen::Dynamic>
  pmx_integrate_ode_group_adams(const F& f,
                                const std::vector<std::vector<T_initial> >& y0,
                                double t0,
                                const std::vector<int>& len,
                                const std::vector<Tt>& ts,
                                const std::vector<std::vector<T_param> >& theta,
                                const std::vector<std::vector<double> >& x_r,
                                const std::vector<std::vector<int> >& x_i,
                                std::ostream* msgs = nullptr,
                                double rtol = 1e-10,
                                double atol = 1e-10,
                                long int max_num_step = 1e6) {  // NOLINT(runtime/int)
    static const char* caller("pmx_integrate_ode_adams");
    dsolve::ode_group_check(y0, t0, len, ts, theta, x_r, x_i, caller);

    dsolve::PMXCvodesIntegrator integrator(rtol, atol, max_num_step);

    torsten::mpi::PMXPopulationIntegrator<F, dsolve::PMXCvodesIntegrator,
                                          dsolve::PMXCvodesFwdSystem_adams_ad> solver(integrator);

    return solver(f, y0, t0, len, ts, theta, x_r, x_i, msgs);
  }

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
   * @param group_theta parameters shared among the ODE group
   * @param theta parameters for each ODE system
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
  Eigen::Matrix<typename torsten::return_t<Tt, T_initial, T_param>::type,
                Eigen::Dynamic, Eigen::Dynamic>
  pmx_integrate_ode_group_adams(const F& f,
                                const std::vector<std::vector<T_initial> >& y0,
                                double t0,
                                const std::vector<int>& len,
                                const std::vector<Tt>& ts,
                                const std::vector<T_param>& group_theta,
                                const std::vector<std::vector<T_param> >& theta,
                                const std::vector<std::vector<double> >& x_r,
                                const std::vector<std::vector<int> >& x_i,
                                std::ostream* msgs = nullptr,
                                double rtol = 1e-10,
                                double atol = 1e-10,
                                long int max_num_step = 1e6) {  // NOLINT(runtime/int)
    static const char* caller("pmx_integrate_ode_adams");
    dsolve::ode_group_check(y0, t0, len, ts, theta, x_r, x_i, caller);

    dsolve::PMXCvodesIntegrator integrator(rtol, atol, max_num_step);

    torsten::mpi::PMXPopulationIntegrator<F, dsolve::PMXCvodesIntegrator,
                                          dsolve::PMXCvodesFwdSystem_adams_ad> solver(integrator);

    return solver(f, y0, t0, len, ts, group_theta, theta, x_r, x_i, msgs);
  }

}
#endif
