#ifndef TORSTEN_ODE_INTEGRATOR_HPP
#define TORSTEN_ODE_INTEGRATOR_HPP

#include <stan/math/torsten/dsolve/dsolve.hpp>
#include <stan/math/rev/mat.hpp>
#include <ostream>
#include <vector>

namespace torsten {
  enum PkOdeIntegratorId {
    StanRk45, StanAdams, StanBdf,
    PkAdams, PkBdf
  };

  template<PkOdeIntegratorId T>
  struct PkOdeIntegrator;

  template<>
  struct PkOdeIntegrator<StanRk45> {
    template <typename F, typename T_initial, typename T_param>
    std::vector<std::vector<typename stan::return_type<T_initial,
                                                       T_param>::type> >
    operator()(const F& f,
               const std::vector<T_initial>& y0,
               double t0,
               const std::vector<double>& ts,
               const std::vector<T_param>& theta,
               const std::vector<double>& x_r,
               const std::vector<int>& x_i,
               std::ostream* msgs,
               double rtol,
               double atol,
               long int max_num_step) {  // NOLINT(runtime/int)                   
      return stan::math::integrate_ode_rk45(f, y0, t0, ts, theta,
                                            x_r, x_i, msgs, rtol,
                                            atol, max_num_step);
    }
  };

  template<>
  struct PkOdeIntegrator<StanAdams> {
    template <typename F, typename T_initial, typename T_param>
    std::vector<std::vector<typename stan::return_type<T_initial,
                                                       T_param>::type> >
    operator()(const F& f,
               const std::vector<T_initial>& y0,
               double t0,
               const std::vector<double>& ts,
               const std::vector<T_param>& theta,
               const std::vector<double>& x_r,
               const std::vector<int>& x_i,
               std::ostream* msgs,
               double rtol,
               double atol,
               long int max_num_step) {  // NOLINT(runtime/int)                   
      return stan::math::integrate_ode_adams(f, y0, t0, ts, theta,
                                             x_r, x_i, msgs, rtol,
                                             atol, max_num_step);
    }
  };

  template<>
  struct PkOdeIntegrator<StanBdf> {
    template <typename F, typename T_initial, typename T_param>
    std::vector<std::vector<typename stan::return_type<T_initial,
                                                       T_param>::type> >
    operator()(const F& f,
               const std::vector<T_initial>& y0,
               double t0,
               const std::vector<double>& ts,
               const std::vector<T_param>& theta,
               const std::vector<double>& x_r,
               const std::vector<int>& x_i,
               std::ostream* msgs,
               double rtol,
               double atol,
               long int max_num_step) {  // NOLINT(runtime/int)                   
      return stan::math::integrate_ode_bdf(f, y0, t0, ts, theta,
                                           x_r, x_i, msgs, rtol,
                                           atol, max_num_step);
    }
  };

  template<>
  struct PkOdeIntegrator<PkAdams> {
    template <typename F, typename Tt, typename T_initial, typename T_param>
    std::vector<std::vector<typename stan::return_type<Tt, T_initial,
                                                       T_param>::type> >
    operator()(const F& f,
               const std::vector<T_initial>& y0,
               double t0,
               const std::vector<Tt>& ts,
               const std::vector<T_param>& theta,
               const std::vector<double>& x_r,
               const std::vector<int>& x_i,
               std::ostream* msgs,
               double rtol,
               double atol,
               long int max_num_step) {  // NOLINT(runtime/int)                   
      return torsten::dsolve::pk_integrate_ode_adams(f, y0, t0, ts, theta,
                                                     x_r, x_i, msgs, rtol,
                                                     atol, max_num_step);
    }
  };

  template<>
  struct PkOdeIntegrator<PkBdf> {
    template <typename F, typename Tt, typename T_initial, typename T_param>
    std::vector<std::vector<typename stan::return_type<Tt, T_initial,
                                                       T_param>::type> >
    operator()(const F& f,
               const std::vector<T_initial>& y0,
               double t0,
               const std::vector<Tt>& ts,
               const std::vector<T_param>& theta,
               const std::vector<double>& x_r,
               const std::vector<int>& x_i,
               std::ostream* msgs,
               double rtol,
               double atol,
               long int max_num_step) {  // NOLINT(runtime/int)                   
      return torsten::dsolve::pk_integrate_ode_bdf(f, y0, t0, ts, theta,
                                                   x_r, x_i, msgs, rtol,
                                                   atol, max_num_step);
    }
  };
}
#endif
