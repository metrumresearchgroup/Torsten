#ifndef STAN_MATH_TORSTEN_DSOLVE_INTEGRATE_ODE_BDF_HPP
#define STAN_MATH_TORSTEN_DSOLVE_INTEGRATE_ODE_BDF_HPP

#include <stan/math/torsten/dsolve/pk_cvodes_integrator.hpp>
#include <ostream>
#include <vector>

namespace torsten {
namespace dsolve {

  template <typename F, typename Tt, typename T_initial, typename T_param>
  std::vector<std::vector<typename stan::return_type<Tt,
                                                     T_initial,
                                                     T_param>::type> >
  pk_integrate_ode_bdf(const F& f,
                         const std::vector<T_initial>& y0,
                         double t0,
                         const std::vector<Tt>& ts,
                         const std::vector<T_param>& theta,
                         const std::vector<double>& x_r,
                         const std::vector<int>& x_i,
                         std::ostream* msgs = nullptr,
                         double rtol = 1e-10,
                         double atol = 1e-10,
                         long int max_num_step = 1e6) {  // NOLINT(runtime/int)
    using torsten::dsolve::pk_cvodes_fwd_system;
    using torsten::dsolve::pk_cvodes_integrator;
    using Ode = pk_cvodes_fwd_system<F, Tt, T_initial, T_param, CV_BDF>;
    const int k = theta.size();
    const int n = y0.size();
    const int m = (Ode::is_var_y0 ? n : 0) + (Ode::is_var_par ? k : 0);

    static auto frhs = Ode::rhs();
    static cvodes_service<typename Ode::Ode> serv(n, m, frhs);

    Ode ode{serv, f, t0, ts, y0, theta, x_r, x_i, msgs};
    pk_cvodes_integrator solver(rtol, atol, max_num_step);
    return solver.integrate(ode);
}

}
}
#endif
