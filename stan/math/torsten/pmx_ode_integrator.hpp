#ifndef STAN_MATH_TORSTEN_PMX_ODE_INTEGRATOR_HPP
#define STAN_MATH_TORSTEN_PMX_ODE_INTEGRATOR_HPP

#include <stan/math/prim/scal/meta/return_type.hpp>
#include <stan/math/prim/arr/fun/value_of.hpp>
#include <stan/math/prim/arr/functor/integrate_ode_rk45.hpp>
#include <stan/math/rev/mat/functor/integrate_ode_adams.hpp>
#include <stan/math/rev/mat/functor/integrate_ode_bdf.hpp>
#include <stan/math/torsten/dsolve.hpp>
#include <ostream>
#include <vector>

namespace torsten {
  enum PMXOdeIntegratorId {
    Analytical,
    StanRk45, StanAdams, StanBdf,
    PkAdams, PkBdf, PkRk45
  };
}

namespace torsten {
  template<PMXOdeIntegratorId It>
  struct PMXOdeIntegrator;

  template<>
  struct PMXOdeIntegrator<StanAdams> {
    const double rtol;
    const double atol;
    const long int max_num_step;
    std::ostream* msgs;

    PMXOdeIntegrator() : rtol(1e-10), atol(1e-10), max_num_step(1e8), msgs(0) {}

    PMXOdeIntegrator(const double rtol0, const double atol0, const long int max_num_step0,
                     std::ostream* msgs0) :
      rtol(rtol0), atol(atol0), max_num_step(max_num_step0), msgs(msgs0)
    {}
    
    /*
     * Stan's ODE solver function doesn't support @c T_t to
     * be @c var
     */
    template <typename F, typename Tt, typename T_initial, typename T_param>
    std::vector<std::vector<typename stan::return_type<T_initial, T_param>::type> > // NOLINT
    operator()(const F& f,
               const std::vector<T_initial>& y0,
               double t0,
               const std::vector<Tt>& ts,
               const std::vector<T_param>& theta,
               const std::vector<double>& x_r,
               const std::vector<int>& x_i) const {
      return stan::math::integrate_ode_adams(f, y0, t0, ts, theta, x_r, x_i, msgs, rtol, atol, max_num_step);
    }
  };

  template<>
  struct PMXOdeIntegrator<StanBdf> {
    const double rtol;
    const double atol;
    const long int max_num_step;
    std::ostream* msgs;

    PMXOdeIntegrator() : rtol(1e-10), atol(1e-10), max_num_step(1e8), msgs(0) {}

    PMXOdeIntegrator(const double rtol0, const double atol0, const long int max_num_step0,
                     std::ostream* msgs0) :
      rtol(rtol0), atol(atol0), max_num_step(max_num_step0), msgs(msgs0)
    {}
    
    /*
     * Stan's ODE solver function doesn't support @c T_t to
     * be @c var
     */
    template <typename F, typename Tt, typename T_initial, typename T_param>
    std::vector<std::vector<typename stan::return_type<T_initial, T_param>::type> > // NOLINT
    operator()(const F& f,
               const std::vector<T_initial>& y0,
               double t0,
               const std::vector<Tt>& ts,
               const std::vector<T_param>& theta,
               const std::vector<double>& x_r,
               const std::vector<int>& x_i) const {
      return stan::math::integrate_ode_bdf(f, y0, t0, ts, theta, x_r, x_i, msgs, rtol, atol, max_num_step);
    }
  };

  template<>
  struct PMXOdeIntegrator<StanRk45> {
    const double rtol;
    const double atol;
    const long int max_num_step;
    std::ostream* msgs;

    PMXOdeIntegrator() : rtol(1e-10), atol(1e-10), max_num_step(1e8), msgs(0) {}

    PMXOdeIntegrator(const double rtol0, const double atol0, const long int max_num_step0,
                     std::ostream* msgs0) :
      rtol(rtol0), atol(atol0), max_num_step(max_num_step0), msgs(msgs0)
    {}
    
    /*
     * Stan's ODE solver function doesn't support @c T_t to
     * be @c var
     */
    template <typename F, typename Tt, typename T_initial, typename T_param>
    std::vector<std::vector<typename stan::return_type<T_initial, T_param>::type> > // NOLINT
    operator()(const F& f,
               const std::vector<T_initial>& y0,
               double t0,
               const std::vector<Tt>& ts,
               const std::vector<T_param>& theta,
               const std::vector<double>& x_r,
               const std::vector<int>& x_i) const {
      return stan::math::integrate_ode_rk45(f, y0, t0, ts, theta, x_r, x_i, msgs, rtol, atol, max_num_step);
    }
  };

  /*
   * specialization for @c PkBdf
   */
  template<>
  struct PMXOdeIntegrator<PkBdf> {
    const double rtol;
    const double atol;
    const long int max_num_step;
    std::ostream* msgs;

    PMXOdeIntegrator() : rtol(1e-10), atol(1e-10), max_num_step(1e8), msgs(0) {}

    PMXOdeIntegrator(const double rtol0, const double atol0, const long int max_num_step0,
                     std::ostream* msgs0) :
      rtol(rtol0), atol(atol0), max_num_step(max_num_step0), msgs(msgs0)
    {}

    /*
     * Torsten's ODE solvers support @c T_t to be @c var
     */
    template <typename F, typename Tt, typename T_initial, typename T_param>
    std::vector<std::vector<typename stan::return_type<Tt, T_initial, T_param>::type> > // NOLINT
    operator()(const F& f,
               const std::vector<T_initial>& y0,
               double t0,
               const std::vector<Tt>& ts,
               const std::vector<T_param>& theta,
               const std::vector<double>& x_r,
               const std::vector<int>& x_i) const {
      return torsten::pmx_integrate_ode_bdf(f, y0, t0, ts, theta, x_r, x_i, msgs, rtol, atol, max_num_step);
    }

    /*
     * For MPI solution we need to return data that consists
     * of solution value and gradients
     */
    template <typename F, typename Tt, typename T_initial, typename T_param>
    Eigen::MatrixXd
    solve_d(const F& f,
            const std::vector<T_initial>& y0,
            double t0,
            const std::vector<Tt>& ts,
            const std::vector<T_param>& theta,
            const std::vector<double>& x_r,
            const std::vector<int>& x_i) const {
      using Ode = torsten::dsolve::PMXCvodesFwdSystem<F, Tt, T_initial, T_param, CV_BDF, AD>;
      const int m = theta.size();
      const int n = y0.size();

      dsolve::PMXOdeService<typename Ode::Ode> serv(n, m);

      Ode ode{serv, f, t0, ts, y0, theta, x_r, x_i, msgs};

      torsten::dsolve::PMXCvodesIntegrator solver(rtol, atol, max_num_step);
      Eigen::MatrixXd res = solver.integrate<Ode, false>(ode);

      return res;
    }
  };

  /*
   * specialization for @c PkAdams
   */
  template<>
  struct PMXOdeIntegrator<PkAdams> {
    const double rtol;
    const double atol;
    const long int max_num_step;
    std::ostream* msgs;

    PMXOdeIntegrator() : rtol(1e-10), atol(1e-10), max_num_step(1e8), msgs(0) {}

    PMXOdeIntegrator(const double rtol0, const double atol0, const long int max_num_step0,
                     std::ostream* msgs0) :
      rtol(rtol0), atol(atol0), max_num_step(max_num_step0), msgs(msgs0)
    {}

    /*
     * Torsten's ODE solvers support @c T_t to be @c var
     */
    template <typename F, typename Tt, typename T_initial, typename T_param>
    std::vector<std::vector<typename stan::return_type<Tt, T_initial, T_param>::type> > // NOLINT
    operator()(const F& f,
               const std::vector<T_initial>& y0,
               double t0,
               const std::vector<Tt>& ts,
               const std::vector<T_param>& theta,
               const std::vector<double>& x_r,
               const std::vector<int>& x_i) const {
      return torsten::pmx_integrate_ode_adams(f, y0, t0, ts, theta, x_r, x_i, msgs, rtol, atol, max_num_step);
    }

    /*
     * For MPI solution we need to return data that consists
     * of solution value and gradients
     */
    template <typename F, typename Tt, typename T_initial, typename T_param>
    Eigen::MatrixXd
    solve_d(const F& f,
            const std::vector<T_initial>& y0,
            double t0,
            const std::vector<Tt>& ts,
            const std::vector<T_param>& theta,
            const std::vector<double>& x_r,
            const std::vector<int>& x_i) const {
      using Ode = torsten::dsolve::PMXCvodesFwdSystem<F, Tt, T_initial, T_param, CV_ADAMS, AD>;
      const int m = theta.size();
      const int n = y0.size();

      dsolve::PMXOdeService<typename Ode::Ode> serv(n, m);

      Ode ode{serv, f, t0, ts, y0, theta, x_r, x_i, msgs};

      torsten::dsolve::PMXCvodesIntegrator solver(rtol, atol, max_num_step);
      Eigen::MatrixXd res = solver.integrate<Ode, false>(ode);

      return res;
    }
  };

  /*
   * specialization for @c PkAdams
   */
  template<>
  struct PMXOdeIntegrator<PkRk45> {
    const double rtol;
    const double atol;
    const long int max_num_step;
    std::ostream* msgs;

    PMXOdeIntegrator() : rtol(1e-10), atol(1e-10), max_num_step(1e8), msgs(0) {}

    PMXOdeIntegrator(const double rtol0, const double atol0, const long int max_num_step0,
                     std::ostream* msgs0) :
      rtol(rtol0), atol(atol0), max_num_step(max_num_step0), msgs(msgs0)
    {}

    /*
     * Torsten's ODE solvers support @c T_t to be @c var
     */
    template <typename F, typename Tt, typename T_initial, typename T_param>
    std::vector<std::vector<typename stan::return_type<Tt, T_initial, T_param>::type> > // NOLINT
    operator()(const F& f,
               const std::vector<T_initial>& y0,
               double t0,
               const std::vector<Tt>& ts,
               const std::vector<T_param>& theta,
               const std::vector<double>& x_r,
               const std::vector<int>& x_i) const {
      return torsten::pmx_integrate_ode_rk45(f, y0, t0, ts, theta, x_r, x_i, msgs, rtol, atol, max_num_step);
    }

    /*
     * For MPI solution we need to return data that consists
     * of solution value and gradients
     */
    template <typename F, typename Tt, typename T_initial, typename T_param>
    Eigen::MatrixXd
    solve_d(const F& f,
            const std::vector<T_initial>& y0,
            double t0,
            const std::vector<Tt>& ts,
            const std::vector<T_param>& theta,
            const std::vector<double>& x_r,
            const std::vector<int>& x_i) const {
      using Ode = torsten::dsolve::PMXOdeintSystem<F, Tt, T_initial, T_param>;
      const int m = theta.size();
      const int n = y0.size();

      dsolve::PMXOdeService<Ode, dsolve::Odeint> serv(n, m);

      Ode ode{serv, f, t0, ts, y0, theta, x_r, x_i, msgs};
      using scheme_t = boost::numeric::odeint::runge_kutta_dopri5<std::vector<double>, double, std::vector<double>, double>;
      dsolve::PMXOdeintIntegrator<scheme_t> solver(rtol, atol, max_num_step);
      Eigen::MatrixXd res = solver.integrate<Ode, false>(ode);

      return res;
    }
  };
}

#endif
