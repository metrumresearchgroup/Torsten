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
    PkAdams, PkBdf
  };
}

namespace torsten {
  namespace internal {
    template<PMXOdeIntegratorId It>
    struct PMXOdeIntegratorDispatcher;

    /* 
     * specification for a Stan's ODE integrator. Since @c
     * ts cannot be param in Stan's integrators, we convert
     * it to data first.
     */
    template<>
    struct PMXOdeIntegratorDispatcher<StanRk45> {

      template <typename F, typename Tt, typename T_initial, typename T_param>
      std::vector<std::vector<typename stan::return_type<T_initial, T_param>::type> > // NOLINT
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
        const std::vector<double>& ts_d(stan::math::value_of(ts));
        return stan::math::integrate_ode_rk45(f, y0, t0, ts_d, theta,
                                              x_r, x_i, msgs, rtol,
                                              atol, max_num_step);
      }
    };

    /* 
     * specification for a Stan's ODE integrator. Since @c
     * ts cannot be param in Stan's integrators, we convert
     * it to data first.
     */
    template<>
    struct PMXOdeIntegratorDispatcher<StanAdams> {
      template <typename F, typename Tt, typename T_initial, typename T_param>
      std::vector<std::vector<typename stan::return_type<T_initial,
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
        const std::vector<double>& ts_d(stan::math::value_of(ts));
        return stan::math::integrate_ode_adams(f, y0, t0, ts_d, theta,
                                               x_r, x_i, msgs, rtol,
                                               atol, max_num_step);
      }
    };

    /* 
     * specification for a Stan's ODE integrator. Since @c
     * ts cannot be param in Stan's integrators, we convert
     * it to data first.
     */
    template<>
    struct PMXOdeIntegratorDispatcher<StanBdf> {
      template <typename F, typename Tt, typename T_initial, typename T_param>
      std::vector<std::vector<typename stan::return_type<T_initial,
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
        const std::vector<double>& ts_d(stan::math::value_of(ts));
        return stan::math::integrate_ode_bdf(f, y0, t0, ts_d, theta,
                                             x_r, x_i, msgs, rtol,
                                             atol, max_num_step);
      }
    };

    /* 
     * specification for a Torsten's ODE integrator, in
     * which @c ts can be a param.
     */
    template<>
    struct PMXOdeIntegratorDispatcher<PkAdams> {
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
        return torsten::dsolve::pmx_integrate_ode_adams(f, y0, t0, ts, theta,
                                                       x_r, x_i, msgs, rtol,
                                                       atol, max_num_step);
      }
    };

    /* 
     * specification for a Torsten's ODE integrator, in
     * which @c ts can be a param.
     */
    template<>
    struct PMXOdeIntegratorDispatcher<PkBdf> {
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
        return torsten::dsolve::pmx_integrate_ode_bdf(f, y0, t0, ts, theta,
                                                     x_r, x_i, msgs, rtol,
                                                     atol, max_num_step);
      }
    };    
  }
}

namespace torsten {
  template<PMXOdeIntegratorId It>
  struct PMXOdeIntegrator {
    const double rtol;
    const double atol;
    const long int max_num_step;
    std::ostream* msgs;

    PMXOdeIntegrator() : rtol(1e-10), atol(1e-10), max_num_step(1e8), msgs(0) {}

    PMXOdeIntegrator(const double rtol0, const double atol0,
                    const long int max_num_step0,
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
      using internal::PMXOdeIntegratorDispatcher;
      return PMXOdeIntegratorDispatcher<It>()(f, y0, t0, ts, theta, x_r, x_i,
                                             msgs, rtol, atol, max_num_step);
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

    PMXOdeIntegrator(const double rtol0, const double atol0,
                    const long int max_num_step0,
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
      using internal::PMXOdeIntegratorDispatcher;
      return PMXOdeIntegratorDispatcher<PkBdf>()(f, y0, t0, ts, theta, x_r, x_i,
                                                msgs, rtol, atol, max_num_step);
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

      dsolve::PMXCvodesService<typename Ode::Ode> serv(n, m);

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

    PMXOdeIntegrator(const double rtol0, const double atol0,
                    const long int max_num_step0,
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
      using internal::PMXOdeIntegratorDispatcher;
      return PMXOdeIntegratorDispatcher<PkAdams>()(f, y0, t0, ts, theta, x_r, x_i,
                                                msgs, rtol, atol, max_num_step);
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

      dsolve::PMXCvodesService<typename Ode::Ode> serv(n, m);

      Ode ode{serv, f, t0, ts, y0, theta, x_r, x_i, msgs};

      torsten::dsolve::PMXCvodesIntegrator solver(rtol, atol, max_num_step);
      Eigen::MatrixXd res = solver.integrate<Ode, false>(ode);

      return res;
    }
  };
}

#endif
