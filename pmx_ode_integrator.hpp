#ifndef STAN_MATH_TORSTEN_PMX_ODE_INTEGRATOR_HPP
#define STAN_MATH_TORSTEN_PMX_ODE_INTEGRATOR_HPP

#include <stan/math/prim/fun/eval.hpp>
#include <stan/math/rev/functor/jacobian.hpp>
#include <stan/math/prim/meta/return_type.hpp>
#include <stan/math/prim/fun/value_of.hpp>
#include <stan/math/prim/functor/integrate_ode_rk45.hpp>
#include <stan/math/rev/functor/integrate_ode_adams.hpp>
#include <stan/math/rev/functor/integrate_ode_bdf.hpp>
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
  using stan::math::integrate_ode_adams;
  using stan::math::integrate_ode_bdf;
  using stan::math::integrate_ode_rk45;
  using torsten::pmx_integrate_ode_adams;
  using torsten::pmx_integrate_ode_bdf;
  using torsten::pmx_integrate_ode_rk45;
  using torsten::dsolve::PMXCvodesFwdSystem;
  using torsten::dsolve::cvodes_def;

#define DEF_STAN_INTEGRATOR(INT_NAME)                                                 \
  template <typename F, typename Tt, typename T_initial, typename T_param>            \
  std::vector<std::vector<typename stan::return_type_t<Tt, T_initial, T_param>> > \
  operator()(const F& f,                                                              \
             const std::vector<T_initial>& y0,                                        \
             double t0,                                                               \
             const std::vector<Tt>& ts,                                               \
             const std::vector<T_param>& theta,                                       \
             const std::vector<double>& x_r,                                          \
             const std::vector<int>& x_i) const {                                     \
    std::vector<double> ts_dbl(stan::math::value_of(ts));                             \
    return INT_NAME(f, y0, t0, ts, theta, x_r, x_i, msgs, rtol, atol, max_num_step);  \
  }

#define DEF_STAN_SINGLE_STEP_INTEGRATOR                                               \
  template <typename F, typename Tt, typename T_initial, typename T_param>            \
  std::vector<std::vector<typename stan::return_type_t<Tt, T_initial, T_param>> > \
  operator()(const F& f,                                                              \
             const std::vector<T_initial>& y0,                                        \
             double t0,                                                               \
             const Tt& t1,                                                            \
             const std::vector<T_param>& theta,                                       \
             const std::vector<double>& x_r,                                          \
             const std::vector<int>& x_i) const {                                     \
    std::vector<double> ts{stan::math::value_of(t1)};                                 \
    return (*this)(f, y0, t0, ts, theta, x_r, x_i);                                   \
  }

#define DEF_TORSTEN_INTEGRATOR(INT_NAME)                                              \
  template <typename F, typename Tt, typename T_initial, typename T_param>            \
  std::vector<std::vector<typename stan::return_type_t<Tt, T_initial, T_param>> > \
  operator()(const F& f,                                                              \
             const std::vector<T_initial>& y0,                                        \
             double t0,                                                               \
             const std::vector<Tt>& ts,                                               \
             const std::vector<T_param>& theta,                                       \
             const std::vector<double>& x_r,                                          \
             const std::vector<int>& x_i) const {                                     \
    return INT_NAME(f, y0, t0, ts, theta, x_r, x_i, rtol, atol, max_num_step, msgs);  \
  }

#define DEF_TORSTEN_SINGLE_STEP_INTEGRATOR                                            \
  template <typename F, typename Tt, typename T_initial, typename T_param>            \
  std::vector<std::vector<typename stan::return_type_t<Tt, T_initial, T_param>> > \
  operator()(const F& f,                                                              \
             const std::vector<T_initial>& y0,                                        \
             double t0,                                                               \
             const Tt& t1,                                                            \
             const std::vector<T_param>& theta,                                       \
             const std::vector<double>& x_r,                                          \
             const std::vector<int>& x_i) const {                                     \
    std::vector<Tt> ts{t1};                                                           \
    return (*this)(f, y0, t0, ts, theta, x_r, x_i);                                   \
  }

#define DEF_TORSTEN_SINGLE_STEP_SOLVE_D                                      \
    template <typename F, typename Tt, typename T_initial, typename T_param> \
    Eigen::MatrixXd                                                          \
    solve_d(const F& f,                                                      \
            const std::vector<T_initial>& y0,                                \
            double t0,                                                       \
            const Tt& t1,                                                    \
            const std::vector<T_param>& theta,                               \
            const std::vector<double>& x_r,                                  \
            const std::vector<int>& x_i) const {                             \
      std::vector<Tt> ts{t1};                                                \
      return this -> solve_d(f, y0, t0, ts, theta, x_r, x_i);                \
    }

  template<PMXOdeIntegratorId It>
  struct PMXOdeIntegrator;

#define DEF_TORSTEN_INTEGRATOR_MEMBER(RTOL, ATOL, MAXSTEP, AS_RTOL, AS_ATOL, AS_MAXSTEP)  \
    const double rtol;                                                                    \
    const double atol;                                                                    \
    const long int max_num_step;                                                          \
    const double as_rtol;                                                                 \
    const double as_atol;                                                                 \
    const long int as_max_num_step;                                                       \
    std::ostream* msgs;                                                                   \
    PMXOdeIntegrator() : rtol(RTOL), atol(ATOL), max_num_step(MAXSTEP),                   \
                         as_rtol(AS_RTOL), as_atol(AS_ATOL), as_max_num_step(AS_MAXSTEP), \
                         msgs(0) {}                                                       \
    PMXOdeIntegrator(double rtol0, double atol0, long int max_num_step0,                  \
                     std::ostream* msgs0) :                                               \
      rtol(rtol0), atol(atol0), max_num_step(max_num_step0),                              \
      as_rtol(AS_RTOL), as_atol(AS_ATOL), as_max_num_step(AS_MAXSTEP),                    \
      msgs(msgs0) {}                                                                      \
    PMXOdeIntegrator(double rtol0, double atol0, long int max_num_step0,                  \
                     double as_rtol0, double as_atol0, long int as_max_num_step0,         \
                     std::ostream* msgs0) :                                               \
      rtol(rtol0), atol(atol0), max_num_step(max_num_step0),                              \
      as_rtol(as_rtol0), as_atol(as_atol0), as_max_num_step(as_max_num_step0),            \
      msgs(msgs0) {}


  template<>
  struct PMXOdeIntegrator<StanAdams> {
    DEF_TORSTEN_INTEGRATOR_MEMBER(1e-10, 1e-10, 1e8, 1e-6, 1e-6, 1e2)
    DEF_STAN_INTEGRATOR(integrate_ode_adams)
    DEF_STAN_SINGLE_STEP_INTEGRATOR
  };

  template<>
  struct PMXOdeIntegrator<StanBdf> {
    DEF_TORSTEN_INTEGRATOR_MEMBER(1e-10, 1e-10, 1e8, 1e-6, 1e-6, 1e2)
    DEF_STAN_INTEGRATOR(integrate_ode_bdf)
    DEF_STAN_SINGLE_STEP_INTEGRATOR
  };

  template<>
  struct PMXOdeIntegrator<StanRk45> {
    DEF_TORSTEN_INTEGRATOR_MEMBER(1e-10, 1e-10, 1e8, 1e-6, 1e-6, 1e2)
    DEF_STAN_INTEGRATOR(integrate_ode_rk45)
    DEF_STAN_SINGLE_STEP_INTEGRATOR
  };

  /**
   * specialization to a dummy integrator used for
   * analytical solutions such as one-cpt & two-cpt
   */
  template<>
  struct PMXOdeIntegrator<Analytical> {};

  /*
   * specialization for @c PkBdf
   */
  template<>
  struct PMXOdeIntegrator<PkBdf> {
    DEF_TORSTEN_INTEGRATOR_MEMBER(1e-10, 1e-10, 1e8, 1e-6, 1e-6, 1e2)
    DEF_TORSTEN_INTEGRATOR(pmx_integrate_ode_bdf)
    DEF_TORSTEN_SINGLE_STEP_INTEGRATOR

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
      using Ode = PMXCvodesFwdSystem<F, Tt, T_initial, T_param, cvodes_def<TORSTEN_CV_SENS, CV_BDF, TORSTEN_CV_ISM>>;
      const int m = theta.size();
      const int n = y0.size();

      dsolve::PMXOdeService<Ode> serv(n, m);

      Ode ode{serv, f, t0, ts, y0, theta, x_r, x_i, msgs};

      torsten::dsolve::PMXCvodesIntegrator solver(rtol, atol, max_num_step);
      Eigen::MatrixXd res = solver.integrate<Ode, false>(ode);

      return res;
    }

    DEF_TORSTEN_SINGLE_STEP_SOLVE_D
  };

  /*
   * specialization for @c PkAdams
   */
  template<>
  struct PMXOdeIntegrator<PkAdams> {
    DEF_TORSTEN_INTEGRATOR_MEMBER(1e-10, 1e-10, 1e8, 1e-6, 1e-6, 1e2)
    DEF_TORSTEN_INTEGRATOR(pmx_integrate_ode_adams)
    DEF_TORSTEN_SINGLE_STEP_INTEGRATOR

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
      using Ode = PMXCvodesFwdSystem<F, Tt, T_initial, T_param, cvodes_def<TORSTEN_CV_SENS, CV_ADAMS, TORSTEN_CV_ISM>>;
      const int m = theta.size();
      const int n = y0.size();

      dsolve::PMXOdeService<Ode> serv(n, m);

      Ode ode{serv, f, t0, ts, y0, theta, x_r, x_i, msgs};

      torsten::dsolve::PMXCvodesIntegrator solver(rtol, atol, max_num_step);
      Eigen::MatrixXd res = solver.integrate<Ode, false>(ode);

      return res;
    }

    DEF_TORSTEN_SINGLE_STEP_SOLVE_D
  };

  /*
   * specialization for @c PkAdams
   */
  template<>
  struct PMXOdeIntegrator<PkRk45> {
    DEF_TORSTEN_INTEGRATOR_MEMBER(1e-10, 1e-10, 1e8, 1e-6, 1e-6, 1e2)
    DEF_TORSTEN_INTEGRATOR(pmx_integrate_ode_rk45)
    DEF_TORSTEN_SINGLE_STEP_INTEGRATOR

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

    DEF_TORSTEN_SINGLE_STEP_SOLVE_D
  };
}

#undef DEF_TORSTEN_INTEGRATOR_MEMBER
#undef DEF_TORSTEN_SINGLE_STEP_SOLVE_D
#undef DEF_TORSTEN_SINGLE_STEP_INTEGRATOR
#undef DEF_TORSTEN_INTEGRATOR
#undef DEF_STAN_SINGLE_STEP_INTEGRATOR
#undef DEF_STAN_INTEGRATOR

#endif
