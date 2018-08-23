#ifndef PK_ODE_SOLVER_HPP
#define PK_ODE_SOLVER_HPP

#include <stan/math/torsten/pk_ode_integrator.hpp>

namespace refactor {

  using boost::math::tools::promote_args;
  using Eigen::Matrix;
  using Eigen::Dynamic;
  using stan::math::to_array_1d;
  using std::vector;
  using boost::math::tools::promote_args;

  /**
   * ODE-based PKPD model solver. It selects the integrator
   * type and applies it onto ODE models. The model's rate
   * is assumed to be data and passed to @c x_r. If the model's
   * rate is parameter, one must first applies @c PKODERateAdaptor
   * and pass it to @c PKODEModelSolver. Most of the time,
   * one should simply use @c PKODERateAdaptor on any ODE
   * model that is to be solved by @c PKODEModelSolver.
   *
   * @tparam F_i integrator function type
   */
  template<PkOdeIntegratorId It>
  class PKODEModelSolver {
    const double rtol_;
    const double atol_;
    const int max_num_steps_;
    std::ostream* msgs_;
  public:
    // constructor
    PKODEModelSolver(double rel_tol,
                     double abs_tol,
                     long int max_num_steps,
                     std::ostream* msgs) :
      rtol_(rel_tol),
      atol_(abs_tol),
      max_num_steps_(max_num_steps),
      msgs_(msgs)
    {}

  /**
   * ODE-based PKPD model solvers.
   *
   * @tparam T_time dt type
   * @tparam T_model ODE model type
   * @param pkmodel ODE-based model
   * @param dt time span
   * @return col vector of ODE solution
   */
    template<typename T_time, typename T_model,
             typename std::enable_if_t<!torsten::has_var_rate<T_model>::value>* = nullptr> // NOLINT
    Eigen::Matrix<torsten::scalar_t<T_model>, Eigen::Dynamic, 1> 
    solve(const T_model &pkmodel, const T_time& dt) const {
      using T_scalar = torsten::scalar_t<T_model>;
      Eigen::Matrix<T_scalar, Eigen::Dynamic, 1> res;
      auto t0 = pkmodel.t0();
      std::vector<T_time> ts{t0 + dt};
      if (ts[0] == t0) {
        res = pkmodel.y0();
      } else {
        auto y = stan::math::to_array_1d(pkmodel.y0());
        auto par = pkmodel.par();
        auto rate = pkmodel.rate();
        std::vector<int> x_i;
        std::vector<std::vector<T_scalar> > res_v =
          PkOdeIntegrator<It>()(pkmodel.f(), y, t0, ts, par, 
                                rate, x_i, msgs_, rtol_, atol_, max_num_steps_);
        res = stan::math::to_vector(res_v[0]);
      }
      return res;
    }

    template<typename T_time, typename T_model,
             typename std::enable_if_t<torsten::has_var_rate<T_model>::value>* = nullptr> // NOLINT
    Eigen::Matrix<torsten::scalar_t<T_model>, Eigen::Dynamic, 1> 
    solve(const T_model &pkmodel, const T_time& dt) const {
      using T_scalar = torsten::scalar_t<T_model>;
      Eigen::Matrix<T_scalar, Eigen::Dynamic, 1> res;
      auto t0 = pkmodel.t0();
      std::vector<T_time> ts{t0 + dt};
      if (ts[0] == t0) {
        res = pkmodel.y0();
      } else {
        auto y = stan::math::to_array_1d(pkmodel.y0());
        auto par = pkmodel.par();
        auto rate = pkmodel.rate();
        std::vector<double> x_r;
        std::vector<int> x_i;
        std::vector<std::vector<T_scalar> > res_v =
          PkOdeIntegrator<It>()(pkmodel.f(), y, t0, ts, par, 
                                x_r, x_i, msgs_, rtol_, atol_, max_num_steps_);
        res = stan::math::to_vector(res_v[0]);
      }
      return res;
    }

    // // WIP, for SS solver
    // template<typename F>
    // Eigen::Matrix<double, Eigen::Dynamic, 1> 
    // solve(const F& f,
    //       const double &t0,
    //       const Eigen::Matrix<double, 1, Eigen::Dynamic>& init,
    //       const double& dt,
    //       const std::vector<double>& pars,
    //       const std::vector<double>& rate){
    //   using stan::math::to_array_1d;
    //   using std::vector;

    //   std::vector<double> EventTime {t0};

    //   assert((size_t) init.cols() == rate.size());

    //   double InitTime = t0 - dt;  // time of previous event

    //   Eigen::Matrix<double, Eigen::Dynamic, 1> pred;
      
    //   run_integrator(f,
    //                  pars,
    //                  rate,
    //                  init,
    //                  InitTime,
    //                  EventTime,
    //                  pred);
    //   return pred;
    // }
  };
}

#endif
