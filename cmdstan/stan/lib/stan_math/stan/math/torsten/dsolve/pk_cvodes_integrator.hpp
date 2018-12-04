#ifndef STAN_MATH_TORSTEN_DSOLVE_CVODES_INTEGRATOR_HPP
#define STAN_MATH_TORSTEN_DSOLVE_CVODES_INTEGRATOR_HPP

#include <stan/math/rev/scal/meta/is_var.hpp>
#include <stan/math/torsten/dsolve/pk_cvodes_fwd_system.hpp>
#include <stan/math/torsten/dsolve/sundials_check.hpp>
#include <cvodes/cvodes.h>
#include <cvodes/cvodes_direct.h>
#include <sunmatrix/sunmatrix_dense.h>
#include <sunlinsol/sunlinsol_dense.h>
#include <nvector/nvector_serial.h>
#include <ostream>
#include <vector>
#include <algorithm>

namespace torsten {
namespace dsolve {
/**
 * CVODES ODE integrator.
 */
  class PKCvodesIntegrator {
    const double rtol_;
    const double atol_;
    const int64_t max_num_steps_;

    // /**
    //  *  ode adjoint sens calculation requires different initialization
    //  *
    //  * @tparam F type of ODE RHS functor
    //  * @tparam Ty0 type of ODE primary unknowns
    //  * @tparam Tpar type of ODE parameters.
    //  * @param[out] ode ODE system
    //  * @param[out] res_y ODE solutions
    //  */
    // template <typename F, typename Ty0, typename Tpar>
    //   // TODO(yizhang): adjoint sensitivity initialization
    // }

    template <typename F, int Lmm>
    void solve(PKCvodesFwdSystem<F, double, double, double, Lmm>& ode,
               std::vector<std::vector<double> >& res_y);

    template <typename F, int Lmm>
    void solve(PKCvodesFwdSystem<F, stan::math::var,
               double, double, Lmm>& ode,
               std::vector<std::vector<stan::math::var> >& res_y);

    template <typename F, typename Ty0, typename Tpar, int Lmm>
    void solve(PKCvodesFwdSystem<F,
               double, Ty0, Tpar, Lmm>& ode,
               std::vector<
               std::vector<stan::math::var> >& res_y);

    template <typename F, typename Ty0, typename Tpar, int Lmm>
    void solve(PKCvodesFwdSystem<F, stan::math::var,
               Ty0, Tpar, Lmm>& ode,
               std::vector<
               std::vector<stan::math::var> >& res_y);

    // TODO(yizhang): adjoint sensitivity solver

  public:
    static constexpr int CVODES_MAX_STEPS = 500;

    /**
     * constructor
     * @param[in] rtol relative tolerance
     * @param[in] atol absolute tolerance
     * @param[in] max_num_steps max nb. of times steps
     */
    PKCvodesIntegrator(const double rtol, const double atol,
                         const int64_t max_num_steps = CVODES_MAX_STEPS)
      : rtol_(rtol), atol_(atol), max_num_steps_(max_num_steps) {
      using stan::math::invalid_argument;
      if (rtol_ <= 0)
        invalid_argument("cvodes_integrator", "relative tolerance,", rtol_, "",
                         ", must be greater than 0");
      if (rtol_ > 1.0E-3)
        invalid_argument("cvodes_integrator", "relative tolerance,", rtol_, "",
                         ", must be less than 1.0E-3");
      if (atol_ <= 0)
        invalid_argument("cvodes_integrator", "absolute tolerance,", atol_, "",
                         ", must be greater than 0");
      if (max_num_steps_ <= 0)
        invalid_argument("cvodes_integrator", "max_num_steps,",
                         max_num_steps_, "",
                         ", must be greater than 0");
    }

    /**
     * Return the solutions for the specified ODE
     * given the specified initial state,
     * initial times, times of desired solution, and parameters and
     * data, writing error and warning messages to the specified
     * stream contained in the ODE system.
     *
     * @tparam ODE type of ODE system
     * @param[in] ode ODE system
     * increasing order, all greater than the initial time.
     * @return a vector of states, each state being a vector of the
     * same size as the state variable, corresponding to a time in ts.
     */
    template <typename Ode>
    typename Ode::return_type integrate(Ode& ode) {
      using Eigen::Dynamic;
      using Eigen::Matrix;
      using Eigen::MatrixXd;
      using Eigen::VectorXd;
      using stan::math::check_finite;
      using stan::math::check_ordered;
      using stan::math::check_nonzero_size;
      using stan::math::check_less;

      auto mem       = ode.mem();
      auto y         = ode.nv_y();
      auto y0        = ode.y0();
      auto ys        = ode.nv_ys();
      const size_t n = ode.n();
      const size_t ns= ode.ns();

      typename Ode::return_type
        res_y(ode.ts().size(), std::vector<typename Ode::scalar_type>(n, 0));

      // Initial condition is from nv_y, which has changed
      // from previous solution, we need to reset it. Likewise
      // we also reset ys.
      for (size_t i = 0; i < n; ++i) {
        NV_Ith_S(y, i) = stan::math::value_of(y0[i]);
      }
      for (size_t is = 0; is < ns; ++is) {
        N_VConst(RCONST(0.0), ys[is]);
      }

      try {
        CHECK_SUNDIALS_CALL(CVodeReInit(mem, ode.t0(), y));
        CHECK_SUNDIALS_CALL(CVodeSStolerances(mem, rtol_, atol_));
        CHECK_SUNDIALS_CALL(CVodeSetUserData(mem, ode.to_user_data()));
        CHECK_SUNDIALS_CALL(CVodeSetMaxNumSteps(mem, max_num_steps_));

        /** if y0 is parameter, the first n sensitivity vector
         * are regarding y0, thus they form a unit matrix.
         **/
        if (Ode::need_fwd_sens) {
          if (Ode::is_var_y0)
            for (size_t i = 0; i < n; ++i) NV_Ith_S(ys[i], i) = 1.0;
          CHECK_SUNDIALS_CALL(CVodeSensInit(mem, ode.ns(), CV_SIMULTANEOUS,
                                            ode.sens_rhs(), ys));
          CHECK_SUNDIALS_CALL(CVodeSensEEtolerances(mem));
        }

        solve(ode, res_y);
      } catch (const std::exception& e) {
        CVodeSensFree(mem);
        throw;
      }

      CVodeSensFree(mem);

      return res_y;
    }
  };  // cvodes integrator

  /**
   * Solve ODE system, no sensitivty
   *
   * @tparam F ODE functor type
   * @param[out] ode ODE system
   * @param[out] res_y ODE solutions
   */
  template <typename F, int Lmm>
  void PKCvodesIntegrator::solve(PKCvodesFwdSystem<F, double,
                                   double, double, Lmm>& ode,
                                   std::vector<std::vector<double> >& res_y) {
    double t1 = ode.t0();
    const std::vector<double>& ts = ode.ts();
    auto mem = ode.mem();
    auto y = ode.nv_y();

    for (size_t i = 0; i < ts.size(); ++i) {
      CHECK_SUNDIALS_CALL(CVode(mem, ts[i], y, &t1, CV_NORMAL));
      for (size_t j = 0; j < ode.n(); ++j) res_y[i][j] = NV_Ith_S(y, j);
    };
  }

  /**
   * Solve ODE system, with only ts sensitivity. In this
   * case we don't have to solve a fwd sensitivty system. As
   * the sensitivity regarding @c ts[i] is just the RHS of ODE
   * for solution at @c ts[i], and is 0 for solution at
   * other time. Since @c fval stores the last RHS
   * evaluation in @c eval_f, we can just use it instead
   * recompute RHS.a
   *
   * @tparam F ODE functor type
   * @param[out] ode ODE system
   * @param[out] res_y ODE solutions
   */
  template <typename F, int Lmm>
  void PKCvodesIntegrator::solve(PKCvodesFwdSystem<F,
                                   stan::math::var, double, double, Lmm>& ode,
                                   std::vector<
                                   std::vector<stan::math::var> >& res_y) {
    using stan::math::value_of;
    using stan::math::var;
    double t1 = ode.t0();
    const std::vector<var>& ts = ode.ts();
    auto mem = ode.mem();
    auto y = ode.nv_y();

    /* 
     * we use @c y_vec for workspace of gradient, most of
     * which are zero when @c ts is the only parameter. We
     * only update the entry that is non-zero when create @c
     * var before overwrite it with zero again.
     */
    static std::vector<double> g(ts.size(), 0.0);

    for (size_t i = 0; i < ts.size(); ++i) {
      double time = value_of(ts[i]);
      CHECK_SUNDIALS_CALL(CVode(mem, time, y, &t1, CV_NORMAL));
      for (size_t j = 0; j < ode.n(); ++j) {
        ode.eval_rhs(time, y);
        g[i] = ode.fval()[j];
        res_y[i][j] = precomputed_gradients(NV_Ith_S(y, j), ts, g);
        g[i] = 0.0;
      }
    }
  }

  /**
   * Solve Ode system with forward sensitivty, return a
   * vector of var with precomputed gradient as sensitivity
   * value, when @c y0 and/or @c theta are parameters but @c ts is not.
   *
   * @tparam Ode ODE system type
   * @param[out] ode ODE system
   * @param[out] res_y ODE solutions
   */
  template <typename F, typename Ty0, typename Tpar, int Lmm>
  void PKCvodesIntegrator::solve(PKCvodesFwdSystem<F,
                                   double, Ty0, Tpar, Lmm>& ode,
                                   std::vector<
                                   std::vector<stan::math::var> >& res_y) {
    using stan::math::precomputed_gradients;
    double t1 = ode.t0();
    const std::vector<double>& ts = ode.ts();
    auto mem = ode.mem();
    auto y = ode.nv_y();
    auto ys = ode.nv_ys();
    const auto n = ode.n();
    const auto ns = ode.ns();
    auto vars = ode.vars();

    static std::vector<double> g(ns);

    for (size_t i = 0; i < ts.size(); ++i) {
      CHECK_SUNDIALS_CALL(CVode(mem, ts[i], y, &t1, CV_NORMAL));
      if (ode.need_fwd_sens) {
        CHECK_SUNDIALS_CALL(CVodeGetSens(mem, &t1, ys));
        for (size_t k = 0; k < n; ++k) {
          for (size_t j = 0; j < ns; ++j) {
            g[j] = NV_Ith_S(ys[j], k);
          }
          res_y[i][k] = precomputed_gradients(NV_Ith_S(y, k), vars, g);
        }
      }
    }
  }

  /**
   * Solve Ode system with forward sensitivty, return a
   * vector of var with precomputed gradient as sensitivity
   * value, when @c ts, @c theta and/or @c y0 are parameters.
   *
   * @tparam Ode ODE system type
   * @param[out] ode ODE system
   * @param[out] res_y ODE solutions
   */
  template <typename F, typename Ty0, typename Tpar, int Lmm>
  void PKCvodesIntegrator::solve(PKCvodesFwdSystem<F,
                                   stan::math::var,
                                   Ty0, Tpar, Lmm>& ode,
                                   std::vector<
                                   std::vector<stan::math::var> >& res_y) {
    using stan::math::precomputed_gradients;
    using stan::math::var;
    double t1 = ode.t0();
    const std::vector<var>& ts = ode.ts();
    auto mem = ode.mem();
    auto y = ode.nv_y();
    auto ys = ode.nv_ys();
    const auto n = ode.n();
    const auto ns = ode.ns();
    auto vars = ode.vars();

    static std::vector<double> g(ns + ts.size(), 0.0);

    for (size_t i = 0; i < ts.size(); ++i) {
      double time = value_of(ts[i]);
      CHECK_SUNDIALS_CALL(CVode(mem, time, y, &t1, CV_NORMAL));
      if (ode.need_fwd_sens) {
        CHECK_SUNDIALS_CALL(CVodeGetSens(mem, &t1, ys));
        for (size_t k = 0; k < n; ++k) {
          for (size_t j = 0; j < ns; ++j) g[j] = NV_Ith_S(ys[j], k);
          ode.eval_rhs(time, y);
          g[i + ns] = ode.fval()[k];
          res_y[i][k] = precomputed_gradients(NV_Ith_S(y, k), vars, g);
          g[i + ns] = 0.0;
        }
      }
    }
  }

}  // namespace dsolve
}  // namespace torsten

#endif
