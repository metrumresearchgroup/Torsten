#ifndef STAN_MATH_TORSTEN_DSOLVE_CVODES_INTEGRATOR_HPP
#define STAN_MATH_TORSTEN_DSOLVE_CVODES_INTEGRATOR_HPP

#include <stan/math/torsten/dsolve/pmx_cvodes_fwd_system.hpp>
#include <stan/math/torsten/dsolve/cvodes_jac.hpp>
#include <type_traits>

namespace torsten {
namespace dsolve {

/**
 * CVODES ODE integrator.
 */
  struct PMXCvodesIntegrator {
    const double rtol_;
    const double atol_;
    const int64_t max_num_steps_;

    /*
     * Observer stores return data that is updated as requested.
     * Usually the returned is of @c array_2d @c var type.
     */
    template<typename Ode, bool GenVar>
    struct SolObserver {
      std::vector<std::vector<typename Ode::scalar_type>> y;

      SolObserver(const Ode& ode) :
        y{ode.ts().size(), std::vector<typename Ode::scalar_type>(ode.n(), 0.0)}
      {}
    };

    /*
     * Observer stores return data that is updated as requested.
     * For MPI results we need return data type so it can be
     * sent over to other nodes before reassembled into @c var
     * type. In this case, the returned matrix contain value
     * and sensitivity.
     */
    template<typename Ode>
    struct SolObserver<Ode, false> {
      Eigen::MatrixXd y;

      SolObserver(const Ode& ode) :
        y(Eigen::MatrixXd::Zero(ode.fwd_system_size() +
                                ode.n() * (Ode::is_var_ts ? ode.ts().size() : 0),
                                ode.ts().size()))
      {}
    };

    // Seqential
    template <typename F, int Lmm, PMXCvodesSensMethod Sm>
    void solve(PMXCvodesFwdSystem<F, double, double, double, Lmm, Sm>& ode,
               std::vector<std::vector<double> >& res_y);

    template <typename F, int Lmm, PMXCvodesSensMethod Sm>
    void solve(PMXCvodesFwdSystem<F, stan::math::var, double, double, Lmm, Sm>& ode,
               std::vector<std::vector<stan::math::var> >& res_y);

    template <typename F, typename Ty0, int Lmm, PMXCvodesSensMethod Sm>
    void solve(PMXCvodesFwdSystem<F, double, Ty0, double, Lmm, Sm>& ode,
               std::vector<std::vector<stan::math::var> >& res_y);

    template <typename F, typename Tpar, int Lmm, PMXCvodesSensMethod Sm>
    void solve(PMXCvodesFwdSystem<F, double, double, Tpar, Lmm, Sm>& ode,
               std::vector<std::vector<stan::math::var> >& res_y);

    template <typename F, typename Ty0, typename Tpar, int Lmm, PMXCvodesSensMethod Sm>
    void solve(PMXCvodesFwdSystem<F, double, Ty0, Tpar, Lmm, Sm>& ode,
               std::vector<std::vector<stan::math::var> >& res_y);

    template <typename F, typename Ty0, typename Tpar, int Lmm, PMXCvodesSensMethod Sm>
    void solve(PMXCvodesFwdSystem<F, stan::math::var, Ty0, Tpar, Lmm, Sm>& ode,
               std::vector<std::vector<stan::math::var> >& res_y);

    // MPI
    template <typename F, int Lmm, PMXCvodesSensMethod Sm>
    void solve(PMXCvodesFwdSystem<F, double, double, double, Lmm, Sm>& ode,
               Eigen::MatrixXd& res_y);

    template <typename F, int Lmm, PMXCvodesSensMethod Sm>
    void solve(PMXCvodesFwdSystem<F, stan::math::var, double, double, Lmm, Sm>& ode,
               Eigen::MatrixXd& res_y);

    template <typename F, typename Ty0, typename Tpar, int Lmm, PMXCvodesSensMethod Sm>
    void solve(PMXCvodesFwdSystem<F, double, Ty0, Tpar, Lmm, Sm>& ode,
               Eigen::MatrixXd& res_y);

    template <typename F, typename Ty0, typename Tpar, int Lmm, PMXCvodesSensMethod Sm>
    void solve(PMXCvodesFwdSystem<F, stan::math::var, Ty0, Tpar, Lmm, Sm>& ode, // NOLINT
                                   Eigen::MatrixXd& res_y);

    // TODO
    // template <typename F, typename Ty0, typename Tpar, int Lmm, PMXCvodesSensMethod Sm>
    // void solve(PMXCvodesFwdSystem<F, stan::math::var, Ty0, Tpar, Lmm, Sm>& ode,
    //            Eigen::MatrixXd& res_y);

    // TODO(yizhang): adjoint sensitivity solver

  public:
    static constexpr int CVODES_MAX_STEPS = 500;

    /**
     * constructor
     * @param[in] rtol relative tolerance
     * @param[in] atol absolute tolerance
     * @param[in] max_num_steps max nb. of times steps
     */
    PMXCvodesIntegrator(const double rtol, const double atol,
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
    template <typename Ode, bool GenVar = true>
    auto integrate(Ode& ode) {
      using std::vector;
      using Eigen::Dynamic;
      using Eigen::Matrix;
      using Eigen::MatrixXd;

      auto mem       = ode.mem();
      auto y         = ode.nv_y();
      auto ys        = ode.nv_ys();
      const size_t n = ode.n();
      const size_t ns= ode.ns();

      SolObserver<Ode, GenVar> observer(ode);

      // Initial condition is from nv_y, which has changed
      // from previous solution, we need to reset it. Likewise
      // we also reset ys.
      for (size_t i = 0; i < n; ++i) {
        NV_Ith_S(y, i) = ode.y0_d()[i];
      }
      for (size_t is = 0; is < ns; ++is) {
        N_VConst(RCONST(0.0), ys[is]);
      }

      try {
        CHECK_SUNDIALS_CALL(CVodeReInit(mem, ode.t0(), y));
        CHECK_SUNDIALS_CALL(CVodeSStolerances(mem, rtol_, atol_));
        CHECK_SUNDIALS_CALL(CVodeSetUserData(mem, ode.to_user_data()));
        CHECK_SUNDIALS_CALL(CVodeSetMaxNumSteps(mem, max_num_steps_));
        CHECK_SUNDIALS_CALL(CVodeSetMaxErrTestFails(mem, 20));
        CHECK_SUNDIALS_CALL(CVodeSetMaxConvFails(mem, 20));
#ifdef TORSTEN_CVS_JAC_AD
        CHECK_SUNDIALS_CALL(CVDlsSetJacFn(mem, cvodes_jac<Ode>()));
#endif

        /** if y0 is parameter, the first n sensitivity vector
         * are regarding y0, thus they form a unit matrix.
         **/
        if (Ode::need_fwd_sens) {
          if (Ode::is_var_y0) for (size_t i = 0; i < n; ++i) NV_Ith_S(ys[i], i) = 1.0;
          CHECK_SUNDIALS_CALL(CVodeSensInit(mem, ode.ns(), CV_STAGGERED, cvodes_sens_rhs<Ode>(), ys));  // NOLINT
          CHECK_SUNDIALS_CALL(CVodeSensEEtolerances(mem));
        }

        // the return type for MPI version is based on @c double,
        // as the return consists of the CVODES solution
        // instead of the assembled @c var vector in the
        // sequential version.
        solve(ode, observer.y);
      } catch (const std::exception& e) {
        CVodeSensFree(mem);
        throw;
      }

      CVodeSensFree(mem);

      return observer.y;
    }
  };  // cvodes integrator

  /**
   * Solve ODE system, no sensitivty, sequential
   *
   * @tparam F ODE functor type
   * @param[out] ode ODE system
   * @param[out] res_y ODE solutions
   */
  template <typename F, int Lmm, PMXCvodesSensMethod Sm>
  void PMXCvodesIntegrator::solve(PMXCvodesFwdSystem<F, double,
                                 double, double, Lmm, Sm>& ode,
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
   * Solve ODE system, no sensitivty, MPI
   *
   * @tparam F ODE functor type
   * @param[out] ode ODE system
   * @param[out] res_y ODE solutions
   */
  template <typename F, int Lmm, PMXCvodesSensMethod Sm>
  void PMXCvodesIntegrator::solve(PMXCvodesFwdSystem<F, double, double, double, Lmm, Sm>& ode,
                                 Eigen::MatrixXd& res_y) {
    double t1 = ode.t0();
    const std::vector<double>& ts = ode.ts();
    auto mem = ode.mem();
    auto y = ode.nv_y();
    const int n = ode.n();

    for (size_t i = 0; i < ts.size(); ++i) {
      CHECK_SUNDIALS_CALL(CVode(mem, ts[i], y, &t1, CV_NORMAL));
      res_y.block(0, i, n, 1) = Eigen::VectorXd::Map(NV_DATA_S(y), n);
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
  template <typename F, int Lmm, PMXCvodesSensMethod Sm>
  void PMXCvodesIntegrator::solve(PMXCvodesFwdSystem<F,
                                 stan::math::var, double, double, Lmm, Sm>& ode,
                                 std::vector<std::vector<stan::math::var> >& res_y) {
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
    std::vector<double> g(ts.size(), 0.0);

    for (size_t i = 0; i < ts.size(); ++i) {
      CHECK_SUNDIALS_CALL(CVode(mem, ts[i].val(), y, &t1, CV_NORMAL));
      ode.eval_rhs(t1, y);
      for (size_t j = 0; j < ode.n(); ++j) {
        g[i] = ode.fval()[j];
        // FIXME: use ts[i] instead of ts
        res_y[i][j] = precomputed_gradients(NV_Ith_S(y, j), ts, g);
        g[i] = 0.0;
      }
    }
  }

  /**
   * Solve ODE system, with only @c ts sensitivity. In this
   * case we don't have to solve a fwd sensitivty system. As
   * the sensitivity regarding @c ts[i] is just the RHS of ODE
   * for solution at @c ts[i], and is 0 for solution at
   * other time. Since @c fval stores the last RHS
   * evaluation in @c eval_f, we can just use it instead
   * recompute RHS. When MPI is used, we don't generate @c var
   * vector but simply return CVODES results consisting of
   * the ODE sensitivity solutions.
   *
   * @tparam F ODE functor type
   * @param[out] ode ODE system
   * @param[out] res_y ODE solutions, each @c vector
   * corresponds a time step, arranged as:
   * @c res_y = [y1, y2, ... dy1/dt1, dy2/dt1...dy1/dtn, dy2/dtn...]
   */
  template <typename F, int Lmm, PMXCvodesSensMethod Sm>
  void PMXCvodesIntegrator::solve(PMXCvodesFwdSystem<F, stan::math::var, double, double, Lmm, Sm>& ode, // NOLINT
                                 Eigen::MatrixXd& res_y) {
    double t1 = ode.t0();
    const std::vector<stan::math::var>& ts = ode.ts();
    auto mem = ode.mem();
    auto y = ode.nv_y();
    const int n = ode.n();

    for (size_t i = 0; i < ts.size(); ++i) {
      CHECK_SUNDIALS_CALL(CVode(mem, ts[i].val(), y, &t1, CV_NORMAL));
      res_y.block(0, i, n, 1) = Eigen::VectorXd::Map(NV_DATA_S(y), n);

      ode.eval_rhs(t1, y);
      res_y.block(n + i * n, i, n, 1) = Eigen::VectorXd::Map(ode.fval().data(), n);
    }
  }

  /**
   * Solve Ode system with forward sensitivty, return a
   * vector of @c var with precomputed gradient as sensitivity
   * value, when @c y0 is parameter but @theta & @c ts are not.
   *
   * @tparam Ode ODE system type
   * @param[out] ode ODE system
   * @param[out] res_y ODE solutions
   */
  template <typename F, typename Ty0, int Lmm, PMXCvodesSensMethod Sm>
  void PMXCvodesIntegrator::solve(PMXCvodesFwdSystem<F, double, Ty0, double, Lmm, Sm>& ode, // NOLINT
                                  std::vector<std::vector<stan::math::var> >& res_y) { // NOLINT
    using stan::math::precomputed_gradients;
    double t1 = ode.t0();
    const std::vector<double>& ts = ode.ts();
    auto mem = ode.mem();
    auto y = ode.nv_y();
    auto ys = ode.nv_ys();
    const auto n = ode.n();
    const auto ns = ode.ns();

    std::vector<double> g(n);

    for (size_t i = 0; i < ts.size(); ++i) {
      CHECK_SUNDIALS_CALL(CVode(mem, ts[i], y, &t1, CV_NORMAL));
      if (ode.need_fwd_sens) {
        CHECK_SUNDIALS_CALL(CVodeGetSens(mem, &t1, ys));
        for (size_t k = 0; k < n; ++k) {
          for (size_t j = 0; j < ns; ++j) {
            g[j] = NV_Ith_S(ys[j], k);
          }
          res_y[i][k] = precomputed_gradients(NV_Ith_S(y, k), ode.y0(), g);
        }
      }
    }
  }

  /**
   * Solve Ode system with forward sensitivty, return a
   * vector of @c var with precomputed gradient as sensitivity
   * value, when @c theta is parameter but @c y0 & @c ts are not.
   *
   * @tparam Ode ODE system type
   * @param[out] ode ODE system
   * @param[out] res_y ODE solutions
   */
  template <typename F, typename Tpar, int Lmm, PMXCvodesSensMethod Sm>
  void PMXCvodesIntegrator::solve(PMXCvodesFwdSystem<F, double, double, Tpar, Lmm, Sm>& ode, // NOLINT
                                  std::vector<std::vector<stan::math::var> >& res_y) { // NOLINT
    using stan::math::precomputed_gradients;
    double t1 = ode.t0();
    const std::vector<double>& ts = ode.ts();
    auto mem = ode.mem();
    auto y = ode.nv_y();
    auto ys = ode.nv_ys();
    const auto n = ode.n();
    const auto ns = ode.ns();

    std::vector<double> g(ns);

    for (size_t i = 0; i < ts.size(); ++i) {
      CHECK_SUNDIALS_CALL(CVode(mem, ts[i], y, &t1, CV_NORMAL));
      if (ode.need_fwd_sens) {
        CHECK_SUNDIALS_CALL(CVodeGetSens(mem, &t1, ys));
        for (size_t k = 0; k < n; ++k) {
          for (size_t j = 0; j < ns; ++j) {
            g[j] = NV_Ith_S(ys[j], k);
          }
          res_y[i][k] = precomputed_gradients(NV_Ith_S(y, k), ode.theta(), g);
        }
      }
    }
  }

  /**
   * Solve Ode system with forward sensitivty, return a
   * vector of @c var with precomputed gradient as sensitivity
   * value, when @c y0 and @c theta are parameters but @c ts is not.
   *
   * @tparam Ode ODE system type
   * @param[out] ode ODE system
   * @param[out] res_y ODE solutions
   */
  template <typename F, typename Ty0, typename Tpar, int Lmm, PMXCvodesSensMethod Sm> // NOLINT
  void PMXCvodesIntegrator::solve(PMXCvodesFwdSystem<F, double, Ty0, Tpar, Lmm, Sm>& ode, // NOLINT
                                 std::vector<std::vector<stan::math::var> >& res_y) { // NOLINT
    using stan::math::precomputed_gradients;
    double t1 = ode.t0();
    const std::vector<double>& ts = ode.ts();
    auto mem = ode.mem();
    auto y = ode.nv_y();
    auto ys = ode.nv_ys();
    const auto n = ode.n();
    const auto ns = ode.ns();
    auto vars = ode.vars();

    std::vector<double> g(ns);

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
   * Solve Ode system with forward sensitivty, 
   * and return a
   * vector of @c double of ODE sensitivity equation
   * solutions for MPI applications,
   * when @c y0 and/or @c theta are parameters but @c ts is not.
   *
   * @tparam Ode ODE system type
   * @param[out] ode ODE system
   * @param[out] res_y ODE solutions
   */
  template <typename F, typename Ty0, typename Tpar, int Lmm, PMXCvodesSensMethod Sm> // NOLINT
  void PMXCvodesIntegrator::solve(PMXCvodesFwdSystem<F, double, Ty0, Tpar, Lmm, Sm>& ode, // NOLINT
                                  Eigen::MatrixXd& res_y) {
    double t1 = ode.t0();
    const std::vector<double>& ts = ode.ts();
    auto mem = ode.mem();
    auto y = ode.nv_y();
    auto ys = ode.nv_ys();
    const auto n = ode.n();
    const auto ns = ode.ns();

    std::vector<double> g(ns);

    for (size_t i = 0; i < ts.size(); ++i) {
      CHECK_SUNDIALS_CALL(CVode(mem, ts[i], y, &t1, CV_NORMAL));
      res_y.block(0, i, n, 1) = Eigen::VectorXd::Map(NV_DATA_S(y), n);

      CHECK_SUNDIALS_CALL(CVodeGetSens(mem, &t1, ys));
      for (size_t j = 0; j < ns; ++j) {
        res_y.block(n + j * n, i, n, 1) = Eigen::VectorXd::Map(NV_DATA_S(ys[j]), n);
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
  template <typename F, typename Ty0, typename Tpar, int Lmm, PMXCvodesSensMethod Sm> // NOLINT
  void PMXCvodesIntegrator::solve(PMXCvodesFwdSystem<F, stan::math::var, Ty0, Tpar, Lmm, Sm>& ode, // NOLINT
                                 std::vector<std::vector<stan::math::var> >& res_y) { // NOLINT
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

    std::vector<double> g(ns + ts.size(), 0.0);

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

  /**
   * Solve Ode system with forward sensitivty, return a
   * vector of var with precomputed gradient as sensitivity
   * value, when @c ts, @c theta and/or @c y0 are parameters.
   *
   * @tparam Ode ODE system type
   * @param[out] ode ODE system
   * @param[out] res_y ODE solutions with sensitivity,
   * arranged as (sol value, grad(y0), grad(theta), grad(ts))
   */
  template <typename F, typename Ty0, typename Tpar, int Lmm, PMXCvodesSensMethod Sm> // NOLINT
  void PMXCvodesIntegrator::solve(PMXCvodesFwdSystem<F, stan::math::var, Ty0, Tpar, Lmm, Sm>& ode, // NOLINT
                                 Eigen::MatrixXd& res_y) {
    double t1 = ode.t0();
    const std::vector<stan::math::var>& ts = ode.ts();
    auto mem = ode.mem();
    auto y = ode.nv_y();
    auto ys = ode.nv_ys();
    const int n = ode.n();
    const int ns = ode.ns();

    for (size_t i = 0; i < ts.size(); ++i) {
      CHECK_SUNDIALS_CALL(CVode(mem, stan::math::value_of(ts[i]), y, &t1, CV_NORMAL));
      res_y.block(0, i, n, 1) = Eigen::VectorXd::Map(NV_DATA_S(y), n);

      CHECK_SUNDIALS_CALL(CVodeGetSens(mem, &t1, ys));
      for (size_t j = 0; j < ns; ++j) {
        res_y.block(n + j * n, i, n, 1) = Eigen::VectorXd::Map(NV_DATA_S(ys[j]), n);
      }

      ode.eval_rhs(t1, y);
      res_y.block(n + (i + ns) * n, i, n, 1) = Eigen::VectorXd::Map(ode.fval().data(), n);      
    }
  }

}  // namespace dsolve
}  // namespace torsten

#endif
