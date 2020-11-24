#ifndef STAN_MATH_TORSTEN_DSOLVE_IDAS_INTEGRATOR_HPP
#define STAN_MATH_TORSTEN_DSOLVE_IDAS_INTEGRATOR_HPP

#include <stan/math/prim/err/check_less.hpp>
#include <stan/math/prim/err/check_finite.hpp>
#include <stan/math/prim/err/check_nonzero_size.hpp>
#include <stan/math/prim/err/check_ordered.hpp>
#include <stan/math/rev/meta/is_var.hpp>
#include <stan/math/torsten/dsolve/pmx_idas_fwd_system.hpp>
#include <stan/math/torsten/dsolve/sundials_check.hpp>
#include <idas/idas.h>
#include <idas/idas_direct.h>
#include <sunmatrix/sunmatrix_dense.h>
#include <sunlinsol/sunlinsol_dense.h>
#include <nvector/nvector_serial.h>
#include <ostream>
#include <vector>
#include <algorithm>

namespace torsten {
namespace dsolve {
/**
 * IDAS DAE integrator.
 */
class PMXIdasIntegrator {
  const double rtol_;
  const double atol_;
  const int64_t max_num_steps_;
  /**
   * Forward decl
   */
  template <typename Dae>
  void init_sensitivity(Dae& dae);

  /**
   * Placeholder for data-only idas_forword_system, no sensitivty
   *
   * @tparam F DAE functor type.
   * @param[in] dae DAE system
   */
  template <typename F>
  void init_sensitivity(PMXIdasFwdSystem<F, double, double, double>& dae) {}

  // /**
  //  *  idas adjoint sens calculation requires different initialization
  //  *
  //  * @tparam F type of DAE RHS functor
  //  * @tparam Tyy type of DAE primary unknowns
  //  * @tparam Typ type of DAE derivative unknowns
  //  * @tparam Tpar type of DAE parameters.
  //  * @param[out] dae DAE system
  //  * @param[in] t0 initial time.
  //  * @param[in] ts times of the desired solutions
  //  * @param[out] res_yy DAE solutions
  //  */
  // template <typename F, typename Tyy, typename Typ, typename Tpar>
  // void init_sensitivity(idas_adjoint_system<F, Tyy, Typ, Tpar>& dae) {
  //   // TODO(yizhang): adjoint sensitivity initialization
  // }

  template <typename F>
  void solve(PMXIdasFwdSystem<F, double, double, double>& dae,
             const double& t0, const std::vector<double>& ts,
             std::vector<std::vector<double> >& res_yy);

  template <typename Dae>
  void solve(Dae& dae, const double& t0, const std::vector<double>& ts,
             typename Dae::return_type& res_yy);

  // TODO(yizhang): adjoint sensitivity solver

 public:
  static constexpr int IDAS_MAX_STEPS = 500;
  /**
   * constructor
   * @param[in] rtol relative tolerance
   * @param[in] atol absolute tolerance
   * @param[in] max_num_steps max nb. of times steps
   */
  PMXIdasIntegrator(const double rtol, const double atol,
                  const int64_t max_num_steps = IDAS_MAX_STEPS)
      : rtol_(rtol), atol_(atol), max_num_steps_(max_num_steps) {
    using stan::math::invalid_argument;
    if (rtol_ <= 0)
      invalid_argument("PMXIdasIntegrator", "relative tolerance,", rtol_, "",
                       ", must be greater than 0");
    if (rtol_ > 1.0E-3)
      invalid_argument("PMXIdasIntegrator", "relative tolerance,", rtol_, "",
                       ", must be less than 1.0E-3");
    if (atol_ <= 0)
      invalid_argument("PMXIdasIntegrator", "absolute tolerance,", atol_, "",
                       ", must be greater than 0");
    if (max_num_steps_ <= 0)
      invalid_argument("PMXIdasIntegrator", "max_num_steps,", max_num_steps_, "",
                       ", must be greater than 0");
  }

  /**
   * Return the solutions for the specified DAE
   * given the specified initial state,
   * initial times, times of desired solution, and parameters and
   * data, writing error and warning messages to the specified
   * stream contained in the DAE system.
   *
   * @tparam DAE type of DAE system
   * @param[in] dae DAE system
   * @param[in] t0 initial time.
   * @param[in] ts times of the desired solutions, in strictly
   * increasing order, all greater than the initial time.
   * @return a vector of states, each state being a vector of the
   * same size as the state variable, corresponding to a time in ts.
   */
  template <typename Dae>
  typename Dae::return_type integrate(Dae& dae, double t0,
                                      const std::vector<double>& ts) {
    using Eigen::Dynamic;
    using Eigen::Matrix;
    using Eigen::MatrixXd;
    using Eigen::VectorXd;
    using stan::math::check_finite;
    using stan::math::check_ordered;
    using stan::math::check_nonzero_size;
    using stan::math::check_less;

    static const char* caller = "PMXIdasIntegrator";
    check_finite(caller, "initial time", t0);
    check_finite(caller, "times", ts);
    check_ordered(caller, "times", ts);
    check_nonzero_size(caller, "times", ts);
    check_less(caller, "initial time", t0, ts.front());

    auto mem = dae.mem();
    auto yy = dae.nv_yy();
    auto yp = dae.nv_yp();
    const size_t n = dae.n();

    typename Dae::return_type res_yy(
        ts.size(), std::vector<typename Dae::scalar_type>(n, 0));

    auto A = SUNDenseMatrix(n, n);
    auto LS = SUNDenseLinearSolver(yy, A);

    try {
      CHECK_SUNDIALS_CALL(IDASetUserData(mem, dae.to_user_data()));

      CHECK_SUNDIALS_CALL(IDAInit(mem, dae.residual(), t0, yy, yp));
      CHECK_SUNDIALS_CALL(IDADlsSetLinearSolver(mem, LS, A));
      CHECK_SUNDIALS_CALL(IDASStolerances(mem, rtol_, atol_));
      CHECK_SUNDIALS_CALL(IDASetMaxNumSteps(mem, max_num_steps_));

      init_sensitivity(dae);

      solve(dae, t0, ts, res_yy);
    } catch (const std::exception& e) {
      SUNLinSolFree(LS);
      SUNMatDestroy(A);
      throw;
    }

    SUNLinSolFree(LS);
    SUNMatDestroy(A);

    return res_yy;
  }
};  // idas integrator

/**
 * Initialize sensitivity calculation and set
 * tolerance. For sensitivity with respect to initial
 * conditions, set sensitivity to identity
 *
 * @tparam Dae DAE system type
 * @param[in/out] dae DAE system
 */
template <typename Dae>
void PMXIdasIntegrator::init_sensitivity(Dae& dae) {
  if (Dae::need_sens) {
    auto mem = dae.mem();
    auto yys = dae.nv_yys();
    auto yps = dae.nv_yps();
    auto n = dae.n();

    if (Dae::is_var_yy0) {
      for (size_t i = 0; i < n; ++i)
        NV_Ith_S(yys[i], i) = 1.0;
    }
    if (Dae::is_var_yp0) {
      for (size_t i = 0; i < n; ++i)
        NV_Ith_S(yps[i + n], i) = 1.0;
    }
    CHECK_SUNDIALS_CALL(IDASensInit(mem, dae.ns(), IDA_SIMULTANEOUS,
                                NULL, yys, yps));
    CHECK_SUNDIALS_CALL(IDASetSensParams(mem, dae.idas_params(), NULL, NULL));
    CHECK_SUNDIALS_CALL(IDASensEEtolerances(mem));
    CHECK_SUNDIALS_CALL(IDAGetSensConsistentIC(mem, yys, yps));
  }
}

/**
 * Solve DAE system, no sensitivty
 *
 * @tparam F DAE functor type
 * @param[out] dae DAE system
 * @param[in] t0 initial time
 * @param[in] ts times of the desired solutions
 * @param[out] res_yy DAE solutions
 */
template <typename F>
void PMXIdasIntegrator::solve(PMXIdasFwdSystem<F, double, double, double>& dae,
                            const double& t0, const std::vector<double>& ts,
                            std::vector<std::vector<double> >& res_yy) {
  double t1 = t0;
  auto mem = dae.mem();
  auto yy = dae.nv_yy();
  auto yp = dae.nv_yp();

  for (size_t i = 0; i < ts.size(); ++i) {
    CHECK_SUNDIALS_CALL(IDASolve(mem, ts[i], &t1, yy, yp, IDA_NORMAL));
    for (size_t j = 0; j < dae.n(); ++j) res_yy[i][j] = NV_Ith_S(yy, j);    
  };
}

/**
 * Solve Dae system with forward sensitivty, return a
 * vector of var with precomputed gradient as sensitivity value
 *
 * @tparam Dae DAE system type
 * @param[out] dae DAE system
 * @param[in] t0 initial time
 * @param[in] ts times of the desired solutions
 * @param[out] res_yy DAE solutions
 */
template <typename Dae>
void PMXIdasIntegrator::solve(Dae& dae, const double& t0,
                            const std::vector<double>& ts,
                            typename Dae::return_type& res_yy) {
  double t1 = t0;
  size_t i = 0;
  auto mem = dae.mem();
  auto yy = dae.nv_yy();
  auto yp = dae.nv_yp();
  auto yys = dae.nv_yys();
  const auto n = dae.n();
  const auto ns = dae.ns();
  auto vars = dae.vars();

  std::vector<stan::math::var> sol_t(n);
  std::vector<double> sol_grad(ns);

  std::for_each(ts.begin(), ts.end(), [&](double t2) {
    CHECK_SUNDIALS_CALL(IDASolve(mem, t2, &t1, yy, yp, IDA_NORMAL));
    CHECK_SUNDIALS_CALL(IDAGetSens(mem, &t1, yys));
    for (size_t k = 0; k < n; ++k) {
      for (size_t j = 0; j < ns; ++j) {
        sol_grad[j] = NV_Ith_S(yys[j], k);
      }
      sol_t[k]
          = stan::math::precomputed_gradients(NV_Ith_S(yy, k), vars, sol_grad);
    }
    res_yy[i] = sol_t;
    ++i;
  });
}
}  // namespace dsolve
}  // namespace torsten

#endif
