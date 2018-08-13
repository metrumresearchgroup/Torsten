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
class pk_cvodes_integrator {
  const double rtol_;
  const double atol_;
  const int64_t max_num_steps_;
  /**
   * Forward decl
   */
  template <typename Ode>
  void init_sensitivity(Ode& ode);

  /**
   * Placeholder for data-only cvodes_forword_system, no sensitivty
   *
   * @tparam F ODE functor type.
   * @param[in] ode ODE system
   */
  template <typename F, int Lmm>
  void init_sensitivity(pk_cvodes_fwd_system<F, double,
                        double, double, Lmm>& ode) {}

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
  void solve(pk_cvodes_fwd_system<F, double, double, double, Lmm>& ode,
             std::vector<std::vector<double> >& res_y);

  template <typename Ode>
  void solve(Ode& ode, typename Ode::return_type& res_y);

  // TODO(yizhang): adjoint sensitivity solver

 public:
  static constexpr int CVODES_MAX_STEPS = 500;

  /**
   * constructor
   * @param[in] rtol relative tolerance
   * @param[in] atol absolute tolerance
   * @param[in] max_num_steps max nb. of times steps
   */
  pk_cvodes_integrator(const double rtol, const double atol,
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

      if (Ode::need_sens) {
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
  void pk_cvodes_integrator::solve(pk_cvodes_fwd_system<F, double,
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
 * Solve Ode system with forward sensitivty, return a
 * vector of var with precomputed gradient as sensitivity value
 *
 * @tparam Ode ODE system type
 * @param[out] ode ODE system
 * @param[out] res_y ODE solutions
 */
template <typename Ode>
void pk_cvodes_integrator::solve(Ode& ode, typename Ode::return_type& res_y) {
  double t1 = ode.t0();
  const std::vector<double>& ts = ode.ts();
  size_t i = 0;
  auto mem = ode.mem();
  auto y = ode.nv_y();
  auto ys = ode.nv_ys();
  const auto n = ode.n();
  const auto ns = ode.ns();
  auto vars = ode.vars();

  std::vector<stan::math::var> sol_t(n);
  std::vector<double> sol_grad(ns);

  std::for_each(ts.begin(), ts.end(), [&](double t2) {
      CHECK_SUNDIALS_CALL(CVode(mem, t2, y, &t1, CV_NORMAL));
      CHECK_SUNDIALS_CALL(CVodeGetSens(mem, &t1, ys));
    for (size_t k = 0; k < n; ++k) {
      for (size_t j = 0; j < ns; ++j) {
        sol_grad[j] = NV_Ith_S(ys[j], k);
      }
      sol_t[k]
          = stan::math::precomputed_gradients(NV_Ith_S(y, k), vars, sol_grad);
    }
    res_y[i] = sol_t;
    ++i;
  });
}
}  // namespace dsolve
}  // namespace torsten

#endif
