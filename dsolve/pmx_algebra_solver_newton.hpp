#ifndef STAN_MATH_TORSTEN_NEWTON_SOLVER_HPP
#define STAN_MATH_TORSTEN_NEWTON_SOLVER_HPP

#include <stan/math/rev/meta.hpp>
#include <stan/math/rev/core.hpp>
#include <stan/math/rev/fun/value_of.hpp>
#include <stan/math/rev/functor/algebra_system.hpp>
#include <stan/math/rev/functor/kinsol_data.hpp>
#include <stan/math/prim/err.hpp>
#include <stan/math/rev/functor/algebra_system.hpp>
#include <stan/math/prim/fun/mdivide_left.hpp>
#include <stan/math/prim/fun/to_array_1d.hpp>
#include <stan/math/prim/fun/to_vector.hpp>
#include <stan/math/prim/fun/value_of.hpp>
#include <stan/math/prim/functor/for_each.hpp>
#include <stan/math/torsten/dsolve/sundials_check.hpp>
#include <stan/math/torsten/dsolve/ode_tuple_functor.hpp>
#include <stan/math/torsten/mpi/precomputed_gradients.hpp>
#include <stan/math/torsten/meta/is_nl_system.hpp>
#include <kinsol/kinsol.h>
#include <sunmatrix/sunmatrix_dense.h>
#include <sunlinsol/sunlinsol_dense.h>
#include <sundials/sundials_context.h>
#include <nvector/nvector_serial.h>
#include <algorithm>
#include <iostream>
#include <string>
#include <vector>

namespace torsten {
  /** 
   * Nonlinear system
   * 
   */
  template <typename F, typename T_init, typename... T_par>
  struct PMXNLSystem : torsten::is_nl_system<F, T_par...> {
    using Nl = PMXNLSystem<F, T_init, T_par...>;
    using scalar_t = typename stan::return_type_t<T_par...>;
    using state_t = Eigen::Matrix<scalar_t, -1, 1>;
    static constexpr bool is_var_x0  = stan::is_var<T_init>::value;
    static constexpr bool is_var_par = stan::is_var<stan::return_type_t<T_par...>>::value;

    const F& f_;
    const dsolve::TupleNlFunc<F> f_tuple_;
    const Eigen::Matrix<T_init, -1, 1>& x0_;
    std::tuple<decltype(stan::math::deep_copy_vars(std::declval<const T_par&>()))...> theta_tuple_;
    // std::tuple<decltype(stan::math::value_of(std::declval<const T_par&>()))...> theta_dbl_tuple_;
    std::tuple<const T_par&...> theta_ref_tuple_;
    const size_t N;
    const size_t M;
    const size_t ns;
    std::ostream* msgs_;

  public:
    PMXNLSystem(const F& f,
                const Eigen::Matrix<T_init, -1, 1>& x0,
                std::ostream* msgs,
                const T_par&... args)
      : f_(f),
        f_tuple_(f_),
        x0_(x0),
        theta_tuple_(stan::math::deep_copy_vars(args)...),
        // theta_dbl_tuple_(stan::math::value_of(args)...),
        theta_ref_tuple_(args...),
        N(x0.size()),
        M(stan::math::count_vars(args...)),
        ns((is_var_x0 ? N : 0) + M),
        msgs_(msgs)
    {
      const char* caller = "PMX nonlinear system";
    }

    /** 
     * User-defined function passed to KINSOL
     * 
     * @param x independent variable
     * @param f nonlinear functor evaluation
     * @param user_data
     * 
     * @return error code
     */
    static int kinsol_nl_system(N_Vector x, N_Vector f, void* user_data) {
      Nl* nl = static_cast<Nl*>(user_data);
      int n = nl -> N;
      Eigen::VectorXd x_eigen(Eigen::Map<Eigen::VectorXd>(NV_DATA_S(x), n));
      Eigen::Map<Eigen::VectorXd>(N_VGetArrayPointer(f), n)
        = stan::math::value_of(nl->f_tuple_(x_eigen, nl->msgs_, nl -> theta_tuple_));
      return 0;
    }

    /** Implements the user-defined Jacobian calculation function passed to KINSOL. */
    static int kinsol_jacobian(N_Vector x, N_Vector f, SUNMatrix Jfx, void *user_data,
                               N_Vector tmp1, N_Vector tmp2) {
      Nl* nl = static_cast<Nl*>(user_data);
      int n = nl -> N;
      Eigen::VectorXd x_eigen(Eigen::Map<Eigen::VectorXd>(NV_DATA_S(x), n));

      stan::math::nested_rev_autodiff nested;
      Eigen::Matrix<stan::math::var, -1, 1> x_var(x_eigen);
      Eigen::Matrix<stan::math::var, -1, 1> f_var(nl->f_tuple_(x_var, nl -> msgs_, nl -> theta_tuple_));
      Eigen::Map<Eigen::VectorXd>(N_VGetArrayPointer(f), n) = stan::math::value_of(f_var);
      for (int i = 0; i < n; ++i) {
        nested.set_zero_all_adjoints();
        grad(f_var(i).vi_);
        for (int j = 0; j < n; ++j) {
          SM_ELEMENT_D(Jfx, i, j) = x_var.adj()[j];
        }
      }
      return 0;
    }

    /**
     * Calculate Jacobian Jxy(Jacobian of unknown x w.r.t. the param y)
     * given the solution. Specifically, for
     *
     * f(x, y) = 0
     *
     * we have (Jpq = Jacobian matrix dq/dq)
     *
     * Jfx * Jxy + Jfy = 0
     *
     * therefore Jxy can be solved from system
     *
     * - Jfx * Jxy = Jfy
     *
     * Jfx and Jfy are obtained through one AD evaluation of f.
     * dense Jacobian Jxy solved through QR decomposition.
     *
     * @param x current solution
     * @return a vector with x's val and Jxy as gradient.
     */
    Eigen::Matrix<stan::math::var, -1, 1> jac_xy(const Eigen::VectorXd& x) {
      stan::math::for_each([](auto&& arg) { stan::math::zero_adjoints(arg); }, theta_tuple_);

      Eigen::MatrixXd jfx(N, N);
      Eigen::MatrixXd jfy(N, M);

      {
        stan::math::nested_rev_autodiff nested;
        Eigen::Matrix<stan::math::var, -1, 1> x_var(x);
        Eigen::Matrix<stan::math::var, -1, 1> f_var(f_tuple_(x_var, msgs_, theta_tuple_));
        Eigen::VectorXd g(M);
        for (auto i = 0; i < N; ++i) {
          if (i > 0) {
            nested.set_zero_all_adjoints();            
          }      
          f_var(i).grad();
          for (auto k = 0; k < N; ++k) {
            jfx(i, k) = x_var(k).adj();
          }
          if (is_var_par) {
            g.fill(0);
            stan::math::apply([&](auto&&... args) {accumulate_adjoints(g.data(), args...);}, theta_tuple_);
            for (auto k = 0; k < M; ++k) {
              jfy(i, k) = g[k];
            }
          }
          stan::math::for_each([](auto&& arg) { stan::math::zero_adjoints(arg); }, theta_tuple_);
        }
      }
      Eigen::MatrixXd jxy = jfx.colPivHouseholderQr().solve(-jfy);

      using stan::math::ChainableStack;
      Eigen::Matrix<stan::math::var, -1, 1> x_sol(N);
      stan::math::vari** varis = torsten::varis_from_ode_pars(theta_ref_tuple_);

      for (size_t i = 0; i < N; i++) {
        double* g = ChainableStack::instance_->memalloc_.alloc_array<double>(M);
        for (size_t k = 0; k < M; k++) {
          *(g + k) = jxy(i, k);        
        }
        x_sol[i] = new stan::math::precomputed_gradients_vari(x[i], M, varis, g);
      }
      return x_sol;
    }
  };

/**
 * KINSOL algebraic system data holder that handles
 * construction & destruction of SUNDIALS data, as well as
 * auxiliary data that will be used for functor evaluation.
 *
 */
struct KinsolNewtonService {
  /** SUNDIALS context */
  sundials::Context sundials_context_;
  /** KINSOL memory block */
  void* mem_;
  /** NVECTOR for unknowns */
  N_Vector nv_x_;
  /** matrix */
  SUNMatrix J;
  /** LINEAR SOLVER */
  SUNLinearSolver LS;
  /** NVECTOR for scaling u */
  N_Vector nv_u_scal_;
  /** NVECTOR for scaling f */
  N_Vector nv_f_scal_;

  /** Constructor when y is data */
  template <typename T, typename T_u, typename T_f>
  KinsolNewtonService(const Eigen::Matrix<T, -1, 1>& x,
                      const std::vector<T_u>& u_scale,
                      const std::vector<T_f>& f_scale)
      : sundials_context_(),
        mem_(KINCreate(sundials_context_)),
        nv_x_(N_VNew_Serial(x.size(), sundials_context_)),
        J(SUNDenseMatrix(x.size(), x.size(), sundials_context_)),
        LS(SUNLinSol_Dense(nv_x_, J, sundials_context_)),
        nv_u_scal_(N_VNew_Serial(x.size(), sundials_context_)),
        nv_f_scal_(N_VNew_Serial(x.size(), sundials_context_))
  {
    for (int i = 0; i < x.size(); ++i) {
      NV_Ith_S(nv_x_, i) = stan::math::value_of(x(i));
      NV_Ith_S(nv_u_scal_, i) = stan::math::value_of(u_scale[i]);
      NV_Ith_S(nv_f_scal_, i) = stan::math::value_of(f_scale[i]);
    }
  }

  ~KinsolNewtonService() {
    SUNLinSolFree(LS);
    SUNMatDestroy(J);
    N_VDestroy_Serial(nv_x_);
    N_VDestroy_Serial(nv_u_scal_);
    N_VDestroy_Serial(nv_f_scal_);
    KINFree(&mem_);
  }
};

/**
 * Newton solver for zero of form
 *
 * F(x; theta) = 0
 *
 * with x as unknowns and theta parameters.
 *
 * @tparam newton_env_type solver environment setup that handles
 *                     workspace & auxiliary data encapsulation & RAII, namely
 *                     the work environment. Currently only support KINSOL's
 *                     dense matrix.
 * @tparam newton_jac_type functor type for calculating the
 *                     jacobian. Currently only support @c
 *                     FixedPointADJac that obtain .
 * @tparam F RHS functor for fixed point iteration.
 * @tparam fp_jac_type functor type for calculating the jacobian
 */

struct NewtonSolver {
  /**
   * Solve nonlinear zero using KINSOL's Newton solvers
   *
   * @param x initial point and final solution.
   * @param env KINSOL solution environment
   * @param f_tol Function tolerance
   * @param max_num_steps max nb. of iterations.
   */
  template<typename Nl_type, typename T_u, typename T_f>
  Eigen::Matrix<std::conditional_t<Nl_type::is_var_par, stan::math::var, double>, -1, 1>
  solve(Nl_type& nl,
        const std::vector<T_u>& u_scale,
        const std::vector<T_f>& f_scale,
        double step_tol, double f_tol, int max_num_steps) {
    KinsolNewtonService env(nl.x0_, u_scale, f_scale);
    int N = nl.N;
    void* mem = env.mem_;

    CHECK_KINSOL_CALL(KINInit(mem, &nl.kinsol_nl_system, env.nv_x_));
    CHECK_KINSOL_CALL(KINSetLinearSolver(mem, env.LS, env.J));
    CHECK_KINSOL_CALL(KINSetNumMaxIters(mem, max_num_steps));
    CHECK_KINSOL_CALL(KINSetScaledStepTol(mem, step_tol));
    CHECK_KINSOL_CALL(KINSetFuncNormTol(mem, f_tol));
    CHECK_KINSOL_CALL(KINSetJacFn(mem, &nl.kinsol_jacobian));
    CHECK_KINSOL_CALL(KINSetUserData(mem, static_cast<void*>(&nl)));

    CHECK_KINSOL_CALL(KINSol(mem, env.nv_x_,
                             KIN_LINESEARCH, env.nv_u_scal_, env.nv_f_scal_));

    Eigen::VectorXd x(N);
    for (int i = 0; i < N; ++i) {
      x(i) = NV_Ith_S(env.nv_x_, i);
    }
    return result(x, nl);
  }

  template<typename Nl_type,
           typename = stan::require_not_var_t<typename Nl_type::scalar_t>>
  Eigen::VectorXd result(const Eigen::VectorXd& x, Nl_type& nl) {
    return x;
  }

  template<typename Nl_type,
           typename = stan::require_var_t<typename Nl_type::scalar_t>>
  Eigen::Matrix<stan::math::var, -1, 1> result(const Eigen::VectorXd& x, Nl_type& nl) {
    return nl.jac_xy(x);
  }
};

/**
 * Return zero <code>x</code> of the specified nonlinear system:
 *
 * F(x; theta) = 0
 *
 * given an initial guess of x, and parameters theta and data. Use the
 * KINSOL solver from the SUNDIALS suite.
 *
 * The user can also specify the scaling controls, the function
 * tolerance, and the maximum number of steps.
 *
 * @tparam F type of equation system function.
 * @tparam T type of initial guess vector. The final solution
 *           type doesn't depend on initial guess type,
 *           but we allow initial guess to be either data or param.
 * @tparam T_u type of scaling vector for unknowns. We allow
 *             it to be @c var because scaling could be parameter
 *             dependent. Internally these params are converted to data
 *             because scaling is applied.
 * @tparam T_f type of scaling vector for residual. We allow
 *             it to be @c var because scaling could be parameter
 *             dependent. Internally these params are converted to data
 *             because scaling is applied.
 *
 * @param[in] f Functor that evaluated the system of equations.
 * @param[in] x Vector of starting values.
 * @param[in] y Parameter vector for the equation system. The function
 *            is overloaded to treat y as a vector of doubles or of a
 *            a template type T.
 * @param[in] x_r Continuous data vector for the equation system.
 * @param[in] x_i Integer data vector for the equation system.
 * @param[in, out] msgs The print stream for warning messages.
 * @param[in] u_scale diagonal scaling matrix elements Du
 *                    such that Du*x has all components roughly the same
 *                    magnitude when x is close to a solution.
 *                    (ref. KINSOL user guide chap.2 sec. "Scaling")
 * @param[in] f_scale diagonal scaling matrix elements such
 *                    that Df*(x-f(x)) has all components roughly the same
 *                    magnitude when x is not too close to a solution.
 *                    (ref. KINSOL user guide chap.2 sec. "Scaling")
 * @param[in] f_tol Function-norm stopping tolerance.
 * @param[in] max_num_steps maximum number of function evaluations.
 * @throw <code>std::invalid_argument</code> if x has size zero.
 * @throw <code>std::invalid_argument</code> if x has non-finite elements.
 * @throw <code>std::invalid_argument</code> if y has non-finite elements.
 * @throw <code>std::invalid_argument</code> if dat has non-finite elements.
 * @throw <code>std::invalid_argument</code> if dat_int has non-finite elements.
 * @throw <code>std::invalid_argument</code> if scaled_step_size is strictly
 * negative.
 * @throw <code>std::invalid_argument</code> if function_tolerance is strictly
 * negative.
 * @throw <code>std::invalid_argument</code> if max_num_steps is not positive.
 * @throw <code>boost::math::evaluation_error</code> (which is a subclass of
 * <code>std::runtime_error</code>) if solver exceeds max_num_steps.
 */
template <typename F, typename T1, typename T_u, typename T_f, typename... Args>
Eigen::Matrix<typename PMXNLSystem<F, T1, Args...>::scalar_t, -1, 1> pmx_algebra_solver_newton_tol (
    const F& f, const Eigen::Matrix<T1, -1, 1>& x,
    const std::vector<T_u>& u_scale,
    const std::vector<T_f>& f_scale,
    double step_tol, double f_tol, int max_num_steps,
    std::ostream* msgs, const Args&... args) {
    PMXNLSystem<F, T1, Args...> nl(f, x, msgs, args...);
    NewtonSolver s;
    return s.solve(nl, u_scale, f_scale, step_tol, f_tol, max_num_steps);
}

//   const std::vector<double> scaling(x.size(), 1.0);
template <typename F, typename T1, typename T_u, typename T_f, typename... Args>
Eigen::Matrix<typename PMXNLSystem<F, T1, Args...>::scalar_t, -1, 1>
pmx_algebra_solver_newton (
    const F& f, const Eigen::Matrix<T1, -1, 1>& x,
    const std::vector<T_u>& u_scale,
    const std::vector<T_f>& f_scale,
    std::ostream* msgs, const Args&... args) {
    PMXNLSystem<F, T1, Args...> nl(f, x, msgs, args...);
    NewtonSolver s;
    Eigen::Matrix<typename PMXNLSystem<F, T1, Args...>::scalar_t, -1, 1> res = s.solve(nl, u_scale, f_scale, 1e-3, 1e-6, 100);
    return res;
}
}

#endif
