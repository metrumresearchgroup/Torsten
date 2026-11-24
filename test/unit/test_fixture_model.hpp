#ifndef TEST_UNIT_TORSTEN_TEST_FIXTURE_MODEL
#define TEST_UNIT_TORSTEN_TEST_FIXTURE_MODEL

#include <gtest/gtest.h>
#include <stan/math/rev/fun/fmax.hpp>
#include <boost/numeric/odeint.hpp>
#include <stan/math/torsten/test/unit/test_macros.hpp>
#include <stan/math/torsten/test/unit/test_functors.hpp>
#include <stan/math/torsten/test/unit/test_util.hpp>
#include <test/unit/math/prim/functor/harmonic_oscillator.hpp>
#include <test/unit/math/prim/functor/lorenz.hpp>
#include <nvector/nvector_serial.h>
#include <iostream>
#include <sstream>
#include <vector>
#include <string>

template<typename T>
struct TorstenPMXTest;

/** 
 * CRTP-based root test fixture for NM-TRAN compatible inputs.
 * The child types derived from this root type will handle specific
 * models and events. 
 * 
 */
template<template<typename> class child_type, typename T>
struct TorstenPMXTest<child_type<T> > : public testing::Test {
  /// solver type for analytical solutions
  using sol1_t = std::tuple_element_t<0, T>;
  /// solver type for numerical ode solutions
  using sol2_t = std::tuple_element_t<1, T>;
  /// time type
  using time_t = std::tuple_element_t<2, T>;
  /// amount type
  using amt_t = std::tuple_element_t<3, T>;
  /// rate type
  using rate_t = std::tuple_element_t<4, T>;
  /// II type
  using ii_t = std::tuple_element_t<5, T>;
  /// param type
  using param_t = std::tuple_element_t<6, T>;
  /// F type
  using biovar_t = std::tuple_element_t<7, T>;
  /// lag time type
  using tlag_t = std::tuple_element_t<8, T>;
  // ODE functor type
  using ode_t = std::tuple_element_t<9, T>;

  stan::math::nested_rev_autodiff nested;
  int nt;
  std::vector<time_t> time;
  std::vector<amt_t> amt;
  std::vector<rate_t> rate;
  std::vector<int> cmt;
  std::vector<int> evid;
  std::vector<ii_t> ii;
  std::vector<int> addl;
  std::vector<int> ss;
  std::vector<std::vector<param_t> > theta;
  std::vector<Eigen::Matrix<param_t, -1, -1> > pMatrix;
  std::vector<std::vector<biovar_t> > biovar;
  std::vector<std::vector<tlag_t> > tlag;

  // solvers
  sol1_t sol1;
  sol2_t sol2;

  // for ODE integrator
  int ncmt;
  double t0;
  std::vector<double> x_r;
  std::vector<int> x_i;
  double rtol;
  double atol;
  int max_num_steps;
  double as_rtol;
  double as_atol;
  int as_max_num_steps;
  std::ostream* msgs;

  TorstenPMXTest() : 
    t0(0.0),
    rtol             {1.E-10},
    atol             {1.E-10},
    max_num_steps    {100000},
    as_rtol          {1.E-4},
    as_atol          {1.E-6},
    as_max_num_steps {100},
    msgs             {nullptr}
  {
    // nested.set_zero_all_adjoints();    
  }

  void reset_events(int n) {
    nt = n;
    time.resize(nt);
    amt .resize(nt);
    rate.resize(nt);
    cmt .resize(nt);
    evid.resize(nt);
    ii  .resize(nt);
    addl.resize(nt);
    ss  .resize(nt);

    std::fill(time.begin(), time.end(), 0);
    std::fill(amt .begin(), amt .end(), 0);
    std::fill(rate.begin(), rate.end(), 0);
    std::fill(cmt .begin(), cmt .end(), ncmt);
    std::fill(evid.begin(), evid.end(), 0);
    std::fill(ii  .begin(), ii  .end(), 0);
    std::fill(addl.begin(), addl.end(), 0);
    std::fill(ss  .begin(), ss  .end(), 0);
  }

  /** 
   * apply numerical ODE solver functor to child test fixture
   * 
   * @param sol solver
   * @param test_ptr pointer to child fixture
   * 
   * @return Events solutions
   */
  template<typename solver_func_t, typename child_test_t,
           stan::require_any_t<std::is_same<solver_func_t, pmx_solve_adams_functor>,
                               std::is_same<solver_func_t, pmx_solve_bdf_functor>,
                               std::is_same<solver_func_t, pmx_solve_rk45_functor> >* = nullptr>
  auto apply_solver(solver_func_t const& sol, child_test_t* test_ptr) {
    typename child_test_t::ode_t f;
    return sol(f, ncmt, time, amt, rate, ii, evid, cmt, addl, ss,
               theta, biovar, tlag,
               rtol, atol, max_num_steps,
               as_rtol, as_atol, as_max_num_steps,
               nullptr);
  }

  /** 
   * apply coupled numerical ODE solver functor to child test fixture
   * 
   * @param sol solver
   * @param test_ptr pointer to child fixture
   * 
   * @return Events solutions
   */
  template<typename solver_func_t, typename child_test_t,
           stan::require_any_t<std::is_same<solver_func_t, pmx_solve_onecpt_rk45_functor>,
                               std::is_same<solver_func_t, pmx_solve_onecpt_bdf_functor> >* = nullptr>
  auto apply_solver(solver_func_t const& sol, child_test_t* test_ptr) {
    typename child_test_t::ode_t f;
    return sol(f, ncmt - 2, time, amt, rate, ii, evid, cmt, addl, ss,
               theta, biovar, tlag,
               rtol, atol, max_num_steps,
               as_rtol, as_atol, as_max_num_steps,
               nullptr);
  }

  /** 
   * apply coupled numerical ODE solver functor to child test fixture
   * 
   * @param sol solver
   * @param test_ptr pointer to child fixture
   * 
   * @return Events solutions
   */
  template<typename solver_func_t, typename child_test_t,
           stan::require_any_t<std::is_same<solver_func_t, pmx_solve_twocpt_rk45_functor>,
                               std::is_same<solver_func_t, pmx_solve_twocpt_bdf_functor> >* = nullptr>
  auto apply_solver(solver_func_t const& sol, child_test_t* test_ptr) {
    typename child_test_t::ode_t f;
    return sol(f, ncmt - 3, time, amt, rate, ii, evid, cmt, addl, ss,
               theta, biovar, tlag,
               rtol, atol, max_num_steps,
               as_rtol, as_atol, as_max_num_steps,
               nullptr);
  }

  /** 
   * apply analytical solution functor to child test fixture
   * 
   * @param sol solver
   * @param test_ptr pointer to child fixture
   * 
   * @return Events solutions
   */
  template<typename solver_func_t, typename child_test_t,
           stan::require_any_t<std::is_same<solver_func_t, pmx_solve_onecpt_functor>,
                               std::is_same<solver_func_t, pmx_solve_twocpt_functor>,
                               std::is_same<solver_func_t, pmx_solve_onecpt_effcpt_functor>,
                               std::is_same<solver_func_t, pmx_solve_twocpt_effcpt_functor>>* = nullptr>
  auto apply_solver(solver_func_t const& sol, child_test_t* test_ptr) {
    return sol(time, amt, rate, ii, evid, cmt, addl, ss, theta, biovar, tlag);
  }

  /** 
   * apply linear ODE solver functor to child test fixture
   * 
   * @param sol solver
   * @param test_ptr pointer to child fixture
   * 
   * @return Events solutions
   */
  template<typename solver_func_t, typename child_test_t,
           stan::require_any_t<std::is_same<solver_func_t, pmx_solve_linode_functor>>* = nullptr>
  auto apply_solver(solver_func_t const& sol, child_test_t* test_ptr) {
    return sol(time, amt, rate, ii, evid,
               cmt, addl, ss, pMatrix,
               biovar, tlag);
  }

  void compare_val(Eigen::MatrixXd const& x) {
    child_type<T>& fixture = static_cast<child_type<T>&>(*this);
    auto res1 = apply_solver(sol1, &fixture);
    EXPECT_MAT_VAL_FLOAT_EQ(res1, x);
  }

  void compare_val(Eigen::MatrixXd const& x, double tol) {
    child_type<T>& fixture = static_cast<child_type<T>&>(*this);
    auto res1 = apply_solver(sol1, &fixture);
    EXPECT_MAT_VAL_NEAR(res1, x, tol);
  }

  void compare_rel_val(Eigen::MatrixXd const& x, double rtol) {
    child_type<T>& fixture = static_cast<child_type<T>&>(*this);
    Eigen::ArrayXXd res = (stan::math::value_of(apply_solver(sol1, &fixture))).array();
    Eigen::ArrayXXd xr = x.array();
    Eigen::ArrayXXd rel = (res + 1) / (xr + 1); // add 1 to remove zero entries
    Eigen::ArrayXXd one = Eigen::ArrayXXd::Zero(xr.rows(), xr.cols()) + 1;
    EXPECT_MAT_VAL_NEAR((rel.matrix()), (one.matrix()), rtol);
  }

  void compare_solvers_val() {
    child_type<T>& fixture = static_cast<child_type<T>&>(*this);
    auto res1 = apply_solver(sol1, &fixture);
    auto res2 = apply_solver(sol2, &fixture);
    EXPECT_MAT_VAL_FLOAT_EQ(res1, res2);
  }

  void compare_solvers_val(double tol) {
    child_type<T>& fixture = static_cast<child_type<T>&>(*this);
    auto res1 = apply_solver(sol1, &fixture);
    auto res2 = apply_solver(sol2, &fixture);
    EXPECT_MAT_VAL_NEAR(res1, res2, tol);
  }

  template<typename T1, typename T2, typename T3,
           stan::require_any_not_t<stan::is_var<T1>, stan::is_var<T2>, stan::is_var<T3>>* = nullptr>
  void compare_mat_adj(Eigen::Matrix<T1, -1, -1>& mat1,
                       Eigen::Matrix<T2, -1, -1>& mat2,
                       std::vector<T3>& p, double tol, const char* diagnostic_msg) {}

  template<typename T1, typename T2, typename T3,
           stan::require_all_t<stan::is_var<T1>, stan::is_var<T2>, stan::is_var<T3>>* = nullptr>
  void compare_mat_adj(Eigen::Matrix<T1, -1, -1>& mat1,
                       Eigen::Matrix<T2, -1, -1>& mat2,
                       std::vector<T3>& p, double tol, const char* diagnostic_msg) {
    EXPECT_MAT_ADJ_NEAR(mat1, mat2, p, nested, tol, diagnostic_msg);
  }

  void compare_solvers_adj(const std::vector<double>& p, double tol, const char*) {}

  void compare_solvers_adj(std::vector<stan::math::var>& p, double tol,
                           const char* diagnostic_msg) {
    child_type<T>& fixture = static_cast<child_type<T>&>(*this);
    auto res1 = apply_solver(sol1, &fixture);
    auto res2 = apply_solver(sol2, &fixture);
    EXPECT_MAT_ADJ_NEAR(res1, res2, p, this -> nested, tol, diagnostic_msg);
  }

#define ADD_TEST_FUNC(NAME)                             \
  template<typename x_type>                             \
  auto test_func_##NAME(std::vector<x_type> const& x) { \
    return test_func_##NAME##_impl(x, sol1);              \
  }

#define ADD_CPT_TEST_FUNC_IMPL(NAME, ...)                               \
  template<typename x_type, typename solver_func_t,                     \
           stan::require_any_t<std::is_same<solver_func_t, pmx_solve_onecpt_functor>, \
                               std::is_same<solver_func_t, pmx_solve_twocpt_functor>, \
                               std::is_same<solver_func_t, pmx_solve_onecpt_effcpt_functor>>* = nullptr> \
  auto test_func_##NAME##_impl(x_type const& x, solver_func_t const& s) { \
    return s(__VA_ARGS__);                                              \
  }

#define ADD_LIN_TEST_FUNC_IMPL(NAME, ...)                               \
  template<typename x_type, typename solver_func_t,                     \
           stan::require_any_t<std::is_same<solver_func_t, pmx_solve_linode_functor>>* = nullptr> \
  auto test_func_##NAME##_impl(x_type const& x, solver_func_t const& s) { \
    return s(__VA_ARGS__);                                              \
  }

#define ADD_ODE_TEST_FUNC_IMPL(NAME, ...)                               \
  template<typename x_type, typename solver_func_t,                     \
           stan::require_any_t<std::is_same<solver_func_t, pmx_solve_adams_functor>, \
                               std::is_same<solver_func_t, pmx_solve_bdf_functor>, \
                               std::is_same<solver_func_t, pmx_solve_rk45_functor> >* = nullptr> \
  auto test_func_##NAME##_impl(x_type const& x, solver_func_t const& s) { \
    ode_t f;                                                            \
    return s(f, ncmt, __VA_ARGS__, rtol, atol, max_num_steps, as_rtol, as_atol, as_max_num_steps, msgs); \
  }

  ADD_CPT_TEST_FUNC_IMPL(time, x, amt, rate, ii, evid, cmt, addl, ss, theta, biovar, tlag);
  ADD_ODE_TEST_FUNC_IMPL(time, x, amt, rate, ii, evid, cmt, addl, ss, theta, biovar, tlag);
  ADD_LIN_TEST_FUNC_IMPL(time, x, amt, rate, ii, evid, cmt, addl, ss, pMatrix, biovar, tlag);
  ADD_TEST_FUNC(time)

  ADD_CPT_TEST_FUNC_IMPL(amt, time, x, rate, ii, evid, cmt, addl, ss, theta, biovar, tlag);
  ADD_ODE_TEST_FUNC_IMPL(amt, time, x, rate, ii, evid, cmt, addl, ss, theta, biovar, tlag);
  ADD_LIN_TEST_FUNC_IMPL(amt, time, x, rate, ii, evid, cmt, addl, ss, pMatrix, biovar, tlag);
  ADD_TEST_FUNC(amt)

  ADD_CPT_TEST_FUNC_IMPL(rate, time, amt, x, ii, evid, cmt, addl, ss, theta, biovar, tlag);
  ADD_ODE_TEST_FUNC_IMPL(rate, time, amt, x, ii, evid, cmt, addl, ss, theta, biovar, tlag);
  ADD_LIN_TEST_FUNC_IMPL(rate, time, amt, x, ii, evid, cmt, addl, ss, pMatrix, biovar, tlag);
  ADD_TEST_FUNC(rate)

  ADD_CPT_TEST_FUNC_IMPL(ii, time, amt, rate, x, evid, cmt, addl, ss, theta, biovar, tlag);
  ADD_ODE_TEST_FUNC_IMPL(ii, time, amt, rate, x, evid, cmt, addl, ss, theta, biovar, tlag);
  ADD_LIN_TEST_FUNC_IMPL(ii, time, amt, rate, x, evid, cmt, addl, ss, pMatrix, biovar, tlag);
  ADD_TEST_FUNC(ii)

  ADD_CPT_TEST_FUNC_IMPL(theta, time, amt, rate, ii, evid, cmt, addl, ss, x, biovar, tlag);
  ADD_ODE_TEST_FUNC_IMPL(theta, time, amt, rate, ii, evid, cmt, addl, ss, x, biovar, tlag);
  template<typename x_type>
  auto test_func_theta(std::vector<x_type> const& x_) {
    std::vector<std::vector<x_type> > x{x_};
    return test_func_theta_impl(x, sol1);
  }

  ADD_CPT_TEST_FUNC_IMPL(biovar, time, amt, rate, ii, evid, cmt, addl, ss, theta, x, tlag);
  ADD_ODE_TEST_FUNC_IMPL(biovar, time, amt, rate, ii, evid, cmt, addl, ss, theta, x, tlag);
  ADD_LIN_TEST_FUNC_IMPL(biovar, time, amt, rate, ii, evid, cmt, addl, ss, pMatrix, x, tlag);
  template<typename x_type>
  auto test_func_biovar(std::vector<x_type> const& x_) {
    std::vector<std::vector<x_type> > x{x_};
    return test_func_biovar_impl(x, sol1);
  }

  ADD_CPT_TEST_FUNC_IMPL(tlag, time, amt, rate, ii, evid, cmt, addl, ss, theta, biovar, x);
  ADD_ODE_TEST_FUNC_IMPL(tlag, time, amt, rate, ii, evid, cmt, addl, ss, theta, biovar, x);
  ADD_LIN_TEST_FUNC_IMPL(tlag, time, amt, rate, ii, evid, cmt, addl, ss, pMatrix, biovar, x);
  template<typename x_type>
  auto test_func_tlag(std::vector<x_type> const& x_) {
    std::vector<std::vector<x_type> > x{x_};
    return test_func_tlag_impl(x, sol1);
  }

#undef ADD_CPT_TEST_FUNC_IMPL
#undef ADD_ODE_TEST_FUNC_IMPL
#undef ADD_TEST_FUNC

#define ADD_FD_TEST(NAME, ARG_VEC)                                      \
  void test_finite_diff_##NAME(double h, double tol) {                  \
    EXPECT_MAT_FUNC_POSITIVE_PARAM_NEAR_FD(test_func_##NAME, ARG_VEC, nested, h, tol, #NAME); \
  }

  ADD_FD_TEST(amt, amt);
  ADD_FD_TEST(time, time);
  ADD_FD_TEST(rate, rate);
  ADD_FD_TEST(ii, ii);
  ADD_FD_TEST(theta, theta[0]);
  ADD_FD_TEST(biovar, biovar[0]);
  ADD_FD_TEST(tlag, tlag[0]);

#undef ADD_FD_TEST

  // test overload signatures
#define ADD_OVERLOAD_TEST(TOL, ...)                                     \
  {                                                                     \
    auto x = s(time, amt, rate, ii, evid, cmt, addl, ss, __VA_ARGS__);  \
    EXPECT_MAT_VAL_FLOAT_EQ(x0, x);                                     \
    compare_mat_adj(x0, x, amt, TOL, "AMT");                            \
    compare_mat_adj(x0, x, rate, TOL, "RATE");                          \
    compare_mat_adj(x0, x, theta[0], TOL, "theta");                     \
  }

  template<typename solver_func_t,
           stan::require_any_t<std::is_same<solver_func_t, pmx_solve_onecpt_functor>,
                               std::is_same<solver_func_t, pmx_solve_twocpt_functor>,
                               std::is_same<solver_func_t, pmx_solve_onecpt_effcpt_functor>>* = nullptr>
  void compare_overload_impl(solver_func_t const& s, double tol) {
      auto x0 = s(time, amt, rate, ii, evid, cmt, addl, ss, theta,    biovar,    tlag);   
      ADD_OVERLOAD_TEST(tol, theta[0], biovar,    tlag);   
      ADD_OVERLOAD_TEST(tol, theta[0], biovar[0], tlag);   
      ADD_OVERLOAD_TEST(tol, theta[0], biovar[0], tlag[0]);
      ADD_OVERLOAD_TEST(tol, theta[0], biovar,    tlag[0]);
      ADD_OVERLOAD_TEST(tol, theta,    biovar[0], tlag);   
      ADD_OVERLOAD_TEST(tol, theta,    biovar[0], tlag[0]);
      ADD_OVERLOAD_TEST(tol, theta,    biovar,    tlag[0]);

      // optional args
      ADD_OVERLOAD_TEST(tol, theta[0], biovar);   
      ADD_OVERLOAD_TEST(tol, theta[0], biovar[0]);
      ADD_OVERLOAD_TEST(tol, theta,    biovar[0]);
      ADD_OVERLOAD_TEST(tol, theta,    biovar);
      ADD_OVERLOAD_TEST(tol, theta[0]);
      ADD_OVERLOAD_TEST(tol, theta);
  }

  template<typename solver_func_t,
           stan::require_any_t<std::is_same<solver_func_t, pmx_solve_onecpt_functor>,
                               std::is_same<solver_func_t, pmx_solve_twocpt_functor>,
                               std::is_same<solver_func_t, pmx_solve_onecpt_effcpt_functor>>* = nullptr>
  void compare_param_overload_impl(solver_func_t const& s, double tol) {
      auto x0 = s(time, amt, rate, ii, evid, cmt, addl, ss, theta,    biovar,    tlag);   
      ADD_OVERLOAD_TEST(tol, theta[0], biovar,    tlag);   
      ADD_OVERLOAD_TEST(tol, theta[0], biovar[0], tlag);   
      ADD_OVERLOAD_TEST(tol, theta[0], biovar[0], tlag[0]);
      ADD_OVERLOAD_TEST(tol, theta[0], biovar,    tlag[0]);
      ADD_OVERLOAD_TEST(tol, theta,    biovar[0], tlag);   
      ADD_OVERLOAD_TEST(tol, theta,    biovar[0], tlag[0]);
      ADD_OVERLOAD_TEST(tol, theta,    biovar,    tlag[0]);
  }
  
  template<typename solver_func_t,
           stan::require_any_t<std::is_same<solver_func_t, pmx_solve_linode_functor>>* = nullptr>
  void compare_param_overload_impl(solver_func_t const& s, double tol) {
      auto x0 = s(time, amt, rate, ii, evid, cmt, addl, ss, pMatrix,   biovar,    tlag);   
      ADD_OVERLOAD_TEST(tol, pMatrix[0], biovar,    tlag);   
      ADD_OVERLOAD_TEST(tol, pMatrix[0], biovar[0], tlag);   
      ADD_OVERLOAD_TEST(tol, pMatrix[0], biovar[0], tlag[0]);
      ADD_OVERLOAD_TEST(tol, pMatrix[0], biovar,    tlag[0]);
      ADD_OVERLOAD_TEST(tol, pMatrix,    biovar[0], tlag);   
      ADD_OVERLOAD_TEST(tol, pMatrix,    biovar[0], tlag[0]);
      ADD_OVERLOAD_TEST(tol, pMatrix,    biovar,    tlag[0]);
  }

#undef ADD_OVERLOAD_TEST

#define ADD_OVERLOAD_TEST(TOL, ...)                                     \
  {                                                                     \
    ode_t f;                                                            \
    auto x = s(f, ncmt, time, amt, rate, ii, evid, cmt, addl, ss, __VA_ARGS__); \
    EXPECT_MAT_VAL_NEAR(x0, x, TOL);                                    \
    compare_mat_adj(x0, x, amt, TOL, "AMT");                            \
    compare_mat_adj(x0, x, rate, TOL, "RATE");                          \
    compare_mat_adj(x0, x, theta[0], TOL, "theta");                     \
  }

  template<typename solver_func_t,
           stan::require_any_t<std::is_same<solver_func_t, pmx_solve_rk45_functor>,
                               std::is_same<solver_func_t, pmx_solve_adams_functor>,
                               std::is_same<solver_func_t, pmx_solve_bdf_functor>>* = nullptr>
  void compare_overload_impl(solver_func_t const& s, double tol) {
    ode_t f;
    auto x0 = s(f, ncmt, time, amt, rate, ii, evid, cmt, addl, ss, theta, biovar, tlag, nullptr);
    ADD_OVERLOAD_TEST(tol, theta[0], biovar,    tlag, nullptr);
    ADD_OVERLOAD_TEST(tol, theta[0], biovar[0], tlag, nullptr);   
    ADD_OVERLOAD_TEST(tol, theta[0], biovar[0], tlag[0], nullptr);
    ADD_OVERLOAD_TEST(tol, theta[0], biovar,    tlag[0], nullptr);
    ADD_OVERLOAD_TEST(tol, theta,    biovar[0], tlag, nullptr);   
    ADD_OVERLOAD_TEST(tol, theta,    biovar[0], tlag[0], nullptr);
    ADD_OVERLOAD_TEST(tol, theta,    biovar,    tlag[0], nullptr);

    // optional args
    ADD_OVERLOAD_TEST(tol, theta[0], biovar, nullptr);   
    ADD_OVERLOAD_TEST(tol, theta[0], biovar[0], nullptr);
    ADD_OVERLOAD_TEST(tol, theta,    biovar[0], nullptr);
    ADD_OVERLOAD_TEST(tol, theta,    biovar, nullptr);
    ADD_OVERLOAD_TEST(tol, theta[0], nullptr);
    ADD_OVERLOAD_TEST(tol, theta, nullptr);

    // ODE controls args
    ADD_OVERLOAD_TEST(tol, theta[0], biovar, 1.e-10, 1.e-10, 100000, 1.e-10, 1.e-10, 100, nullptr);
    ADD_OVERLOAD_TEST(tol, theta[0], biovar, 1.e-10, 1.e-10, 100000, nullptr);
    ADD_OVERLOAD_TEST(tol, theta[0], biovar[0], 1.e-10, 1.e-10, 100000, 1.e-10, 1.e-10, 100, nullptr);
    ADD_OVERLOAD_TEST(tol, theta,    biovar[0], 1.e-10, 1.e-10, 100000, nullptr);
    ADD_OVERLOAD_TEST(tol, theta,    biovar, 1.e-10, 1.e-10, 100000, 1.e-10, 1.e-10, 100, nullptr);
    ADD_OVERLOAD_TEST(tol, theta,    biovar, 1.e-10, 1.e-10, 100000, nullptr);
    ADD_OVERLOAD_TEST(tol, theta[0], 1.e-10, 1.e-10, 100000, 1.e-10, 1.e-10, 100, nullptr);
    ADD_OVERLOAD_TEST(tol, theta, 1.e-10, 1.e-10, 100000, nullptr);
  }

  template<typename solver_func_t,
           stan::require_any_t<std::is_same<solver_func_t, pmx_solve_rk45_functor>,
                               std::is_same<solver_func_t, pmx_solve_adams_functor>,
                               std::is_same<solver_func_t, pmx_solve_bdf_functor>>* = nullptr>
  void compare_param_overload_impl(solver_func_t const& s, double tol) {
    ode_t f;
    auto x0 = s(f, ncmt, time, amt, rate, ii, evid, cmt, addl, ss, theta, biovar, tlag, nullptr);

    // optional args
    ADD_OVERLOAD_TEST(tol, theta[0], biovar, nullptr);   
    ADD_OVERLOAD_TEST(tol, theta[0], biovar[0], nullptr);
    ADD_OVERLOAD_TEST(tol, theta,    biovar[0], nullptr);
    ADD_OVERLOAD_TEST(tol, theta,    biovar, nullptr);
    ADD_OVERLOAD_TEST(tol, theta[0], nullptr);
    ADD_OVERLOAD_TEST(tol, theta, nullptr);
  }

#undef ADD_OVERLOAD_TEST

  void compare_overload(double tol) {
    compare_overload_impl(sol1, tol);
  }

  void compare_param_overload(double tol) {
    compare_param_overload_impl(sol1, tol);
  }
};

#endif
