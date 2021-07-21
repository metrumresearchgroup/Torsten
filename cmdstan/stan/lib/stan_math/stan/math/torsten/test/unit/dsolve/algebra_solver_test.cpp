#include <stan/math/rev/core.hpp>
#include <stan/math/torsten/dsolve/pmx_algebra_solver_newton.hpp>
#include <test/unit/math/rev/fun/util.hpp>
#include <test/unit/math/rev/functor/util_algebra_solver.hpp>
#include <test/unit/util.hpp>
#include <gtest/gtest.h>
#include <iostream>
#include <fstream>
#include <vector>
#include <string>

using torsten::pmx_algebra_solver_newton;
using torsten::pmx_algebra_solver_newton_tol;

TEST(newton_solver, simple_eq) {
  const std::vector<double> scaling(2, 1.0);
  simple_eq_functor f;
  torsten::nl_system_adaptor<simple_eq_functor> f_nl(f);

  stan::math::nested_rev_autodiff nested;  
  Eigen::Matrix<stan::math::var, -1, 1> x_var(2);
  x_var << 1, 1;
  Eigen::Matrix<stan::math::var, -1, 1>  y(3);
  y << 5, 4, 2;
  std::vector<double> x_r;
  std::vector<int> x_i;
  Eigen::Matrix<stan::math::var, -1, -1> sol = pmx_algebra_solver_newton(f_nl, x_var, scaling, scaling, nullptr, y,x_r,x_i);

  EXPECT_FLOAT_EQ(sol(0).val(), (y(0) * y(1)).val());
  EXPECT_FLOAT_EQ(sol(1).val(), y(2).val());

  sol(0).grad();
  EXPECT_FLOAT_EQ(y(0).adj(), y(1).val());
  EXPECT_FLOAT_EQ(y(1).adj(), y(0).val());

  // soluton is independent from initial guess
  EXPECT_FLOAT_EQ(x_var(0).adj(), 0.0);
  EXPECT_FLOAT_EQ(x_var(1).adj(), 0.0);

  nested.set_zero_all_adjoints();  
  sol(1).grad();
  EXPECT_FLOAT_EQ(y(0).adj(), 0.0);
  EXPECT_FLOAT_EQ(y(1).adj(), 0.0);
  EXPECT_FLOAT_EQ(y(2).adj(), 1.0);

  // soluton is independent from initial guess
  EXPECT_FLOAT_EQ(x_var(0).adj(), 0.0);
  EXPECT_FLOAT_EQ(x_var(1).adj(), 0.0);
}

struct nonlinear_eq_func {
  template <typename T0, typename T1>
  inline Eigen::Matrix<stan::return_type_t<T0, T1>, Eigen::Dynamic, 1>
  operator()(const Eigen::Matrix<T0, Eigen::Dynamic, 1>& x,
             const Eigen::Matrix<T1, Eigen::Dynamic, 1>& y,
             const std::vector<double>& dat, const std::vector<int>& dat_int,
             std::ostream* pstream__) const {
    Eigen::Matrix<stan::return_type_t<T0, T1>, Eigen::Dynamic, 1> z(3);
    z(0) = x(2) * x(2) - y(2) * (x(0) - x(1) + 1);
    z(1) = x(0) * x(1) - y(0) * y(1);
    z(2) = x(2) / x(0) + y(2) / y(0);
    return z;
  }
};

struct nonlinear_eq_variadic_func {
  template <typename T0, typename T1>
  inline Eigen::Matrix<stan::return_type_t<T0, T1>, Eigen::Dynamic, 1>
  operator()(const Eigen::Matrix<T0, Eigen::Dynamic, 1>& x,
             std::ostream* pstream__,
             const T1& y0, const T1& y1, const T1& y2) const {
    Eigen::Matrix<stan::return_type_t<T0, T1>, Eigen::Dynamic, 1> z(3);
    z(0) = x(2) * x(2) - y2 * (x(0) - x(1) + 1);
    z(1) = x(0) * x(1) - y0 * y1;
    z(2) = x(2) / x(0) + y2 / y0;
    return z;
  }
};

TEST(newton_solver, nonlinear_eq) {
  const std::vector<double> scaling(3, 1.0);
  nonlinear_eq_func f;
  torsten::nl_system_adaptor<nonlinear_eq_func> f_nl(f);

  stan::math::nested_rev_autodiff nested;

  Eigen::Matrix<stan::math::var, -1, 1> x_var(3);
  x_var << -3, -3, -3;
  Eigen::Matrix<stan::math::var, -1, 1> y(3);
  y << 4, 6, 3;
  std::vector<double> x_r;
  std::vector<int> x_i;

  Eigen::Matrix<stan::math::var, -1, 1> sol = pmx_algebra_solver_newton(f_nl, x_var, scaling, scaling, nullptr, y, x_r, x_i);

  EXPECT_FLOAT_EQ(sol(0).val(), -4.0);
  EXPECT_FLOAT_EQ(sol(1).val(), -6.0);
  EXPECT_FLOAT_EQ(sol(2).val(), 3.0);

  nested.set_zero_all_adjoints();  
  sol(0).grad();
  EXPECT_FLOAT_EQ(y(0).adj(), -0.75);
  EXPECT_FLOAT_EQ(y(1).adj(), -0.25);
  EXPECT_FLOAT_EQ(y(2).adj(), 0.25);

  nested.set_zero_all_adjoints();  
  sol(1).grad();
  EXPECT_FLOAT_EQ(y(0).adj(), -0.375);
  EXPECT_FLOAT_EQ(y(1).adj(), -0.625);
  EXPECT_FLOAT_EQ(y(2).adj(), -0.375);

  nested.set_zero_all_adjoints();  
  sol(2).grad();
  EXPECT_FLOAT_EQ(y(0).adj(), -0.1875);
  EXPECT_FLOAT_EQ(y(1).adj(),  0.1875);
  EXPECT_FLOAT_EQ(y(2).adj(),  0.8125);
}

TEST(newton_solver, nonlinear_eq_variadic) {
  const std::vector<double> scaling(3, 1.0);
  nonlinear_eq_variadic_func f_nl;

  stan::math::nested_rev_autodiff nested;

  Eigen::Matrix<stan::math::var, -1, 1> x_var(3);
  x_var << -3, -3, -3;
  stan::math::var y0 = 4, y1 = 6, y2 = 3;
  std::vector<double> x_r;
  std::vector<int> x_i;

  Eigen::Matrix<stan::math::var, -1, 1> sol =
  pmx_algebra_solver_newton(f_nl, x_var, scaling, scaling, nullptr, y0, y1, y2);

  EXPECT_FLOAT_EQ(sol(0).val(), -4.0);
  EXPECT_FLOAT_EQ(sol(1).val(), -6.0);
  EXPECT_FLOAT_EQ(sol(2).val(), 3.0);

  nested.set_zero_all_adjoints();  
  sol(0).grad();
  EXPECT_FLOAT_EQ(y0.adj(), -0.75);
  EXPECT_FLOAT_EQ(y1.adj(), -0.25);
  EXPECT_FLOAT_EQ(y2.adj(), 0.25);

  nested.set_zero_all_adjoints();  
  sol(1).grad();
  EXPECT_FLOAT_EQ(y0.adj(), -0.375);
  EXPECT_FLOAT_EQ(y1.adj(), -0.625);
  EXPECT_FLOAT_EQ(y2.adj(), -0.375);

  nested.set_zero_all_adjoints();  
  sol(2).grad();
  EXPECT_FLOAT_EQ(y0.adj(), -0.1875);
  EXPECT_FLOAT_EQ(y1.adj(),  0.1875);
  EXPECT_FLOAT_EQ(y2.adj(),  0.8125);
}

TEST(newton_solver, nonlinear_eq_dbl) {
  const std::vector<double> scaling(3, 1.0);
  nonlinear_eq_func f;
  torsten::nl_system_adaptor<nonlinear_eq_func> f_nl(f);

  stan::math::nested_rev_autodiff nested;

  Eigen::Matrix<double, -1, 1> x(3);
  x << -3, -3, -3;
  Eigen::Matrix<stan::math::var, -1, 1> x_var(stan::math::to_var(x));
  Eigen::Matrix<double, -1, 1> y(3);
  y << 4, 6, 3;
  std::vector<double> x_r;
  std::vector<int> x_i;

  Eigen::Matrix<double, -1, 1> sol1 = pmx_algebra_solver_newton(f_nl, x_var, scaling, scaling, nullptr, y, x_r, x_i);
  Eigen::Matrix<double, -1, 1> sol2 = pmx_algebra_solver_newton(f_nl, x, scaling, scaling, nullptr, y, x_r, x_i);

  EXPECT_FLOAT_EQ(sol1(0), -4.0);
  EXPECT_FLOAT_EQ(sol1(1), -6.0);
  EXPECT_FLOAT_EQ(sol1(2), 3.0);

  EXPECT_FLOAT_EQ(sol2(0), -4.0);
  EXPECT_FLOAT_EQ(sol2(1), -6.0);
  EXPECT_FLOAT_EQ(sol2(2), 3.0);
}

TEST(newton_solver, no_solution) {
  const std::vector<double> scaling(3, 1.0);
  unsolvable_eq_functor f;
  torsten::nl_system_adaptor<unsolvable_eq_functor> f_nl(f);

  stan::math::nested_rev_autodiff nested;

  Eigen::Matrix<double, -1, 1> x(2);
  x << 1, 1;
  Eigen::Matrix<stan::math::var, -1, 1> x_var(stan::math::to_var(x));
  Eigen::Matrix<double, -1, 1> y(2);
  y << 1, 1;
  std::vector<double> x_r;
  std::vector<int> x_i;

  EXPECT_THROW_MSG(pmx_algebra_solver_newton(f_nl, x_var, scaling, scaling, nullptr, y, x_r, x_i),
                   std::domain_error, "KIN_LSETUP_FAIL");
}

TEST(newton_solver, error) {
  const std::vector<double> scaling(3, 1.0);
  nonlinear_eq_variadic_func f_nl;

  stan::math::nested_rev_autodiff nested;

  Eigen::Matrix<stan::math::var, -1, 1> x_var(3);
  x_var << -3, -3, -3;
  stan::math::var y0 = 4, y1 = 6, y2 = 3;
  std::vector<double> x_r;
  std::vector<int> x_i;

  Eigen::Matrix<stan::math::var, -1, 1> sol =
    pmx_algebra_solver_newton_tol(f_nl, x_var, scaling, scaling,
                                  1.e-4, 1.e-7, 50, nullptr, y0, y1, y2);

  EXPECT_FLOAT_EQ(sol(0).val(), -4.0);
  EXPECT_FLOAT_EQ(sol(1).val(), -6.0);
  EXPECT_FLOAT_EQ(sol(2).val(), 3.0);

  nested.set_zero_all_adjoints();  
  sol(0).grad();
  EXPECT_FLOAT_EQ(y0.adj(), -0.75);
  EXPECT_FLOAT_EQ(y1.adj(), -0.25);
  EXPECT_FLOAT_EQ(y2.adj(), 0.25);

  nested.set_zero_all_adjoints();  
  sol(1).grad();
  EXPECT_FLOAT_EQ(y0.adj(), -0.375);
  EXPECT_FLOAT_EQ(y1.adj(), -0.625);
  EXPECT_FLOAT_EQ(y2.adj(), -0.375);

  nested.set_zero_all_adjoints();  
  sol(2).grad();
  EXPECT_FLOAT_EQ(y0.adj(), -0.1875);
  EXPECT_FLOAT_EQ(y1.adj(),  0.1875);
  EXPECT_FLOAT_EQ(y2.adj(),  0.8125);

  EXPECT_THROW_MSG(
    pmx_algebra_solver_newton_tol(f_nl, x_var, scaling, scaling,
                                  1.e-4, 1.e-7, 5, nullptr, y0, y1, y2),
    std::domain_error, "KIN_MAXITER_REACHED");
}
