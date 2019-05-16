#include <stan/math.hpp>
#include <stan/math/rev/core.hpp>
#include <stan/math/rev/mat/fun/typedefs.hpp>
#include <test/unit/math/rev/mat/fun/util.hpp>
#include <stan/math/torsten/dsolve/pmx_integrate_ode_rk45.hpp>
#include <test/unit/math/torsten/pmx_ode_test_fixture.hpp>
#include <test/unit/math/torsten/test_util.hpp>
#include <stan/math/torsten/mpi/envionment.hpp>
#include <test/unit/util.hpp>
#include <gtest/gtest.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <memory>
#include <iostream>
#include <fstream>
#include <limits>
#include <vector>
#include <chrono>
#include <ctime>

using stan::math::integrate_ode_rk45;
using torsten::pmx_integrate_ode_group_rk45;
using stan::math::matrix_v;
using stan::math::var;
using std::vector;

TEST_F(TorstenOdeTest_chem, group_rk45_ivp) {
  torsten::mpi::Envionment::init();

  const int np = 7;

  vector<int> len(np, ts.size());
  vector<double> ts_m;
  ts_m.reserve(np * ts.size());
  for (int i = 0; i < np; ++i) ts_m.insert(ts_m.end(), ts.begin(), ts.end());

  vector<vector<double> > y0_m (np, y0);
  vector<vector<double> > theta_m (np, theta);
  vector<vector<double> > x_r_m (np, x_r);
  vector<vector<int> > x_i_m (np, x_i);

  vector<vector<double> > y = integrate_ode_rk45(f, y0, t0, ts, theta , x_r, x_i); // NOLINT
  Eigen::MatrixXd y_m = pmx_integrate_ode_group_rk45(f, y0_m, t0, len, ts_m, theta_m , x_r_m, x_i_m);

  EXPECT_EQ(y_m.cols(), ts_m.size());
  int icol = 0;
  for (size_t i = 0; i < np; ++i) {
    Eigen::MatrixXd y_i = y_m.block(0, icol, y0.size(), len[i]);
    torsten::test::test_val(y, y_i);
    icol += len[i];
  }
}

TEST_F(TorstenOdeTest_chem, group_rk45_fwd_sensitivity_theta) {
  const int np = 2;

  std::vector<var> theta_var1 = stan::math::to_var(theta);
  std::vector<var> theta_var2 = stan::math::to_var(theta);
  theta_var2[0] *= 1.2;

  vector<vector<var> > theta_var_m(2);
  theta_var_m[0] = stan::math::to_var(theta);
  theta_var_m[1] = stan::math::to_var(theta);
  theta_var_m[1][0] *= 1.2;

  vector<int> len(np, ts.size());
  vector<double> ts_m;
  ts_m.reserve(np * ts.size());
  for (int i = 0; i < np; ++i) ts_m.insert(ts_m.end(), ts.begin(), ts.end());

  vector<vector<double> > y0_m {y0, y0};

  vector<vector<double> > x_r_m {x_r, x_r};
  vector<vector<int> > x_i_m {x_i, x_i};

  vector<vector<var> > y1 = stan::math::integrate_ode_rk45(f, y0, t0, ts, theta_var1, x_r, x_i);
  vector<vector<var> > y2 = stan::math::integrate_ode_rk45(f, y0, t0, ts, theta_var2, x_r, x_i);
  vector<vector<var> > y3 = pmx_integrate_ode_rk45(f, y0, t0, ts, theta_var1, x_r, x_i);
  vector<vector<var> > y4 = pmx_integrate_ode_rk45(f, y0, t0, ts, theta_var2, x_r, x_i);
  Eigen::Matrix<var, -1, -1> y_m = pmx_integrate_ode_group_rk45(f, y0_m, t0, len, ts_m, theta_var_m , x_r_m, x_i_m);

  EXPECT_EQ(y_m.cols(), ts_m.size());

  matrix_v y;

  int icol = 0;
  y = y_m.block(0, icol, y0.size(), len[0]);
  torsten::test::test_grad(theta_var1, theta_var_m[0], y1, y, 1e-7, 1e-6);
  icol += len[0];
  y = y_m.block(0, icol, y0.size(), len[1]);
  torsten::test::test_grad(theta_var2, theta_var_m[1], y2, y, 1e-7, 1e-6);

  icol = 0;
  y = y_m.block(0, icol, y0.size(), len[0]);
  torsten::test::test_grad(theta_var1, theta_var_m[0], y3, y, 1e-15, 1e-15);
  icol += len[0];
  y = y_m.block(0, icol, y0.size(), len[1]);
  torsten::test::test_grad(theta_var2, theta_var_m[1], y4, y, 1e-15, 1e-15);
}

TEST_F(TorstenOdeTest_neutropenia, group_pmx_rk45_fwd_sensitivity_theta) {
  // size of population
  const int np = 7;

  vector<var> theta_var = stan::math::to_var(theta);

  vector<int> len(np, ts.size());
  vector<double> ts_m;
  ts_m.reserve(np * ts.size());
  for (int i = 0; i < np; ++i) ts_m.insert(ts_m.end(), ts.begin(), ts.end());

  vector<vector<double> > y0_m (np, y0);
  vector<vector<var> > theta_var_m (np, stan::math::to_var(theta));
  vector<vector<double> > x_r_m (np, x_r);
  vector<vector<int> > x_i_m (np, x_i);

  vector<vector<var> > y = pmx_integrate_ode_rk45(f, y0, t0, ts, theta_var, x_r, x_i);
  Eigen::Matrix<var, -1, -1> y_m = pmx_integrate_ode_group_rk45(f, y0_m, t0, len, ts_m, theta_var_m , x_r_m, x_i_m);

  EXPECT_EQ(y_m.cols(), ts_m.size());
  int icol = 0;
  for (int i = 0; i < np; ++i) {
    matrix_v y_i = y_m.block(0, icol, y0.size(), len[i]);
    torsten::test::test_grad(theta_var, theta_var_m[i], y, y_i, 1e-16, 1e-16);
    icol += len[i];
  }
}

TEST_F(TorstenOdeTest_neutropenia, group_rk45_fwd_sensitivity_theta) {
  // size of population
  const int np = 7;

  vector<var> theta_var = stan::math::to_var(theta);

  vector<int> len(np, ts.size());
  vector<double> ts_m;
  ts_m.reserve(np * ts.size());
  for (int i = 0; i < np; ++i) ts_m.insert(ts_m.end(), ts.begin(), ts.end());

  vector<vector<double> > y0_m (np, y0);
  vector<vector<var> > theta_var_m (np, stan::math::to_var(theta));
  vector<vector<double> > x_r_m (np, x_r);
  vector<vector<int> > x_i_m (np, x_i);

  vector<vector<var> > y = integrate_ode_rk45(f, y0, t0, ts, theta_var, x_r, x_i);
  Eigen::Matrix<var, -1, -1> y_m = pmx_integrate_ode_group_rk45(f, y0_m, t0, len, ts_m, theta_var_m , x_r_m, x_i_m);

  EXPECT_EQ(y_m.cols(), ts_m.size());
  int icol = 0;
  for (int i = 0; i < np; ++i) {
    matrix_v y_i = y_m.block(0, icol, y0.size(), len[i]);
    torsten::test::test_grad(theta_var, theta_var_m[i], y, y_i, 1e-8, 1e-8);
    icol += len[i];
  }
}

TEST_F(TorstenOdeTest_sho, group_rk45_fwd_sensitivity_theta) {
  const int np = 5;

  vector<var> theta_var = stan::math::to_var(theta);

  vector<int> len(np, ts.size());
  vector<double> ts_m;
  ts_m.reserve(np * ts.size());
  for (int i = 0; i < np; ++i) ts_m.insert(ts_m.end(), ts.begin(), ts.end());

  vector<vector<double> > y0_m (np, y0);
  vector<vector<var> > theta_var_m (np, stan::math::to_var(theta));
  vector<vector<double> > x_r_m (np, x_r);
  vector<vector<int> > x_i_m (np, x_i);

  vector<vector<var> > y = integrate_ode_rk45(f, y0, t0, ts, theta_var, x_r, x_i);
  Eigen::Matrix<var, -1, -1> y_m = pmx_integrate_ode_group_rk45(f, y0_m, t0, len, ts_m, theta_var_m , x_r_m, x_i_m);

  EXPECT_EQ(y_m.cols(), ts_m.size());
  int icol = 0;
  for (int i = 0; i < np; ++i) {
    matrix_v y_i = y_m.block(0, icol, y0.size(), len[i]);
    torsten::test::test_grad(theta_var, theta_var_m[i], y, y_i, 1e-8, 1e-6);
    icol += len[i];
  }
}

TEST_F(TorstenOdeTest_chem, group_pmx_rk45_fwd_sensitivity_theta) {
  const int np = 5;

  vector<var> theta_var = stan::math::to_var(theta);

  vector<int> len(np, ts.size());
  vector<double> ts_m;
  ts_m.reserve(np * ts.size());
  for (int i = 0; i < np; ++i) ts_m.insert(ts_m.end(), ts.begin(), ts.end());

  vector<vector<double> > y0_m (np, y0);
  vector<vector<var> > theta_var_m (np, stan::math::to_var(theta));
  vector<vector<double> > x_r_m (np, x_r);
  vector<vector<int> > x_i_m (np, x_i);

  vector<vector<var> > y = pmx_integrate_ode_rk45(f, y0, t0, ts, theta_var, x_r, x_i);
  Eigen::Matrix<var, -1, -1> y_m = pmx_integrate_ode_group_rk45(f, y0_m, t0, len, ts_m, theta_var_m , x_r_m, x_i_m);

  EXPECT_EQ(y_m.cols(), ts_m.size());
  int icol = 0;
  for (int i = 0; i < np; ++i) {
    matrix_v y_i = y_m.block(0, icol, y0.size(), len[i]);
    torsten::test::test_grad(theta_var, theta_var_m[i], y, y_i, 1e-16, 1e-16);
    icol += len[i];
  }
}

TEST_F(TorstenOdeTest_sho, group_pmx_rk45_fwd_sensitivity_y0) {
  const int np = 10;

  vector<var> y0_var = stan::math::to_var(y0);

  vector<int> len(np, ts.size());
  vector<double> ts_m;
  ts_m.reserve(np * ts.size());
  for (int i = 0; i < np; ++i) ts_m.insert(ts_m.end(), ts.begin(), ts.end());

  vector<vector<double> > theta_m (np, theta);
  vector<vector<var> > y0_var_m (np, stan::math::to_var(y0));
  vector<vector<double> > x_r_m (np, x_r);
  vector<vector<int> > x_i_m (np, x_i);

  vector<vector<var> > y = pmx_integrate_ode_rk45(f, y0_var, t0, ts, theta, x_r, x_i);
  Eigen::Matrix<var, -1, -1> y_m = pmx_integrate_ode_group_rk45(f, y0_var_m, t0, len, ts_m, theta_m , x_r_m, x_i_m);

  EXPECT_EQ(y_m.cols(), ts_m.size());
  int icol = 0;
  for (int i = 0; i < np; ++i) {
    matrix_v y_i = y_m.block(0, icol, y0.size(), len[i]);
    torsten::test::test_grad(y0_var, y0_var_m[i], y, y_i, 1e-16, 1e-16);
    icol += len[i];
  }
}

TEST_F(TorstenOdeTest_sho, group_rk45_fwd_sensitivity_y0) {
  const int np = 5;

  vector<var> y0_var = stan::math::to_var(y0);

  vector<int> len(np, ts.size());
  vector<double> ts_m;
  ts_m.reserve(np * ts.size());
  for (int i = 0; i < np; ++i) ts_m.insert(ts_m.end(), ts.begin(), ts.end());

  vector<vector<double> > theta_m (np, theta);
  vector<vector<var> > y0_var_m (np, stan::math::to_var(y0));
  vector<vector<double> > x_r_m (np, x_r);
  vector<vector<int> > x_i_m (np, x_i);

  vector<vector<var> > y = integrate_ode_rk45(f, y0_var, t0, ts, theta, x_r, x_i);
  Eigen::Matrix<var, -1, -1> y_m = pmx_integrate_ode_group_rk45(f, y0_var_m, t0, len, ts_m, theta_m , x_r_m, x_i_m);

  EXPECT_EQ(y_m.cols(), ts_m.size());
  int icol = 0;
  for (int i = 0; i < np; ++i) {
    matrix_v y_i = y_m.block(0, icol, y0.size(), len[i]);
    torsten::test::test_grad(y0_var, y0_var_m[i], y, y_i, 1e-16, 1e-16);
    icol += len[i];
  }
}

TEST_F(TorstenOdeTest_sho, group_rk45_fwd_sensitivity_y0_theta) {
  const int np = 10;

  vector<var> y0_var = stan::math::to_var(y0);
  vector<var> theta_var = stan::math::to_var(theta);
  vector<var> vars = y0_var;
  vars.insert(vars.end(), theta_var.begin(), theta_var.end());

  vector<int> len(np, ts.size());
  vector<double> ts_m;
  ts_m.reserve(np * ts.size());
  for (int i = 0; i < np; ++i) ts_m.insert(ts_m.end(), ts.begin(), ts.end());

  vector<vector<var> > theta_var_m (np, stan::math::to_var(theta));
  vector<vector<var> > y0_var_m (np, stan::math::to_var(y0));
  vector<vector<double> > x_r_m (np, x_r);
  vector<vector<int> > x_i_m (np, x_i);

  vector<vector<var> > y = integrate_ode_rk45(f, y0_var, t0, ts, theta_var, x_r, x_i);
  Eigen::Matrix<var, -1, -1> y_m = pmx_integrate_ode_group_rk45(f, y0_var_m, t0, len, ts_m, theta_var_m , x_r_m, x_i_m);

  EXPECT_EQ(y_m.cols(), ts_m.size());
  int icol = 0;
  for (int i = 0; i < np; ++i) {
    matrix_v y_i = y_m.block(0, icol, y0.size(), len[i]);
    vector<var> vars_i = y0_var_m[i];
    vars_i.insert(vars_i.end(), theta_var_m[i].begin(), theta_var_m[i].end());
    torsten::test::test_grad(vars, vars_i, y, y_i, 1e-14, 1e-14);
    icol += len[i];
  }
}

TEST_F(TorstenOdeTest_chem, group_rk45_fwd_sensitivity_y0_theta) {
  const int np = 10;

  vector<var> y0_var = stan::math::to_var(y0);
  vector<var> theta_var = stan::math::to_var(theta);
  vector<var> vars = y0_var;
  vars.insert(vars.end(), theta_var.begin(), theta_var.end());

  vector<int> len(np, ts.size());
  vector<double> ts_m;
  ts_m.reserve(np * ts.size());
  for (int i = 0; i < np; ++i) ts_m.insert(ts_m.end(), ts.begin(), ts.end());

  vector<vector<var> > theta_var_m (np, stan::math::to_var(theta));
  vector<vector<var> > y0_var_m (np, stan::math::to_var(y0));
  vector<vector<double> > x_r_m (np, x_r);
  vector<vector<int> > x_i_m (np, x_i);

  vector<vector<var> > y = integrate_ode_rk45(f, y0_var, t0, ts, theta_var, x_r, x_i);
  Eigen::Matrix<var, -1, -1> y_m = pmx_integrate_ode_group_rk45(f, y0_var_m, t0, len, ts_m, theta_var_m , x_r_m, x_i_m);

  EXPECT_EQ(y_m.cols(), ts_m.size());
  int icol = 0;
  for (int i = 0; i < np; ++i) {
    matrix_v y_i = y_m.block(0, icol, y0.size(), len[i]);
    vector<var> vars_i = y0_var_m[i];
    vars_i.insert(vars_i.end(), theta_var_m[i].begin(), theta_var_m[i].end());
    torsten::test::test_grad(vars, vars_i, y, y_i, 2e-14, 1e-12);
    icol += len[i];
  }
}

TEST_F(TorstenOdeTest_neutropenia, group_rk45_fwd_sensitivity_y0_theta) {
  const int np = 10;

  vector<var> y0_var = stan::math::to_var(y0);
  vector<var> theta_var = stan::math::to_var(theta);
  vector<var> vars = y0_var;
  vars.insert(vars.end(), theta_var.begin(), theta_var.end());

  vector<int> len(np, ts.size());
  vector<double> ts_m;
  ts_m.reserve(np * ts.size());
  for (int i = 0; i < np; ++i) ts_m.insert(ts_m.end(), ts.begin(), ts.end());

  vector<vector<var> > theta_var_m (np, stan::math::to_var(theta));
  vector<vector<var> > y0_var_m (np, stan::math::to_var(y0));
  vector<vector<double> > x_r_m (np, x_r);
  vector<vector<int> > x_i_m (np, x_i);

  vector<vector<var> > y = integrate_ode_rk45(f, y0_var, t0, ts, theta_var, x_r, x_i);
  Eigen::Matrix<var, -1, -1> y_m = pmx_integrate_ode_group_rk45(f, y0_var_m, t0, len, ts_m, theta_var_m , x_r_m, x_i_m);

  EXPECT_EQ(y_m.cols(), ts_m.size());
  int icol = 0;
  for (int i = 0; i < np; ++i) {
    matrix_v y_i = y_m.block(0, icol, y0.size(), len[i]);
    vector<var> vars_i = y0_var_m[i];
    vars_i.insert(vars_i.end(), theta_var_m[i].begin(), theta_var_m[i].end());
    torsten::test::test_grad(vars, vars_i, y, y_i, 1e-12, 1e-9);
    icol += len[i];
  }
}
