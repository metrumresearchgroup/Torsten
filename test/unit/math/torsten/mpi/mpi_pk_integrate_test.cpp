#include <stan/math.hpp>
#include <stan/math/rev/core.hpp>
#include <stan/math/rev/mat/fun/typedefs.hpp>
#include <test/unit/math/rev/mat/fun/util.hpp>
#include <stan/math/torsten/dsolve/pmx_integrate_ode_adams.hpp>
#include <stan/math/torsten/dsolve/pmx_integrate_ode_bdf.hpp>
#include <test/unit/math/torsten/pmx_ode_test_fixture.hpp>
#include <test/unit/math/torsten/test_util.hpp>
#include <stan/math/rev/mat/functor/integrate_ode_bdf.hpp>
#include <nvector/nvector_serial.h>
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

using torsten::dsolve::PMXCvodesFwdSystem;
using stan::math::integrate_ode_bdf;
using stan::math::integrate_ode_adams;
using torsten::pmx_integrate_ode_group_bdf;
using torsten::pmx_integrate_ode_group_adams;
using stan::math::matrix_v;
using stan::math::var;
using std::vector;

TEST_F(TorstenOdeTest_chem, cvodes_ivp_system_bdf_mpi) {
  torsten::mpi::Envionment::init();

  const int np = 10;

  vector<int> len(np, ts.size());
  vector<double> ts_m;
  ts_m.reserve(np * ts.size());
  for (int i = 0; i < np; ++i) ts_m.insert(ts_m.end(), ts.begin(), ts.end());

  vector<vector<double> > y0_m (np, y0);
  vector<vector<double> > theta_m (np, theta);
  vector<vector<double> > x_r_m (np, x_r);
  vector<vector<int> > x_i_m (np, x_i);

  vector<vector<double> > y = integrate_ode_bdf(f, y0, t0, ts, theta , x_r, x_i); // NOLINT
  Eigen::MatrixXd y_m = pmx_integrate_ode_group_bdf(f, y0_m, t0, len, ts_m, theta_m , x_r_m, x_i_m);

  EXPECT_EQ(y_m.cols(), ts_m.size());
  int icol = 0;
  for (size_t i = 0; i < np; ++i) {
    Eigen::MatrixXd y_i = y_m.block(0, icol, y0.size(), len[i]);
    torsten::test::test_val(y, y_i);
    icol += len[i];
  }
}

TEST_F(TorstenOdeTest_neutropenia, fwd_sensitivity_group_theta_adams) {
  const int np = 10;
  std::vector<double> ts0 {ts};

  // the first two param are group params
  vector<var> theta_var_all(theta.begin(), theta.end());
  vector<var> theta_var(theta_var_all.begin() + 2, theta_var_all.end());
  vector<var> group_theta_var(theta_var_all.begin(), theta_var_all.begin() + 2);
  assert(theta_var.size() + group_theta_var.size() == theta.size());

  vector<int> len(np, ts0.size());
  vector<double> ts_m;
  ts_m.reserve(np * ts0.size());
  for (int i = 0; i < np; ++i) ts_m.insert(ts_m.end(), ts0.begin(), ts0.end());
  
  vector<vector<double> > y0_m (np, y0);
  vector<vector<var> > theta_var_m (np, theta_var);
  vector<vector<var> > theta_var_all_m (np, theta_var_all);
  vector<vector<double> > x_r_m (np, x_r);
  vector<vector<int> > x_i_m (np, x_i);

  Eigen::Matrix<var, -1, -1> y_m1 = pmx_integrate_ode_group_adams(f, y0_m, t0, len, ts_m, group_theta_var,theta_var_m , x_r_m, x_i_m);
  Eigen::Matrix<var, -1, -1> y_m2 = pmx_integrate_ode_group_adams(f, y0_m, t0, len, ts_m, theta_var_all_m , x_r_m, x_i_m);

  torsten::test::test_grad(theta_var_all, y_m1, y_m2, 1e-10, 1e-10);
}

TEST_F(TorstenOdeTest_neutropenia, ivp_group_theta_adams) {
  const int np = 10;
  std::vector<double> ts0 {ts};

  // the first two param are group params
  vector<double> theta_i(theta.begin() + 2, theta.end());
  vector<double> group_theta(theta.begin(), theta.begin() + 2);
  assert(theta_i.size() + group_theta.size() == theta.size());

  vector<int> len(np, ts0.size());
  vector<double> ts_m;
  ts_m.reserve(np * ts0.size());
  for (int i = 0; i < np; ++i) ts_m.insert(ts_m.end(), ts0.begin(), ts0.end());
  
  vector<vector<double> > y0_m (np, y0);
  vector<vector<double> > theta_m (np, theta_i);
  vector<vector<double> > theta_all_m (np, theta);
  vector<vector<double> > x_r_m (np, x_r);
  vector<vector<int> > x_i_m (np, x_i);

  Eigen::MatrixXd y_m1 = pmx_integrate_ode_group_adams(f, y0_m, t0, len, ts_m, group_theta, theta_m, x_r_m, x_i_m);
  Eigen::MatrixXd y_m2 = pmx_integrate_ode_group_adams(f, y0_m, t0, len, ts_m, theta_all_m, x_r_m, x_i_m);

  torsten::test::test_val(y_m1, y_m2);
}

TEST_F(TorstenOdeTest_neutropenia, fwd_sensitivity_group_theta_bdf) {
  const int np = 10;
  std::vector<double> ts0 {ts};

  // the first two param are group params
  vector<var> theta_var_all(theta.begin(), theta.end());
  vector<var> theta_var(theta_var_all.begin() + 2, theta_var_all.end());
  vector<var> group_theta_var(theta_var_all.begin(), theta_var_all.begin() + 2);
  assert(theta_var.size() + group_theta_var.size() == theta.size());

  vector<int> len(np, ts0.size());
  vector<double> ts_m;
  ts_m.reserve(np * ts0.size());
  for (int i = 0; i < np; ++i) ts_m.insert(ts_m.end(), ts0.begin(), ts0.end());
  
  vector<vector<double> > y0_m (np, y0);
  vector<vector<var> > theta_var_m (np, theta_var);
  vector<vector<var> > theta_var_all_m (np, theta_var_all);
  vector<vector<double> > x_r_m (np, x_r);
  vector<vector<int> > x_i_m (np, x_i);

  Eigen::Matrix<var, -1, -1> y_m1 = pmx_integrate_ode_group_bdf(f, y0_m, t0, len, ts_m, group_theta_var,theta_var_m , x_r_m, x_i_m);
  Eigen::Matrix<var, -1, -1> y_m2 = pmx_integrate_ode_group_bdf(f, y0_m, t0, len, ts_m, theta_var_all_m , x_r_m, x_i_m);

  torsten::test::test_grad(theta_var_all, y_m1, y_m2, 1e-10, 1e-10);
}

TEST_F(TorstenOdeTest_neutropenia, ivp_group_theta_bdf) {
  const int np = 10;
  std::vector<double> ts0 {ts};

  // the first two param are group params
  vector<double> theta_i(theta.begin() + 2, theta.end());
  vector<double> group_theta(theta.begin(), theta.begin() + 2);
  assert(theta_i.size() + group_theta.size() == theta.size());

  vector<int> len(np, ts0.size());
  vector<double> ts_m;
  ts_m.reserve(np * ts0.size());
  for (int i = 0; i < np; ++i) ts_m.insert(ts_m.end(), ts0.begin(), ts0.end());
  
  vector<vector<double> > y0_m (np, y0);
  vector<vector<double> > theta_m (np, theta_i);
  vector<vector<double> > theta_all_m (np, theta);
  vector<vector<double> > x_r_m (np, x_r);
  vector<vector<int> > x_i_m (np, x_i);

  Eigen::MatrixXd y_m1 = pmx_integrate_ode_group_bdf(f, y0_m, t0, len, ts_m, group_theta, theta_m, x_r_m, x_i_m);
  Eigen::MatrixXd y_m2 = pmx_integrate_ode_group_bdf(f, y0_m, t0, len, ts_m, theta_all_m, x_r_m, x_i_m);

  torsten::test::test_val(y_m1, y_m2);
}

TEST_F(TorstenOdeTest_chem, fwd_sensitivity_theta_AD_bdf_mpi) {
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

  vector<vector<var> > y1 = stan::math::integrate_ode_bdf(f, y0, t0, ts, theta_var1, x_r, x_i);
  vector<vector<var> > y2 = stan::math::integrate_ode_bdf(f, y0, t0, ts, theta_var2, x_r, x_i);
  Eigen::Matrix<var, -1, -1> y_m = pmx_integrate_ode_group_bdf(f, y0_m, t0, len, ts_m, theta_var_m , x_r_m, x_i_m);

  EXPECT_EQ(y_m.cols(), ts_m.size());
  int icol = 0;
  matrix_v y;
  y = y_m.block(0, icol, y0.size(), len[0]);
  torsten::test::test_grad(theta_var1, theta_var_m[0], y1, y, 1e-7, 1e-7);

  icol += len[0];
  y = y_m.block(0, icol, y0.size(), len[1]);
  torsten::test::test_grad(theta_var2, theta_var_m[1], y2, y, 1e-7, 1e-7);
}

TEST_F(TorstenOdeTest_chem, fwd_sensitivity_theta_AD_adams_mpi) {
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

  vector<vector<var> > y1 = stan::math::integrate_ode_adams(f, y0, t0, ts, theta_var1, x_r, x_i);
  vector<vector<var> > y2 = stan::math::integrate_ode_adams(f, y0, t0, ts, theta_var2, x_r, x_i);
  Eigen::Matrix<var, -1, -1> y_m = pmx_integrate_ode_group_adams(f, y0_m, t0, len, ts_m, theta_var_m , x_r_m, x_i_m);

  EXPECT_EQ(y_m.cols(), ts_m.size());
  int icol = 0;
  matrix_v y;
  y = y_m.block(0, icol, y0.size(), len[0]);
  torsten::test::test_grad(theta_var1, theta_var_m[0], y1, y, 1e-7, 1e-7);

  icol += len[0];
  y = y_m.block(0, icol, y0.size(), len[1]);
  torsten::test::test_grad(theta_var2, theta_var_m[1], y2, y, 1e-7, 1e-7);
}

TEST_F(TorstenOdeTest_chem, fwd_sensitivity_theta_AD_bdf_mpi_performance) {
  const int np = 10;
  std::vector<double> ts0 {ts};
  ts0.push_back(400);

  vector<var> theta_var = stan::math::to_var(theta);

  vector<int> len(np, ts0.size());
  vector<double> ts_m;
  ts_m.reserve(np * ts0.size());
  for (int i = 0; i < np; ++i) ts_m.insert(ts_m.end(), ts0.begin(), ts0.end());

  vector<vector<double> > y0_m (np, y0);
  vector<vector<var> > theta_var_m (np, stan::math::to_var(theta));

  vector<vector<double> > x_r_m (np, x_r);
  vector<vector<int> > x_i_m (np, x_i);

  vector<vector<var> > y = stan::math::integrate_ode_bdf(f, y0, t0, ts0, theta_var, x_r, x_i);
  Eigen::Matrix<var, -1, -1> y_m = pmx_integrate_ode_group_bdf(f, y0_m, t0, len, ts_m, theta_var_m , x_r_m, x_i_m);

  EXPECT_EQ(y_m.cols(), ts_m.size());
  int icol = 0;
  for (int i = 0; i < np; ++i) {
    matrix_v y_i = y_m.block(0, icol, y0.size(), len[i]);
    torsten::test::test_grad(theta_var, theta_var_m[i], y, y_i, 1e-7, 3e-7);
    icol += len[i];
  }
}

TEST_F(TorstenOdeTest_chem, fwd_sensitivity_theta_AD_adams_mpi_performance) {
  const int np = 10;
  std::vector<double> ts0 {ts};
  ts0.push_back(400);

  vector<var> theta_var = stan::math::to_var(theta);

  vector<int> len(np, ts0.size());
  vector<double> ts_m;
  ts_m.reserve(np * ts0.size());
  for (int i = 0; i < np; ++i) ts_m.insert(ts_m.end(), ts0.begin(), ts0.end());

  vector<vector<double> > y0_m (np, y0);
  vector<vector<var> > theta_var_m (np, stan::math::to_var(theta));

  vector<vector<double> > x_r_m (np, x_r);
  vector<vector<int> > x_i_m (np, x_i);

  vector<vector<var> > y = stan::math::integrate_ode_adams(f, y0, t0, ts0, theta_var, x_r, x_i);
  Eigen::Matrix<var, -1, -1> y_m = pmx_integrate_ode_group_adams(f, y0_m, t0, len, ts_m, theta_var_m , x_r_m, x_i_m);

  EXPECT_EQ(y_m.cols(), ts_m.size());
  int icol = 0;
  for (int i = 0; i < np; ++i) {
    matrix_v y_i = y_m.block(0, icol, y0.size(), len[i]);
    torsten::test::test_grad(theta_var, theta_var_m[i], y, y_i, 1e-7, 3e-7);
    icol += len[i];
  }
}

TEST_F(TorstenOdeTest_lorenz, fwd_sensitivity_theta_AD_bdf_mpi_performance) {
  const int np = 10;
  std::vector<double> ts0 {ts};

  vector<var> theta_var = stan::math::to_var(theta);

  vector<int> len(np, ts0.size());
  vector<double> ts_m;
  ts_m.reserve(np * ts0.size());
  for (int i = 0; i < np; ++i) ts_m.insert(ts_m.end(), ts0.begin(), ts0.end());

  vector<vector<double> > y0_m (np, y0);
  vector<vector<var> > theta_var_m (np, stan::math::to_var(theta));

  vector<vector<double> > x_r_m (np, x_r);
  vector<vector<int> > x_i_m (np, x_i);

  vector<vector<var> > y = stan::math::integrate_ode_bdf(f, y0, t0, ts0, theta_var, x_r, x_i);
  Eigen::Matrix<var, -1, -1> y_m = pmx_integrate_ode_group_bdf(f, y0_m, t0, len, ts_m, theta_var_m , x_r_m, x_i_m);

  EXPECT_EQ(y_m.cols(), ts_m.size());
  int icol = 0;
  for (int i = 0; i < np; ++i) {
    matrix_v y_i = y_m.block(0, icol, y0.size(), len[i]);
    torsten::test::test_grad(theta_var, theta_var_m[i], y, y_i, 1e-7, 3e-7);
    icol += len[i];
  }
}

TEST_F(TorstenOdeTest_lorenz, fwd_sensitivity_theta_AD_adams_mpi_performance) {
  const int np = 10;
  std::vector<double> ts0 {ts};

  vector<var> theta_var = stan::math::to_var(theta);

  vector<int> len(np, ts0.size());
  vector<double> ts_m;
  ts_m.reserve(np * ts0.size());
  for (int i = 0; i < np; ++i) ts_m.insert(ts_m.end(), ts0.begin(), ts0.end());

  vector<vector<double> > y0_m (np, y0);
  vector<vector<var> > theta_var_m (np, stan::math::to_var(theta));

  vector<vector<double> > x_r_m (np, x_r);
  vector<vector<int> > x_i_m (np, x_i);

  vector<vector<var> > y = stan::math::integrate_ode_adams(f, y0, t0, ts0, theta_var, x_r, x_i);
  Eigen::Matrix<var, -1, -1> y_m = pmx_integrate_ode_group_adams(f, y0_m, t0, len, ts_m, theta_var_m , x_r_m, x_i_m);

  EXPECT_EQ(y_m.cols(), ts_m.size());
  int icol = 0;
  for (int i = 0; i < np; ++i) {
    matrix_v y_i = y_m.block(0, icol, y0.size(), len[i]);
    torsten::test::test_grad(theta_var, theta_var_m[i], y, y_i, 1e-6, 5e-5);
    icol += len[i];
  }
}

TEST_F(TorstenOdeTest_neutropenia, fwd_sensitivity_theta_AD_bdf_mpi) {
  torsten::mpi::Envionment::init();

  // size of population
  const int np = 10;

  vector<var> theta_var = stan::math::to_var(theta);

  vector<int> len(np, ts.size());
  vector<double> ts_m;
  ts_m.reserve(np * ts.size());
  for (int i = 0; i < np; ++i) ts_m.insert(ts_m.end(), ts.begin(), ts.end());

  vector<vector<double> > y0_m (np, y0);
  vector<vector<var> > theta_var_m (np, stan::math::to_var(theta));
  vector<vector<double> > x_r_m (np, x_r);
  vector<vector<int> > x_i_m (np, x_i);

  vector<vector<var> > y = stan::math::integrate_ode_bdf(f, y0, t0, ts, theta_var, x_r, x_i);
  Eigen::Matrix<var, -1, -1> y_m = pmx_integrate_ode_group_bdf(f, y0_m, t0, len, ts_m, theta_var_m , x_r_m, x_i_m);

  EXPECT_EQ(y_m.cols(), ts_m.size());
  int icol = 0;
  for (int i = 0; i < np; ++i) {
    matrix_v y_i = y_m.block(0, icol, y0.size(), len[i]);
    torsten::test::test_grad(theta_var, theta_var_m[i], y, y_i, 1e-7, 2e-7);
    icol += len[i];
  }
}

TEST_F(TorstenOdeTest_neutropenia, fwd_sensitivity_theta_AD_adams_mpi) {
  torsten::mpi::Envionment::init();

  // size of population
  const int np = 10;

  vector<var> theta_var = stan::math::to_var(theta);

  vector<int> len(np, ts.size());
  vector<double> ts_m;
  ts_m.reserve(np * ts.size());
  for (int i = 0; i < np; ++i) ts_m.insert(ts_m.end(), ts.begin(), ts.end());

  vector<vector<double> > y0_m (np, y0);
  vector<vector<var> > theta_var_m (np, stan::math::to_var(theta));
  vector<vector<double> > x_r_m (np, x_r);
  vector<vector<int> > x_i_m (np, x_i);

  vector<vector<var> > y = stan::math::integrate_ode_adams(f, y0, t0, ts, theta_var, x_r, x_i);
  Eigen::Matrix<var, -1, -1> y_m = pmx_integrate_ode_group_adams(f, y0_m, t0, len, ts_m, theta_var_m , x_r_m, x_i_m);

  EXPECT_EQ(y_m.cols(), ts_m.size());
  int icol = 0;
  for (int i = 0; i < np; ++i) {
    matrix_v y_i = y_m.block(0, icol, y0.size(), len[i]);
    torsten::test::test_grad(theta_var, theta_var_m[i], y, y_i, 1e-7, 2e-7);
    icol += len[i];
  }
}

TEST_F(TorstenOdeTest_neutropenia, mpi_rank_exception_data_only) {
  const int np = 5;

  vector<int> len(np, ts.size());
  vector<double> ts_m;
  ts_m.reserve(np * ts.size());
  for (int i = 0; i < np; ++i) ts_m.insert(ts_m.end(), ts.begin(), ts.end());

  vector<vector<double> > y0_m (np, y0);
  vector<vector<double> > theta_m (np, theta);
  vector<vector<double> > x_r_m (np, x_r);
  vector<vector<int> > x_i_m (np, x_i);

  torsten::mpi::Envionment::init();

#ifdef TORSTEN_MPI
  MPI_Comm comm;
  comm = MPI_COMM_WORLD;
  int rank, size;
  MPI_Comm_size(comm, &size);
  MPI_Comm_rank(comm, &rank);
#endif

  int id = 1;

  theta_m[id][0] = std::numeric_limits<double>::infinity();
  EXPECT_THROW(torsten::pmx_integrate_ode_group_bdf(f, y0_m, t0, len, ts_m, theta_m , x_r_m, x_i_m),
               std::exception);
  EXPECT_THROW(torsten::pmx_integrate_ode_group_adams(f, y0_m, t0, len, ts_m, theta_m , x_r_m, x_i_m),
               std::exception);
#ifdef TORSTEN_MPI
  MPI_Barrier(comm);
#endif

  theta_m[id][0] = -1.0E+30;
  theta_m[id][1] = -1.0E+30;
  theta_m[id][4] = -1.0E+30;
  EXPECT_THROW(torsten::pmx_integrate_ode_group_bdf(f, y0_m, t0, len, ts_m, theta_m , x_r_m, x_i_m,
                                                     0, 1e-6, 1e-6, 1e2),
               std::exception);
  EXPECT_THROW(torsten::pmx_integrate_ode_group_adams(f, y0_m, t0, len, ts_m, theta_m , x_r_m, x_i_m,
                                                       0, 1e-6, 1e-6, 1e2),
               std::exception);
#ifdef TORSTEN_MPI
  MPI_Barrier(comm);
#endif

  theta_m[id][0] = std::numeric_limits<double>::quiet_NaN();
  EXPECT_THROW(torsten::pmx_integrate_ode_group_bdf(f, y0_m, t0, len, ts_m, theta_m , x_r_m, x_i_m),
               std::exception);
  EXPECT_THROW(torsten::pmx_integrate_ode_group_bdf(f, y0_m, t0, len, ts_m, theta_m , x_r_m, x_i_m),
               std::exception);
#ifdef TORSTEN_MPI
  MPI_Barrier(comm);
#endif

  id = 1;
  theta_m[id][0] = std::numeric_limits<double>::infinity();
  id = 4;
  theta_m[id][0] = -1.0E+30;
  theta_m[id][1] = -1.0E+30;
  theta_m[id][4] = -1.0E+30;
  EXPECT_THROW(torsten::pmx_integrate_ode_group_bdf(f, y0_m, t0, len, ts_m, theta_m , x_r_m, x_i_m),
               std::exception);
#ifdef TORSTEN_MPI
  MPI_Barrier(comm);
#endif

  id = 1;
  theta_m[id][0] = -1.0E+30;
  theta_m[id][1] = -1.0E+30;
  theta_m[id][4] = -1.0E+30;
  id = 4;
  theta_m[id][0] = std::numeric_limits<double>::infinity();
  EXPECT_THROW(torsten::pmx_integrate_ode_group_bdf(f, y0_m, t0, len, ts_m, theta_m , x_r_m, x_i_m,
                                                     0, 1e-6, 1e-6, 1e2),
               std::exception);
#ifdef TORSTEN_MPI
  MPI_Barrier(comm);
#endif
}

TEST_F(TorstenOdeTest_neutropenia, mpi_rank_exception_par_var) {
  const int np = 5;

  vector<var> theta_var(stan::math::to_var(theta));

  vector<int> len(np, ts.size());
  vector<double> ts_m;
  ts_m.reserve(np * ts.size());
  for (int i = 0; i < np; ++i) ts_m.insert(ts_m.end(), ts.begin(), ts.end());

  vector<vector<double> > y0_m (np, y0);
  vector<vector<var> > theta_m_v (np, theta_var);
  vector<vector<double> > x_r_m (np, x_r);
  vector<vector<int> > x_i_m (np, x_i);

  torsten::mpi::Envionment::init();

#ifdef TORSTEN_MPI
  MPI_Comm comm;
  comm = MPI_COMM_WORLD;
  int rank, size;
  MPI_Comm_size(comm, &size);
  MPI_Comm_rank(comm, &rank);
#endif

  int id = 1;

  theta_m_v[id][0] = std::numeric_limits<double>::infinity();
  EXPECT_THROW(torsten::pmx_integrate_ode_group_bdf(f, y0_m, t0, len, ts_m, theta_m_v , x_r_m, x_i_m),
               std::exception);
  EXPECT_THROW(torsten::pmx_integrate_ode_group_adams(f, y0_m, t0, len, ts_m, theta_m_v , x_r_m, x_i_m),
               std::exception);
#ifdef TORSTEN_MPI
  MPI_Barrier(comm);
#endif

  theta_m_v[id][0] = -1.0E+30;
  theta_m_v[id][1] = -1.0E+30;
  theta_m_v[id][4] = -1.0E+30;
  EXPECT_THROW(torsten::pmx_integrate_ode_group_bdf(f, y0_m, t0, len, ts_m, theta_m_v , x_r_m, x_i_m,
                                                     0, 1e-6, 1e-6, 1e2),
               std::exception);
  EXPECT_THROW(torsten::pmx_integrate_ode_group_adams(f, y0_m, t0, len, ts_m, theta_m_v , x_r_m, x_i_m,
                                                       0, 1e-6, 1e-6, 1e2),
               std::exception);
#ifdef TORSTEN_MPI
  MPI_Barrier(comm);
#endif

  theta_m_v[id][0] = std::numeric_limits<double>::quiet_NaN();
  EXPECT_THROW(torsten::pmx_integrate_ode_group_bdf(f, y0_m, t0, len, ts_m, theta_m_v , x_r_m, x_i_m),
               std::exception);
  EXPECT_THROW(torsten::pmx_integrate_ode_group_bdf(f, y0_m, t0, len, ts_m, theta_m_v , x_r_m, x_i_m),
               std::exception);
#ifdef TORSTEN_MPI
  MPI_Barrier(comm);
#endif

  id = 1;
  theta_m_v[id][0] = std::numeric_limits<double>::infinity();
  id = 4;
  theta_m_v[id][0] = -1.0E+30;
  theta_m_v[id][1] = -1.0E+30;
  theta_m_v[id][4] = -1.0E+30;
  EXPECT_THROW(torsten::pmx_integrate_ode_group_bdf(f, y0_m, t0, len, ts_m, theta_m_v , x_r_m, x_i_m),
               std::exception);
#ifdef TORSTEN_MPI
  MPI_Barrier(comm);
#endif

  id = 1;
  theta_m_v[id][0] = -1.0E+30;
  theta_m_v[id][1] = -1.0E+30;
  theta_m_v[id][4] = -1.0E+30;
  id = 4;
  theta_m_v[id][0] = std::numeric_limits<double>::infinity();
  EXPECT_THROW(torsten::pmx_integrate_ode_group_bdf(f, y0_m, t0, len, ts_m, theta_m_v , x_r_m, x_i_m,
                                                     0, 1e-6, 1e-6, 1e2),
               std::exception);
#ifdef TORSTEN_MPI
  MPI_Barrier(comm);
#endif
}
