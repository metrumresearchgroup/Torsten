#include <stan/math.hpp>
#include <stan/math/rev/core.hpp>
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
#include <vector>
#include <limits>
#include <vector>
#include <chrono>
#include <ctime>



TEST_F(TorstenOdeTest_chem, cvodes_ivp_system_bdf_mpi) {
  using torsten::dsolve::PMXCvodesFwdSystem;
  using stan::math::integrate_ode_bdf;
  using torsten::dsolve::pmx_integrate_ode_bdf;
  using std::vector;

  torsten::mpi::Envionment::init();

  const int np = 10;
  vector<vector<double> > ts_m(np, ts);
  vector<vector<double> > y0_m (np, y0);
  vector<vector<double> > theta_m (np, theta);
  vector<vector<double> > x_r_m (np, x_r);
  vector<vector<int> > x_i_m (np, x_i);

  vector<vector<double> > y = integrate_ode_bdf(f, y0, t0, ts, theta , x_r, x_i); // NOLINT
  vector<Eigen::MatrixXd> y_m = pmx_integrate_ode_bdf(f, y0_m, t0, ts_m, theta_m , x_r_m, x_i_m);
  EXPECT_EQ(y_m.size(), theta_m.size());
  for (size_t i = 0; i < y_m.size(); ++i) {
    torsten::test::test_val(y, y_m[i]);
  }
}

TEST_F(TorstenOdeTest_chem, fwd_sensitivity_theta_AD_bdf_mpi) {
  using torsten::dsolve::PMXCvodesFwdSystem;
  using torsten::dsolve::pmx_integrate_ode_bdf;
  using stan::math::var;
  using std::vector;

  std::vector<var> theta_var1 = stan::math::to_var(theta);
  std::vector<var> theta_var2 = stan::math::to_var(theta);
  theta_var2[0] *= 1.2;

  vector<vector<var> > theta_var_m(2);
  theta_var_m[0] = stan::math::to_var(theta);
  theta_var_m[1] = stan::math::to_var(theta);
  theta_var_m[1][0] *= 1.2;

  vector<vector<double> > ts_m {ts, ts};
  vector<vector<double> > y0_m {y0, y0};

  vector<vector<double> > x_r_m {x_r, x_r};
  vector<vector<int> > x_i_m {x_i, x_i};

  vector<vector<var> > y1 = stan::math::integrate_ode_bdf(f, y0, t0, ts, theta_var1, x_r, x_i);
  vector<vector<var> > y2 = stan::math::integrate_ode_bdf(f, y0, t0, ts, theta_var2, x_r, x_i);
  vector<Eigen::Matrix<var, -1, -1> > y_m = pmx_integrate_ode_bdf(f, y0_m, t0, ts_m, theta_var_m , x_r_m, x_i_m);

  torsten::test::test_grad(theta_var1, theta_var_m[0], y1, y_m[0], 1e-7, 1e-7);
  torsten::test::test_grad(theta_var2, theta_var_m[1], y2, y_m[1], 1e-7, 1e-7);
}

TEST_F(TorstenOdeTest_chem, fwd_sensitivity_theta_AD_adams_mpi) {
  using torsten::dsolve::PMXCvodesFwdSystem;
  using torsten::dsolve::pmx_integrate_ode_adams;
  using stan::math::var;
  using std::vector;

  std::vector<var> theta_var1 = stan::math::to_var(theta);
  std::vector<var> theta_var2 = stan::math::to_var(theta);
  theta_var2[0] *= 1.2;

  vector<vector<var> > theta_var_m(2);
  theta_var_m[0] = stan::math::to_var(theta);
  theta_var_m[1] = stan::math::to_var(theta);
  theta_var_m[1][0] *= 1.2;

  vector<vector<double> > ts_m {ts, ts};
  vector<vector<double> > y0_m {y0, y0};

  vector<vector<double> > x_r_m {x_r, x_r};
  vector<vector<int> > x_i_m {x_i, x_i};

  vector<vector<var> > y1 = stan::math::integrate_ode_adams(f, y0, t0, ts, theta_var1, x_r, x_i);
  vector<vector<var> > y2 = stan::math::integrate_ode_adams(f, y0, t0, ts, theta_var2, x_r, x_i);
  vector<Eigen::Matrix<var, -1, -1> > y_m = pmx_integrate_ode_adams(f, y0_m, t0, ts_m, theta_var_m , x_r_m, x_i_m);

  torsten::test::test_grad(theta_var1, theta_var_m[0], y1, y_m[0], 1e-7, 5e-6);
  torsten::test::test_grad(theta_var2, theta_var_m[1], y2, y_m[1], 1e-7, 5e-6);
}

TEST_F(TorstenOdeTest_chem, fwd_sensitivity_theta_AD_bdf_mpi_performance) {
  using torsten::dsolve::PMXCvodesFwdSystem;
  using torsten::dsolve::pmx_integrate_ode_bdf;
  using stan::math::var;
  using std::vector;

  // size of population
  const int np = 10;
  std::vector<double> ts0 {ts};
  ts0.push_back(400);

  vector<var> theta_var = stan::math::to_var(theta);

  vector<vector<double> > ts_m (np, ts0);
  vector<vector<double> > y0_m (np, y0);
  vector<vector<var> > theta_var_m (np, stan::math::to_var(theta));

  vector<vector<double> > x_r_m (np, x_r);
  vector<vector<int> > x_i_m (np, x_i);

  vector<vector<var> > y = stan::math::integrate_ode_bdf(f, y0, t0, ts0, theta_var, x_r, x_i);
  vector<Eigen::Matrix<var, -1, -1> > y_m = pmx_integrate_ode_bdf(f, y0_m, t0, ts_m, theta_var_m , x_r_m, x_i_m);

  for (int i = 0; i < np; ++i) {
    torsten::test::test_grad(theta_var, theta_var_m[i], y, y_m[i], 1e-7, 3e-7);
  }
}

TEST_F(TorstenOdeTest_chem, fwd_sensitivity_theta_AD_adams_mpi_performance) {
  using torsten::dsolve::PMXCvodesFwdSystem;
  using torsten::dsolve::pmx_integrate_ode_adams;
  using stan::math::var;
  using std::vector;

  // size of population
  const int np = 10;
  std::vector<double> ts0 {ts};
  ts0.push_back(400);

  vector<var> theta_var = stan::math::to_var(theta);

  vector<vector<double> > ts_m (np, ts0);
  vector<vector<double> > y0_m (np, y0);
  vector<vector<var> > theta_var_m (np, stan::math::to_var(theta));

  vector<vector<double> > x_r_m (np, x_r);
  vector<vector<int> > x_i_m (np, x_i);

  vector<vector<var> > y = stan::math::integrate_ode_adams(f, y0, t0, ts0, theta_var, x_r, x_i);
  vector<Eigen::Matrix<var, -1, -1> > y_m = pmx_integrate_ode_adams(f, y0_m, t0, ts_m, theta_var_m , x_r_m, x_i_m);

  for (int i = 0; i < np; ++i) {
    torsten::test::test_grad(theta_var, theta_var_m[i], y, y_m[i], 1e-7, 5e-6);
  }
}

TEST_F(TorstenOdeTest_lorenz, fwd_sensitivity_theta_AD_bdf_mpi_performance) {
  using torsten::dsolve::PMXCvodesFwdSystem;
  using torsten::dsolve::pmx_integrate_ode_bdf;
  using stan::math::var;
  using std::vector;

  // size of population
  const int np = 10;
  std::vector<double> ts0 {ts};

  vector<var> theta_var = stan::math::to_var(theta);

  vector<vector<double> > ts_m (np, ts0);
  vector<vector<double> > y0_m (np, y0);
  vector<vector<var> > theta_var_m (np, stan::math::to_var(theta));

  vector<vector<double> > x_r_m (np, x_r);
  vector<vector<int> > x_i_m (np, x_i);

  vector<vector<var> > y = stan::math::integrate_ode_bdf(f, y0, t0, ts0, theta_var, x_r, x_i);
  vector<Eigen::Matrix<var, -1, -1> > y_m = pmx_integrate_ode_bdf(f, y0_m, t0, ts_m, theta_var_m , x_r_m, x_i_m);

  for (int i = 0; i < np; ++i) {
    torsten::test::test_grad(theta_var, theta_var_m[i], y, y_m[i], 1e-7, 1e-6);
  }
}

TEST_F(TorstenOdeTest_lorenz, fwd_sensitivity_theta_AD_adams_mpi_performance) {
  using torsten::dsolve::PMXCvodesFwdSystem;
  using torsten::dsolve::pmx_integrate_ode_adams;
  using stan::math::var;
  using std::vector;

  // size of population
  const int np = 10;
  std::vector<double> ts0 {ts};

  vector<var> theta_var = stan::math::to_var(theta);

  vector<vector<double> > ts_m (np, ts0);
  vector<vector<double> > y0_m (np, y0);
  vector<vector<var> > theta_var_m (np, stan::math::to_var(theta));

  vector<vector<double> > x_r_m (np, x_r);
  vector<vector<int> > x_i_m (np, x_i);

  vector<vector<var> > y = stan::math::integrate_ode_adams(f, y0, t0, ts0, theta_var, x_r, x_i);
  vector<Eigen::Matrix<var, -1, -1> > y_m = pmx_integrate_ode_adams(f, y0_m, t0, ts_m, theta_var_m , x_r_m, x_i_m);

  for (int i = 0; i < np; ++i) {
    torsten::test::test_grad(theta_var, theta_var_m[i], y, y_m[i], 1e-6, 5e-5);
  }
}

TEST_F(TorstenOdeTest_neutropenia, fwd_sensitivity_theta_AD_bdf_mpi) {
  using torsten::dsolve::PMXCvodesFwdSystem;
  using torsten::dsolve::pmx_integrate_ode_bdf;
  using stan::math::var;
  using std::vector;

  torsten::mpi::Envionment::init();

  // size of population
  const int np = 10;

  vector<var> theta_var = stan::math::to_var(theta);

  vector<vector<double> > ts_m (np, ts);
  vector<vector<double> > y0_m (np, y0);
  vector<vector<var> > theta_var_m (np, stan::math::to_var(theta));
  vector<vector<double> > x_r_m (np, x_r);
  vector<vector<int> > x_i_m (np, x_i);

  vector<vector<var> > y = stan::math::integrate_ode_bdf(f, y0, t0, ts, theta_var, x_r, x_i);
  vector<Eigen::Matrix<var, -1, -1> > y_m = pmx_integrate_ode_bdf(f, y0_m, t0, ts_m, theta_var_m , x_r_m, x_i_m);

  for (int i = 0; i < np; ++i) {
    torsten::test::test_grad(theta_var, theta_var_m[i], y, y_m[i], 1e-7, 1e-7);    
  }
}

TEST_F(TorstenOdeTest_neutropenia, fwd_sensitivity_theta_AD_adams_mpi) {
  using torsten::dsolve::PMXCvodesFwdSystem;
  using torsten::dsolve::pmx_integrate_ode_adams;
  using stan::math::var;
  using std::vector;

  torsten::mpi::Envionment::init();

  // size of population
  const int np = 10;

  vector<var> theta_var = stan::math::to_var(theta);

  vector<vector<double> > ts_m (np, ts);
  vector<vector<double> > y0_m (np, y0);
  vector<vector<var> > theta_var_m (np, stan::math::to_var(theta));
  vector<vector<double> > x_r_m (np, x_r);
  vector<vector<int> > x_i_m (np, x_i);

  vector<vector<var> > y = stan::math::integrate_ode_adams(f, y0, t0, ts, theta_var, x_r, x_i);
  vector<Eigen::Matrix<var, -1, -1> > y_m = pmx_integrate_ode_adams(f, y0_m, t0, ts_m, theta_var_m , x_r_m, x_i_m);

  for (int i = 0; i < np; ++i) {
    torsten::test::test_grad(theta_var, theta_var_m[i], y, y_m[i], 1e-7, 5e-6);    
  }
}

TEST_F(TorstenOdeTest_neutropenia, mpi_rank_exception_data_only) {
  using torsten::dsolve::PMXCvodesFwdSystem;
  using stan::math::var;
  using std::vector;

  // size of population
  const int np = 5;

  vector<vector<double> > ts_m (np, ts);
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
  EXPECT_THROW(torsten::dsolve::pmx_integrate_ode_bdf(f, y0_m, t0, ts_m, theta_m , x_r_m, x_i_m),
               std::exception);
  EXPECT_THROW(torsten::dsolve::pmx_integrate_ode_adams(f, y0_m, t0, ts_m, theta_m , x_r_m, x_i_m),
               std::exception);
#ifdef TORSTEN_MPI
  MPI_Barrier(comm);
#endif

  theta_m[id][0] = -1.0E+30;
  theta_m[id][1] = -1.0E+30;
  theta_m[id][4] = -1.0E+30;
  EXPECT_THROW(torsten::dsolve::pmx_integrate_ode_bdf(f, y0_m, t0, ts_m, theta_m , x_r_m, x_i_m,
                                                     0, 1e-6, 1e-6, 1e2),
               std::exception);
  EXPECT_THROW(torsten::dsolve::pmx_integrate_ode_adams(f, y0_m, t0, ts_m, theta_m , x_r_m, x_i_m,
                                                       0, 1e-6, 1e-6, 1e2),
               std::exception);
#ifdef TORSTEN_MPI
  MPI_Barrier(comm);
#endif

  theta_m[id][0] = std::numeric_limits<double>::quiet_NaN();
  EXPECT_THROW(torsten::dsolve::pmx_integrate_ode_bdf(f, y0_m, t0, ts_m, theta_m , x_r_m, x_i_m),
               std::exception);
  EXPECT_THROW(torsten::dsolve::pmx_integrate_ode_bdf(f, y0_m, t0, ts_m, theta_m , x_r_m, x_i_m),
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
  EXPECT_THROW(torsten::dsolve::pmx_integrate_ode_bdf(f, y0_m, t0, ts_m, theta_m , x_r_m, x_i_m),
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
  EXPECT_THROW(torsten::dsolve::pmx_integrate_ode_bdf(f, y0_m, t0, ts_m, theta_m , x_r_m, x_i_m,
                                                     0, 1e-6, 1e-6, 1e2),
               std::exception);
#ifdef TORSTEN_MPI
  MPI_Barrier(comm);
#endif
}

TEST_F(TorstenOdeTest_neutropenia, mpi_rank_exception_par_var) {
  using torsten::dsolve::PMXCvodesFwdSystem;
  using stan::math::var;
  using std::vector;

  // size of population
  const int np = 5;

  vector<var> theta_var(stan::math::to_var(theta));
  vector<vector<double> > ts_m (np, ts);
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
  EXPECT_THROW(torsten::dsolve::pmx_integrate_ode_bdf(f, y0_m, t0, ts_m, theta_m_v , x_r_m, x_i_m),
               std::exception);
  EXPECT_THROW(torsten::dsolve::pmx_integrate_ode_adams(f, y0_m, t0, ts_m, theta_m_v , x_r_m, x_i_m),
               std::exception);
#ifdef TORSTEN_MPI
  MPI_Barrier(comm);
#endif

  theta_m_v[id][0] = -1.0E+30;
  theta_m_v[id][1] = -1.0E+30;
  theta_m_v[id][4] = -1.0E+30;
  EXPECT_THROW(torsten::dsolve::pmx_integrate_ode_bdf(f, y0_m, t0, ts_m, theta_m_v , x_r_m, x_i_m,
                                                     0, 1e-6, 1e-6, 1e2),
               std::exception);
  EXPECT_THROW(torsten::dsolve::pmx_integrate_ode_adams(f, y0_m, t0, ts_m, theta_m_v , x_r_m, x_i_m,
                                                       0, 1e-6, 1e-6, 1e2),
               std::exception);
#ifdef TORSTEN_MPI
  MPI_Barrier(comm);
#endif

  theta_m_v[id][0] = std::numeric_limits<double>::quiet_NaN();
  EXPECT_THROW(torsten::dsolve::pmx_integrate_ode_bdf(f, y0_m, t0, ts_m, theta_m_v , x_r_m, x_i_m),
               std::exception);
  EXPECT_THROW(torsten::dsolve::pmx_integrate_ode_bdf(f, y0_m, t0, ts_m, theta_m_v , x_r_m, x_i_m),
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
  EXPECT_THROW(torsten::dsolve::pmx_integrate_ode_bdf(f, y0_m, t0, ts_m, theta_m_v , x_r_m, x_i_m),
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
  EXPECT_THROW(torsten::dsolve::pmx_integrate_ode_bdf(f, y0_m, t0, ts_m, theta_m_v , x_r_m, x_i_m,
                                                     0, 1e-6, 1e-6, 1e2),
               std::exception);
#ifdef TORSTEN_MPI
  MPI_Barrier(comm);
#endif
}
