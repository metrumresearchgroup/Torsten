#ifdef TORSTEN_MPI

#include <stan/math.hpp>
#include <stan/math/rev/core.hpp>
#include <test/unit/math/rev/fun/util.hpp>
#include <stan/math/torsten/dsolve/pmx_cvodes_integrator.hpp>
#include <stan/math/prim/fun/eval.hpp>
#include <stan/math/rev/functor/jacobian.hpp>
#include <stan/math/torsten/mpi.hpp>
#include <stan/math/torsten/dsolve/pmx_integrate_ode_adams.hpp>
#include <stan/math/torsten/dsolve/pmx_integrate_ode_bdf.hpp>
#include <stan/math/torsten/dsolve/pmx_integrate_ode_group_adams.hpp>
#include <stan/math/torsten/dsolve/pmx_integrate_ode_group_bdf.hpp>
#include <stan/math/torsten/test/unit/pmx_ode_test_fixture.hpp>
#include <test/unit/math/prim/functor/harmonic_oscillator.hpp>
#include <stan/math/rev/functor/integrate_ode_bdf.hpp>
#include <nvector/nvector_serial.h>
#include <boost/mpi.hpp>
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

TORSTEN_MPI_SESSION_INIT;

TEST_F(TorstenOdeTest_chem, fwd_sensitivity_theta_bdf_mpi) {
  
  using torsten::pmx_integrate_ode_group_bdf;
  using stan::math::var;
  using std::vector;

  stan::math::mpi::Envionment::init();

  // size of population
  const int np = 100;
  std::vector<double> ts0 {ts};

  vector<var> theta_var = stan::math::to_var(theta);

  vector<int> len(np, ts0.size());
  vector<double> ts_m;
  ts_m.reserve(np * ts0.size());
  for (int i = 0; i < np; ++i) ts_m.insert(ts_m.end(), ts0.begin(), ts0.end());

  vector<vector<double> > y0_m (np, y0);
  vector<vector<var> > theta_var_m (np, theta_var);
  vector<vector<double> > x_r_m (np, x_r);
  vector<vector<int> > x_i_m (np, x_i);

  Eigen::Matrix<var, -1, -1> y_m = pmx_integrate_ode_group_bdf(f, y0_m, t0, len, ts_m, theta_var_m , x_r_m, x_i_m);
}

TEST_F(TorstenOdeTest_lorenz, fwd_sensitivity_theta_bdf_mpi) {
  
  using torsten::pmx_integrate_ode_group_bdf;
  using stan::math::var;
  using std::vector;

  stan::math::mpi::Envionment::init();

  // size of population
  const int np = 100;
  std::vector<double> ts0 {ts};

  vector<var> theta_var = stan::math::to_var(theta);

  vector<int> len(np, ts0.size());
  vector<double> ts_m;
  ts_m.reserve(np * ts0.size());
  for (int i = 0; i < np; ++i) ts_m.insert(ts_m.end(), ts0.begin(), ts0.end());

  vector<vector<double> > y0_m (np, y0);
  vector<vector<var> > theta_var_m (np, theta_var);
  vector<vector<double> > x_r_m (np, x_r);
  vector<vector<int> > x_i_m (np, x_i);

  Eigen::Matrix<var, -1, -1> y_m = pmx_integrate_ode_group_bdf(f, y0_m, t0, len, ts_m, theta_var_m , x_r_m, x_i_m);
}

TEST_F(TorstenOdeTest_neutropenia, fwd_sensitivity_theta_bdf_mpi) {
  
  using torsten::pmx_integrate_ode_group_bdf;
  using stan::math::var;
  using std::vector;

  stan::math::mpi::Envionment::init();

  // size of population
  const int np = 100;
  std::vector<double> ts0 {ts};

  vector<var> theta_var = stan::math::to_var(theta);

  vector<int> len(np, ts0.size());
  vector<double> ts_m;
  ts_m.reserve(np * ts0.size());
  for (int i = 0; i < np; ++i) ts_m.insert(ts_m.end(), ts0.begin(), ts0.end());

  vector<vector<double> > y0_m (np, y0);
  vector<vector<var> > theta_var_m (np, theta_var);
  vector<vector<double> > x_r_m (np, x_r);
  vector<vector<int> > x_i_m (np, x_i);

  Eigen::Matrix<var, -1, -1> y_m = pmx_integrate_ode_group_bdf(f, y0_m, t0, len, ts_m, theta_var_m , x_r_m, x_i_m);
}

TEST_F(TorstenOdeTest_chem, fwd_sensitivity_theta_adams_mpi) {
  
  using torsten::pmx_integrate_ode_group_adams;
  using stan::math::var;
  using std::vector;

  stan::math::mpi::Envionment::init();

  // size of population
  const int np = 100;
  std::vector<double> ts0 {ts};

  vector<var> theta_var = stan::math::to_var(theta);

  vector<int> len(np, ts0.size());
  vector<double> ts_m;
  ts_m.reserve(np * ts0.size());
  for (int i = 0; i < np; ++i) ts_m.insert(ts_m.end(), ts0.begin(), ts0.end());

  vector<vector<double> > y0_m (np, y0);
  vector<vector<var> > theta_var_m (np, theta_var);
  vector<vector<double> > x_r_m (np, x_r);
  vector<vector<int> > x_i_m (np, x_i);

  Eigen::Matrix<var, -1, -1> y_m = pmx_integrate_ode_group_adams(f, y0_m, t0, len, ts_m, theta_var_m , x_r_m, x_i_m);
}

TEST_F(TorstenOdeTest_lorenz, fwd_sensitivity_theta_adams_mpi) {
  
  using torsten::pmx_integrate_ode_group_adams;
  using stan::math::var;
  using std::vector;

  stan::math::mpi::Envionment::init();

  // size of population
  const int np = 100;
  std::vector<double> ts0 {ts};

  vector<var> theta_var = stan::math::to_var(theta);

  vector<int> len(np, ts0.size());
  vector<double> ts_m;
  ts_m.reserve(np * ts0.size());
  for (int i = 0; i < np; ++i) ts_m.insert(ts_m.end(), ts0.begin(), ts0.end());

  vector<vector<double> > y0_m (np, y0);
  vector<vector<var> > theta_var_m (np, theta_var);
  vector<vector<double> > x_r_m (np, x_r);
  vector<vector<int> > x_i_m (np, x_i);

  Eigen::Matrix<var, -1, -1> y_m = pmx_integrate_ode_group_adams(f, y0_m, t0, len, ts_m, theta_var_m , x_r_m, x_i_m);
}

TEST_F(TorstenOdeTest_neutropenia, fwd_sensitivity_theta_adams_mpi) {
  
  using torsten::pmx_integrate_ode_group_adams;
  using stan::math::var;
  using std::vector;

  stan::math::mpi::Envionment::init();

  // size of population
  const int np = 100;
  std::vector<double> ts0 {ts};

  vector<var> theta_var = stan::math::to_var(theta);

  vector<int> len(np, ts0.size());
  vector<double> ts_m;
  ts_m.reserve(np * ts0.size());
  for (int i = 0; i < np; ++i) ts_m.insert(ts_m.end(), ts0.begin(), ts0.end());

  vector<vector<double> > y0_m (np, y0);
  vector<vector<var> > theta_var_m (np, theta_var);
  vector<vector<double> > x_r_m (np, x_r);
  vector<vector<int> > x_i_m (np, x_i);

  Eigen::Matrix<var, -1, -1> y_m = pmx_integrate_ode_group_adams(f, y0_m, t0, len, ts_m, theta_var_m , x_r_m, x_i_m);
}

#endif
