#include <stan/math/rev.hpp>
#include <gtest/gtest.h>
#include <stan/math/torsten/dsolve/pmx_integrate_ode_adams.hpp>
#include <stan/math/torsten/dsolve/pmx_integrate_ode_group_adams.hpp>
#include <stan/math.hpp>
#include <stan/math/rev/core.hpp>
#include <stan/math/torsten/mpi.hpp>
#include <test/unit/math/rev/fun/util.hpp>
#include <stan/math/torsten/test/unit/pmx_ode_test_fixture.hpp>
#include <test/unit/math/prim/functor/harmonic_oscillator.hpp>
#include <stan/math/rev/functor/integrate_ode_bdf.hpp>
#include <nvector/nvector_serial.h>
#include <test/unit/util.hpp>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <iostream>
#include <fstream>
#include <vector>
#include <limits>

using stan::math::integrate_ode_adams;
using torsten::pmx_integrate_ode_adams;
using torsten::dsolve::PMXVariadicOdeSystem;
using torsten::dsolve::PMXOdeintIntegrator;
using torsten::pmx_integrate_ode_group_adams;
using dsolve::OdeObserver;
using dsolve::OdeDataObserver;
using stan::math::var;
using std::vector;

#if defined(STAN_LANG_MPI) || defined(TORSTEN_MPI)
TORSTEN_MPI_SESSION_INIT;
#endif

TEST_F(TorstenOdeTest_sho, cvodes_adams_ivp_system_matrix_result) {
  std::vector<std::vector<double> > y1(integrate_ode_adams(f, y0, t0, ts, theta, x_r, x_i, msgs, atol, rtol, max_num_steps));

  using Ode = PMXVariadicOdeSystem<harm_osc_ode_fun_eigen, double, double, std::vector<double>,std::vector<double>,std::vector<int> >;
  Ode ode{f_eigen, t0, ts, y0_vec, msgs, theta, x_r, x_i};
  using scheme_t = torsten::dsolve::odeint_scheme_rk45;
  PMXOdeintIntegrator<scheme_t> solver(rtol, atol, max_num_steps);
  OdeDataObserver<Ode> observer(ode);
  solver.integrate(ode, observer);

  torsten::test::test_val(y1, observer.y);
}

TEST_F(TorstenOdeTest_lorenz, cvodes_adams_ivp_system) {
  std::vector<std::vector<double> > y1(integrate_ode_adams(f, y0, t0, ts, theta , x_r, x_i));
  std::vector<std::vector<double> > y2(pmx_integrate_ode_adams(f, y0, t0, ts, theta , x_r, x_i));
  torsten::test::test_val(y1, y2);
}

TEST_F(TorstenOdeTest_lorenz, cvodes_adams_ivp_system_matrix_result) {
 std::vector<std::vector<double> > y1(integrate_ode_adams(f, y0, t0, ts, theta, x_r, x_i, msgs, atol, rtol, max_num_steps));

  using Ode = PMXVariadicOdeSystem<lorenz_ode_eigen_fun, double, double, std::vector<double>,std::vector<double>,std::vector<int> >;
  Ode ode{f_eigen, t0, ts, y0_vec, msgs, theta, x_r, x_i};
  using scheme_t = torsten::dsolve::odeint_scheme_rk45;
  PMXOdeintIntegrator<scheme_t> solver(rtol, atol, max_num_steps);
  OdeDataObserver<Ode> observer(ode);
  solver.integrate(ode, observer);

  torsten::test::test_val(y1, observer.y);
}

TEST_F(TorstenOdeTest_chem, cvodes_adams_fwd_sensitivity_theta) {
  std::vector<var> theta_var1 = stan::math::to_var(theta);
  std::vector<var> theta_var2 = stan::math::to_var(theta);

  std::vector<std::vector<stan::math::var>> y1 = integrate_ode_adams(f, y0, t0, ts, theta_var1, x_r, x_i);
  std::vector<std::vector<stan::math::var>> y2 = pmx_integrate_ode_adams(f, y0, t0, ts, theta_var2, x_r, x_i);
  torsten::test::test_grad(theta_var1, theta_var2, y1, y2, 1.E-6, 3.E-6);
}

TEST_F(TorstenOdeTest_sho, fwd_sensitivity_ts) {
  using stan::math::value_of;
  std::vector<var> ts_var = stan::math::to_var(ts);

  std::vector<std::vector<var> > y, y1, y2;
  y = pmx_integrate_ode_adams(f, y0, t0, ts_var, theta, x_r, x_i, msgs);

  std::vector<double> g(ts.size()), fval(y0.size());
  for (size_t i = 0; i < ts.size(); ++i) {
    fval = f(value_of(ts[i]), value_of(y[i]), theta, x_r, x_i, msgs);
    for (size_t j = 0; j < y0.size(); ++j) {
      stan::math::set_zero_all_adjoints();
      y[i][j].grad(ts_var, g);
      for (size_t k = 0; k < ts.size(); ++k) {
        if (k == i) {
          EXPECT_FLOAT_EQ(g[k], fval[j]);
        } else {
          EXPECT_FLOAT_EQ(g[k], 0.0);
        }
      }
    }
  }
}

TEST_F(TorstenOdeTest_sho, adams_theta_var_matrix_result) {
  std::vector<var> theta_var = stan::math::to_var(theta);
  vector<vector<var> > y1(integrate_ode_adams(f, y0, t0, ts, theta_var, x_r, x_i, msgs, atol, rtol, max_num_steps));

  using Ode = PMXVariadicOdeSystem<harm_osc_ode_fun_eigen, double, double, std::vector<var>,std::vector<double>,std::vector<int> >;
  Ode ode{f_eigen, t0, ts,y0_vec, msgs, theta_var, x_r, x_i};
  using scheme_t = torsten::dsolve::odeint_scheme_rk45;
  PMXOdeintIntegrator<scheme_t> solver(rtol, atol, max_num_steps);
  OdeDataObserver<Ode> observer(ode);
  solver.integrate(ode, observer);

  torsten::test::test_grad(theta_var, y1, observer.y, 1.e-8, 1.e-8);
}

TEST_F(TorstenOdeTest_lorenz, cvodes_adams_fwd_sensitivity_theta) {
  std::vector<var> theta_var1 = stan::math::to_var(theta);
  std::vector<var> theta_var2 = stan::math::to_var(theta);

  std::vector<std::vector<stan::math::var>> y1 = integrate_ode_adams(f, y0, t0, ts, theta_var1, x_r, x_i);
  std::vector<std::vector<stan::math::var>> y2 = pmx_integrate_ode_adams(f, y0, t0, ts, theta_var2, x_r, x_i);
  torsten::test::test_grad(theta_var1, theta_var2, y1, y2, 1.E-6, 1.E-5);
}

TEST_F(TorstenOdeTest_sho, cvodes_adams_fwd_sensitivity_theta) {
  std::vector<var> theta_var1 = stan::math::to_var(theta);
  std::vector<var> theta_var2 = stan::math::to_var(theta);

  std::vector<std::vector<stan::math::var>> y1 = integrate_ode_adams(f, y0, t0, ts, theta_var1, x_r, x_i);
  std::vector<std::vector<stan::math::var>> y2 = pmx_integrate_ode_adams(f, y0, t0, ts, theta_var2, x_r, x_i);
  torsten::test::test_grad(theta_var1, theta_var2, y1, y2, 1.E-6, 6.E-6);
}

TEST_F(TorstenOdeTest_chem, cvodes_adams_fwd_sensitivity_y0) {
  std::vector<var> y0_var1 = stan::math::to_var(y0);
  std::vector<var> y0_var2 = stan::math::to_var(y0);

  std::vector<std::vector<stan::math::var>> y1 = integrate_ode_adams(f, y0_var1, t0, ts, theta, x_r, x_i);
  std::vector<std::vector<stan::math::var>> y2 = pmx_integrate_ode_adams(f, y0_var2, t0, ts, theta, x_r, x_i);
  torsten::test::test_grad(y0_var1, y0_var2, y1, y2, 5.E-6, 5.E-6);
}

TEST_F(TorstenOdeTest_lorenz, cvodes_adams_fwd_sensitivity_y0) {
  std::vector<var> y0_var1 = stan::math::to_var(y0);
  std::vector<var> y0_var2 = stan::math::to_var(y0);

  std::vector<std::vector<stan::math::var>> y1 = integrate_ode_adams(f, y0_var1, t0, ts, theta, x_r, x_i);
  std::vector<std::vector<stan::math::var>> y2 = pmx_integrate_ode_adams(f, y0_var2, t0, ts, theta, x_r, x_i);
  torsten::test::test_grad(y0_var1, y0_var2, y1, y2, 1.E-6, 1.E-6);
}

TEST_F(TorstenOdeTest_sho, cvodes_adams_fwd_sensitivity_y0) {
  std::vector<var> y0_var1 = stan::math::to_var(y0);
  std::vector<var> y0_var2 = stan::math::to_var(y0);

  std::vector<std::vector<stan::math::var>> y1 = integrate_ode_adams(f, y0_var1, t0, ts, theta, x_r, x_i);
  std::vector<std::vector<stan::math::var>> y2 = pmx_integrate_ode_adams(f, y0_var2, t0, ts, theta, x_r, x_i);
  torsten::test::test_grad(y0_var1, y0_var2, y1, y2, 5.E-6, 5.E-6);
}

TEST_F(TorstenOdeTest_sho, fwd_sensitivity_theta_y0) {
  std::vector<var> theta_var1 = stan::math::to_var(theta);
  std::vector<var> y0_var1 = stan::math::to_var(y0);
  std::vector<var> theta_var2 = stan::math::to_var(theta);
  std::vector<var> y0_var2 = stan::math::to_var(y0);

  std::vector<std::vector<stan::math::var>> y1 = integrate_ode_adams(f, y0_var1, t0, ts, theta_var1, x_r, x_i);
  std::vector<std::vector<stan::math::var>> y2 = pmx_integrate_ode_adams(f, y0_var2, t0, ts, theta_var2, x_r, x_i);
  torsten::test::test_grad(y0_var1, y0_var2, y1, y2, 1.E-6, 1.E-6);
  torsten::test::test_grad(theta_var1, theta_var2, y1, y2, 1.E-6, 6.E-6);
}

TEST_F(TorstenOdeTest_lorenz, fwd_sensitivity_theta_y0) {
  std::vector<var> theta_var1 = stan::math::to_var(theta);
  std::vector<var> y0_var1 = stan::math::to_var(y0);
  std::vector<var> theta_var2 = stan::math::to_var(theta);
  std::vector<var> y0_var2 = stan::math::to_var(y0);

  std::vector<std::vector<stan::math::var>> y1 = integrate_ode_adams(f, y0_var1, t0, ts, theta_var1, x_r, x_i);
  std::vector<std::vector<stan::math::var>> y2 = pmx_integrate_ode_adams(f, y0_var2, t0, ts, theta_var2, x_r, x_i);
  torsten::test::test_grad(y0_var1, y0_var2, y1, y2, 1.E-6, 1.E-5);
  torsten::test::test_grad(theta_var1, theta_var2, y1, y2, 1.E-6, 1.E-5);
}

TEST_F(TorstenOdeTest_chem, fwd_sensitivity_theta_y0) {
  std::vector<var> theta_var1 = stan::math::to_var(theta);
  std::vector<var> y0_var1 = stan::math::to_var(y0);
  std::vector<var> theta_var2 = stan::math::to_var(theta);
  std::vector<var> y0_var2 = stan::math::to_var(y0);

  std::vector<std::vector<stan::math::var>> y1 = integrate_ode_adams(f, y0_var1, t0, ts, theta_var1, x_r, x_i);
  std::vector<std::vector<stan::math::var>> y2 = pmx_integrate_ode_adams(f, y0_var2, t0, ts, theta_var2, x_r, x_i);
  torsten::test::test_grad(y0_var1, y0_var2, y1, y2, 1.E-8, 6.E-6);
  torsten::test::test_grad(theta_var1, theta_var2, y1, y2, 1.E-8, 6.E-6);
}

// test sequential run of group solver
TEST_F(TorstenOdeTest_neutropenia, group_adams_fwd_sensitivity_theta) {
  // size of population
  const int np = 2;

  vector<var> theta_var = stan::math::to_var(theta);

  vector<int> len(np, ts.size());
  vector<double> ts_m;
  ts_m.reserve(np * ts.size());
  for (int i = 0; i < np; ++i) ts_m.insert(ts_m.end(), ts.begin(), ts.end());

  vector<vector<double> > y0_m (np, y0);
  vector<vector<var> > theta_var_m (np, stan::math::to_var(theta));
  vector<vector<double> > x_r_m (np, x_r);
  vector<vector<int> > x_i_m (np, x_i);

  vector<vector<var> > y = integrate_ode_adams(f, y0, t0, ts, theta_var, x_r, x_i);
  Eigen::Matrix<var, -1, -1> y_m1 = pmx_integrate_ode_group_adams(f, y0_m, t0, len, ts_m, theta_var_m , x_r_m, x_i_m);
  Eigen::Matrix<var, -1, -1> y_m2 = pmx_integrate_ode_group_adams(f, y0_m, t0, len, ts_m, theta_var_m , x_r_m, x_i_m);

  EXPECT_EQ(y_m1.cols(), ts_m.size());
  EXPECT_EQ(y_m2.cols(), ts_m.size());
  int icol = 0;
  for (int i = 0; i < np; ++i) {
    stan::math::matrix_v y_i = y_m1.block(0, icol, y0.size(), len[i]);
    torsten::test::test_grad(theta_var, theta_var_m[i], y, y_i, 1e-8, 2e-8);
    icol += len[i];
  }

  icol = 0;
  for (int i = 0; i < np; ++i) {
    stan::math::matrix_v y_i = y_m2.block(0, icol, y0.size(), len[i]);
    torsten::test::test_grad(theta_var, theta_var_m[i], y, y_i, 1e-8, 2e-8);
    icol += len[i];
  }
}
