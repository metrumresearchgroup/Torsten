#include <stan/math/rev/arr.hpp>
#include <gtest/gtest.h>
#include <stan/math/torsten/dsolve/pmx_integrate_ode_rk45.hpp>
#include <stan/math.hpp>
#include <stan/math/rev/core.hpp>
#include <test/unit/math/rev/mat/fun/util.hpp>
#include <test/unit/math/torsten/pmx_ode_test_fixture.hpp>
#include <test/unit/math/prim/arr/functor/harmonic_oscillator.hpp>
#include <stan/math/rev/mat/functor/integrate_ode_bdf.hpp>
#include <nvector/nvector_serial.h>
#include <test/unit/util.hpp>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <iostream>
#include <fstream>
#include <vector>
#include <limits>

using stan::math::integrate_ode_rk45;
using torsten::pmx_integrate_ode_rk45;
using stan::math::var;
using std::vector;

TEST_F(TorstenOdeTest_sho, odeint_rk45_ivp_system) {
  std::vector<std::vector<double> > y1(integrate_ode_rk45(f, y0, t0, ts, theta , x_r, x_i));
  std::vector<std::vector<double> > y2(pmx_integrate_ode_rk45(f, y0, t0, ts, theta , x_r, x_i));
  torsten::test::test_val(y1, y2);
}

TEST_F(TorstenOdeTest_sho, odeint_rk45_ivp_system_matrix_result) {
  std::vector<std::vector<double> > y1(integrate_ode_rk45(f, y0, t0, ts, theta, x_r, x_i, msgs, atol, rtol, max_num_steps));

  using Ode = dsolve::PMXOdeintSystem<harm_osc_ode_fun, double, double, double>;
  dsolve::PMXOdeService<Ode, dsolve::Odeint> serv(y0.size(), theta.size());
  Ode ode{serv, f, t0, ts, y0, theta, x_r, x_i, msgs};
  using scheme_t = boost::numeric::odeint::runge_kutta_dopri5<std::vector<double>, double, std::vector<double>, double>;
  dsolve::PMXOdeintIntegrator<scheme_t> solver(rtol, atol, max_num_steps);
  Eigen::MatrixXd y2 = solver.integrate<Ode, false>(ode);

  torsten::test::test_val(y1, y2);
}

TEST_F(TorstenOdeTest_lorenz, odeint_rk45_ivp_system) {
  std::vector<std::vector<double> > y1(integrate_ode_rk45(f, y0, t0, ts, theta , x_r, x_i));
  std::vector<std::vector<double> > y2(pmx_integrate_ode_rk45(f, y0, t0, ts, theta , x_r, x_i));
  torsten::test::test_val(y1, y2);
}

TEST_F(TorstenOdeTest_lorenz, odeint_rk45_ivp_system_matrix_result) {
 std::vector<std::vector<double> > y1(integrate_ode_rk45(f, y0, t0, ts, theta, x_r, x_i, msgs, atol, rtol, max_num_steps));

  using Ode = dsolve::PMXOdeintSystem<lorenz_ode_fun, double, double, double>;
  dsolve::PMXOdeService<Ode, dsolve::Odeint> serv(y0.size(), theta.size());
  Ode ode{serv, f, t0, ts, y0, theta, x_r, x_i, msgs};
  using scheme_t = boost::numeric::odeint::runge_kutta_dopri5<std::vector<double>, double, std::vector<double>, double>;
  dsolve::PMXOdeintIntegrator<scheme_t> solver(rtol, atol, max_num_steps);
  Eigen::MatrixXd y2 = solver.integrate<Ode, false>(ode);

  torsten::test::test_val(y1, y2);
}

TEST_F(TorstenOdeTest_chem, odeint_rk45_ivp_system) {
  std::vector<std::vector<double> > y1(integrate_ode_rk45(f, y0, t0, ts, theta , x_r, x_i));
  std::vector<std::vector<double> > y2(pmx_integrate_ode_rk45(f, y0, t0, ts, theta , x_r, x_i));
  torsten::test::test_val(y1, y2);
}

TEST_F(TorstenOdeTest_chem, odeint_rk45_fwd_sensitivity_theta) {
  std::vector<var> theta_var1 = stan::math::to_var(theta);
  std::vector<var> theta_var2 = stan::math::to_var(theta);

  std::vector<std::vector<stan::math::var>> y1 = integrate_ode_rk45(f, y0, t0, ts, theta_var1, x_r, x_i);
  std::vector<std::vector<stan::math::var>> y2 = pmx_integrate_ode_rk45(f, y0, t0, ts, theta_var2, x_r, x_i);
  torsten::test::test_grad(theta_var1, theta_var2, y1, y2, 1.E-10, 1.E-6);
}

TEST_F(TorstenOdeTest_sho, fwd_sensitivity_ts) {
  using stan::math::value_of;
  std::vector<var> ts_var = stan::math::to_var(ts);

  std::vector<std::vector<var> > y, y1, y2;
  y = pmx_integrate_ode_rk45(f, y0, t0, ts_var, theta, x_r, x_i, msgs);

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

TEST_F(TorstenOdeTest_sho, rk45_theta_var_matrix_result) {
  std::vector<var> theta_var = stan::math::to_var(theta);
  vector<vector<var> > y1(integrate_ode_rk45(f, y0, t0, ts, theta_var, x_r, x_i, msgs, atol, rtol, max_num_steps));

  using Ode = dsolve::PMXOdeintSystem<harm_osc_ode_fun, double, double, var>;
  dsolve::PMXOdeService<Ode, dsolve::Odeint> serv(y0.size(), theta.size());
  Ode ode{serv, f, t0, ts, y0, theta_var, x_r, x_i, msgs};
  using scheme_t = boost::numeric::odeint::runge_kutta_dopri5<std::vector<double>, double, std::vector<double>, double>;
  dsolve::PMXOdeintIntegrator<scheme_t> solver(rtol, atol, max_num_steps);
  Eigen::MatrixXd y2 = solver.integrate<Ode, false>(ode);

  torsten::test::test_grad(theta_var, y1, y2, 1.e-8, 1.e-8);
}

TEST_F(TorstenOdeTest_lorenz, odeint_rk45_fwd_sensitivity_theta) {
  std::vector<var> theta_var1 = stan::math::to_var(theta);
  std::vector<var> theta_var2 = stan::math::to_var(theta);

  std::vector<std::vector<stan::math::var>> y1 = integrate_ode_rk45(f, y0, t0, ts, theta_var1, x_r, x_i);
  std::vector<std::vector<stan::math::var>> y2 = pmx_integrate_ode_rk45(f, y0, t0, ts, theta_var2, x_r, x_i);
  torsten::test::test_grad(theta_var1, theta_var2, y1, y2, 1.E-10, 1.E-8);
}

TEST_F(TorstenOdeTest_sho, odeint_rk45_fwd_sensitivity_theta) {
  std::vector<var> theta_var1 = stan::math::to_var(theta);
  std::vector<var> theta_var2 = stan::math::to_var(theta);

  std::vector<std::vector<stan::math::var>> y1 = integrate_ode_rk45(f, y0, t0, ts, theta_var1, x_r, x_i);
  std::vector<std::vector<stan::math::var>> y2 = pmx_integrate_ode_rk45(f, y0, t0, ts, theta_var2, x_r, x_i);
  torsten::test::test_grad(theta_var1, theta_var2, y1, y2, 1.E-10, 1.E-8);
}

TEST_F(TorstenOdeTest_chem, odeint_rk45_fwd_sensitivity_y0) {
  std::vector<var> y0_var1 = stan::math::to_var(y0);
  std::vector<var> y0_var2 = stan::math::to_var(y0);

  std::vector<std::vector<stan::math::var>> y1 = integrate_ode_rk45(f, y0_var1, t0, ts, theta, x_r, x_i);
  std::vector<std::vector<stan::math::var>> y2 = pmx_integrate_ode_rk45(f, y0_var2, t0, ts, theta, x_r, x_i);
  torsten::test::test_grad(y0_var1, y0_var2, y1, y2, 1.E-10, 1.E-8);
}

TEST_F(TorstenOdeTest_lorenz, odeint_rk45_fwd_sensitivity_y0) {
  std::vector<var> y0_var1 = stan::math::to_var(y0);
  std::vector<var> y0_var2 = stan::math::to_var(y0);

  std::vector<std::vector<stan::math::var>> y1 = integrate_ode_rk45(f, y0_var1, t0, ts, theta, x_r, x_i);
  std::vector<std::vector<stan::math::var>> y2 = pmx_integrate_ode_rk45(f, y0_var2, t0, ts, theta, x_r, x_i);
  torsten::test::test_grad(y0_var1, y0_var2, y1, y2, 1.E-10, 1.E-8);
}

TEST_F(TorstenOdeTest_sho, odeint_rk45_fwd_sensitivity_y0) {
  std::vector<var> y0_var1 = stan::math::to_var(y0);
  std::vector<var> y0_var2 = stan::math::to_var(y0);

  std::vector<std::vector<stan::math::var>> y1 = integrate_ode_rk45(f, y0_var1, t0, ts, theta, x_r, x_i);
  std::vector<std::vector<stan::math::var>> y2 = pmx_integrate_ode_rk45(f, y0_var2, t0, ts, theta, x_r, x_i);
  torsten::test::test_grad(y0_var1, y0_var2, y1, y2, 1.E-10, 1.E-8);
}

TEST_F(TorstenOdeTest_sho, fwd_sensitivity_theta_y0) {
  std::vector<var> theta_var1 = stan::math::to_var(theta);
  std::vector<var> y0_var1 = stan::math::to_var(y0);
  std::vector<var> theta_var2 = stan::math::to_var(theta);
  std::vector<var> y0_var2 = stan::math::to_var(y0);

  std::vector<std::vector<stan::math::var>> y1 = integrate_ode_rk45(f, y0_var1, t0, ts, theta_var1, x_r, x_i);
  std::vector<std::vector<stan::math::var>> y2 = pmx_integrate_ode_rk45(f, y0_var2, t0, ts, theta_var2, x_r, x_i);
  torsten::test::test_grad(y0_var1, y0_var2, y1, y2, 1.E-10, 1.E-8);
  torsten::test::test_grad(theta_var1, theta_var2, y1, y2, 1.E-10, 1.E-8);
}

TEST_F(TorstenOdeTest_lorenz, fwd_sensitivity_theta_y0) {
  std::vector<var> theta_var1 = stan::math::to_var(theta);
  std::vector<var> y0_var1 = stan::math::to_var(y0);
  std::vector<var> theta_var2 = stan::math::to_var(theta);
  std::vector<var> y0_var2 = stan::math::to_var(y0);

  std::vector<std::vector<stan::math::var>> y1 = integrate_ode_rk45(f, y0_var1, t0, ts, theta_var1, x_r, x_i);
  std::vector<std::vector<stan::math::var>> y2 = pmx_integrate_ode_rk45(f, y0_var2, t0, ts, theta_var2, x_r, x_i);
  torsten::test::test_grad(y0_var1, y0_var2, y1, y2, 1.E-10, 1.E-8);
  torsten::test::test_grad(theta_var1, theta_var2, y1, y2, 1.E-10, 1.E-8);
}

TEST_F(TorstenOdeTest_chem, fwd_sensitivity_theta_y0) {
  std::vector<var> theta_var1 = stan::math::to_var(theta);
  std::vector<var> y0_var1 = stan::math::to_var(y0);
  std::vector<var> theta_var2 = stan::math::to_var(theta);
  std::vector<var> y0_var2 = stan::math::to_var(y0);

  std::vector<std::vector<stan::math::var>> y1 = integrate_ode_rk45(f, y0_var1, t0, ts, theta_var1, x_r, x_i);
  std::vector<std::vector<stan::math::var>> y2 = pmx_integrate_ode_rk45(f, y0_var2, t0, ts, theta_var2, x_r, x_i);
  torsten::test::test_grad(y0_var1, y0_var2, y1, y2, 1.E-10, 1.E-8);
  torsten::test::test_grad(theta_var1, theta_var2, y1, y2, 1.E-10, 1.E-8);
}

// test sequential run of group solver
TEST_F(TorstenOdeTest_neutropenia, group_rk45_fwd_sensitivity_theta) {
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

  vector<vector<var> > y = pmx_integrate_ode_rk45(f, y0, t0, ts, theta_var, x_r, x_i);
  Eigen::Matrix<var, -1, -1> y_m1 = pmx_integrate_ode_group_rk45(f, y0_m, t0, len, ts_m, theta_var_m , x_r_m, x_i_m);
  Eigen::Matrix<var, -1, -1> y_m2 = pmx_integrate_ode_group_rk45(f, y0_m, t0, len, ts_m, theta_var_m , x_r_m, x_i_m);

  EXPECT_EQ(y_m1.cols(), ts_m.size());
  EXPECT_EQ(y_m2.cols(), ts_m.size());
  int icol = 0;
  for (int i = 0; i < np; ++i) {
    stan::math::matrix_v y_i = y_m1.block(0, icol, y0.size(), len[i]);
    torsten::test::test_grad(theta_var, theta_var_m[i], y, y_i, 1e-16, 1e-16);
    icol += len[i];
  }

  icol = 0;
  for (int i = 0; i < np; ++i) {
    stan::math::matrix_v y_i = y_m2.block(0, icol, y0.size(), len[i]);
    torsten::test::test_grad(theta_var, theta_var_m[i], y, y_i, 1e-16, 1e-16);
    icol += len[i];
  }
}
