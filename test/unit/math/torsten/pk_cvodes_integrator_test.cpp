#include <stan/math.hpp>
#include <stan/math/rev/core.hpp>
#include <test/unit/math/rev/mat/fun/util.hpp>
#include <stan/math/torsten/dsolve/pk_cvodes_fwd_system.hpp>
#include <stan/math/torsten/dsolve/pk_cvodes_integrator.hpp>
#include <stan/math/torsten/dsolve/pk_integrate_ode_adams.hpp>
#include <stan/math/torsten/dsolve/pk_integrate_ode_bdf.hpp>
#include <test/unit/math/torsten/pk_ode_test_fixture.hpp>
#include <test/unit/math/prim/arr/functor/harmonic_oscillator.hpp>
#include <stan/math/rev/mat/functor/integrate_ode_bdf.hpp>
#include <nvector/nvector_serial.h>
#include <test/unit/util.hpp>
#include <gtest/gtest.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <iostream>
#include <fstream>
#include <vector>
#include <limits>
#include <chrono>
#include <ctime>

TEST_F(TorstenOdeTest_sho, cvodes_ivp_system) {
  using torsten::dsolve::pk_cvodes_fwd_system;
  using torsten::dsolve::pk_cvodes_integrator;
  using torsten::dsolve::cvodes_service;

  pk_cvodes_integrator solver(rtol, atol, 1000);

  using Ode1 = pk_cvodes_fwd_system<F, double, double, double, CV_BDF>;
  auto f1 = Ode1::rhs();
  cvodes_service<typename Ode1::Ode> s1(2, 1, f1);
  Ode1 ode{s1, f, t0, ts, y0, theta, x_r, x_i, msgs};
  std::vector<std::vector<double> > y = solver.integrate(ode);
  std::vector<std::vector<double> > y1 =
    stan::math::integrate_ode_bdf(f, y0, t0, ts, theta , x_r, x_i);
  for (size_t i = 0; i < y1.size(); ++i) {
    for (size_t j = 0; j < y1[0].size(); ++j) {
      EXPECT_FLOAT_EQ(y[i][j], y1[i][j]);
    }
  }

  using Ode2 = pk_cvodes_fwd_system<F, double, double, double, CV_ADAMS>;
  auto f2 = Ode2::rhs();
  cvodes_service<typename Ode2::Ode> s2(2, 1, f2);
  Ode2 ode2{s2, f, t0, ts, y0, theta, x_r, x_i, msgs};
  y = solver.integrate(ode2);
  y1 = stan::math::integrate_ode_adams(f, y0, t0, ts, theta , x_r, x_i);
  for (size_t i = 0; i < y1.size(); ++i) {
    for (size_t j = 0; j < y1[0].size(); ++j) {
      EXPECT_FLOAT_EQ(y[i][j], y1[i][j]);
    }
  }
}

TEST_F(TorstenOdeTest_lorenz, cvodes_ivp_system) {
  using torsten::dsolve::pk_cvodes_fwd_system;
  using torsten::dsolve::pk_cvodes_integrator;
  using torsten::dsolve::cvodes_service;

  pk_cvodes_integrator solver(rtol, atol, max_num_steps);

  using Ode1 = pk_cvodes_fwd_system<F, double, double, double, CV_BDF>;
  auto f1 = Ode1::rhs();
  cvodes_service<typename Ode1::Ode> s1(3, 3, f1);
  Ode1 ode{s1, f, t0, ts, y0, theta, x_r, x_i, msgs};
  std::vector<std::vector<double> > y = solver.integrate(ode);
  std::vector<std::vector<double> > y1 =
    stan::math::integrate_ode_bdf(f, y0, t0, ts, theta , x_r, x_i);
  for (size_t i = 0; i < y1.size(); ++i) {
    for (size_t j = 0; j < y1[0].size(); ++j) {
      EXPECT_FLOAT_EQ(y[i][j], y1[i][j]);
    }
  }

  using Ode2 = pk_cvodes_fwd_system<F, double, double, double, CV_ADAMS>;
  auto f2 = Ode2::rhs();
  cvodes_service<typename Ode2::Ode> s2(3, 3, f2);
  Ode2 ode2{s2, f, t0, ts, y0, theta, x_r, x_i, msgs};
  y = solver.integrate(ode2);
  y1 = stan::math::integrate_ode_adams(f, y0, t0, ts, theta , x_r, x_i);
  for (size_t i = 0; i < y1.size(); ++i) {
    for (size_t j = 0; j < y1[0].size(); ++j) {
      EXPECT_FLOAT_EQ(y[i][j], y1[i][j]);
    }
  }
}

TEST_F(TorstenOdeTest_chem, fwd_sensitivity_theta) {
  using torsten::dsolve::pk_cvodes_system;
  using torsten::dsolve::pk_cvodes_fwd_system;
  using torsten::dsolve::pk_cvodes_integrator;
  using torsten::dsolve::cvodes_service;

  using stan::math::var;

  pk_cvodes_integrator solver(rtol, atol, 1e8);
  std::vector<var> theta_var = stan::math::to_var(theta);

  std::vector<std::vector<var> > y_a, y_b, y1, y2;
  using Ode1 = pk_cvodes_fwd_system<F, double, double, var, CV_ADAMS>;
  using Ode2 = pk_cvodes_fwd_system<F, double, double, var, CV_BDF>;
  auto f1 = Ode1::rhs();
  auto f2 = Ode2::rhs();
  cvodes_service<typename Ode1::Ode> s1(3, 3, f1);
  cvodes_service<typename Ode2::Ode> s2(3, 3, f2);
  Ode1 ode1(s1, f, t0, ts, y0, theta_var, x_r, x_i, msgs);
  Ode2 ode2(s2, f, t0, ts, y0, theta_var, x_r, x_i, msgs);

  y_a = solver.integrate(ode1);
  y_b = solver.integrate(ode2);

  y1 = stan::math::integrate_ode_adams(f, y0, t0, ts, theta_var, x_r, x_i);
  y2 = stan::math::integrate_ode_bdf(f, y0, t0, ts, theta_var, x_r, x_i);

  test_cvodes_fwd_sens(theta_var, y_a, y1, 1.E-8, 1.E-5);
  test_cvodes_fwd_sens(theta_var, y_b, y2, 1.E-8, 1.E-5);
}

TEST_F(TorstenOdeTest_chem, fwd_sensitivity_y0) {
  using torsten::dsolve::pk_cvodes_system;
  using torsten::dsolve::pk_cvodes_fwd_system;
  using torsten::dsolve::pk_cvodes_integrator;
  using torsten::dsolve::cvodes_service;

  using stan::math::var;

  pk_cvodes_integrator solver(rtol, atol, 1e8);
  std::vector<var> y0_var = stan::math::to_var(y0);

  std::vector<std::vector<var> > y_a, y_b, y1, y2;
  using Ode1 = pk_cvodes_fwd_system<F, double, var, double, CV_ADAMS>;
  using Ode2 = pk_cvodes_fwd_system<F, double, var, double, CV_BDF>;
  auto f1 = Ode1::rhs();
  auto f2 = Ode2::rhs();
  cvodes_service<typename Ode1::Ode> s1(3, 3, f1);
  cvodes_service<typename Ode2::Ode> s2(3, 3, f2);
  Ode1 ode1(s1, f, t0, ts, y0_var, theta, x_r, x_i, msgs);
  Ode2 ode2(s2, f, t0, ts, y0_var, theta, x_r, x_i, msgs);

  y_a = solver.integrate(ode1);
  y_b = solver.integrate(ode2);

  y1 = stan::math::integrate_ode_adams(f, y0_var, t0, ts, theta, x_r, x_i);
  y2 = stan::math::integrate_ode_bdf(f, y0_var, t0, ts, theta, x_r, x_i);

  test_cvodes_fwd_sens(y0_var, y_a, y1, 1.E-8, 1.E-5);
  test_cvodes_fwd_sens(y0_var, y_b, y2, 1.E-8, 1.E-5);
}

TEST_F(TorstenOdeTest_chem, fwd_sensitivity_theta_y0) {
  using torsten::dsolve::pk_cvodes_system;
  using torsten::dsolve::pk_cvodes_fwd_system;
  using torsten::dsolve::pk_cvodes_integrator;
  using torsten::dsolve::cvodes_service;

  using stan::math::var;

  pk_cvodes_integrator solver(rtol, atol, 1e8);
  std::vector<var> theta_var = stan::math::to_var(theta);
  std::vector<var> y0_var = stan::math::to_var(y0);

  std::vector<std::vector<var> > y_a, y_b, y1, y2;
  using Ode1 = pk_cvodes_fwd_system<F, double, var, var, CV_ADAMS>;
  using Ode2 = pk_cvodes_fwd_system<F, double, var, var, CV_BDF>;
  auto f1 = Ode1::rhs();
  auto f2 = Ode2::rhs();
  cvodes_service<typename Ode1::Ode> s1(3, 3, f1);
  cvodes_service<typename Ode2::Ode> s2(3, 3, f2);
  Ode1 ode1(s1, f, t0, ts, y0_var, theta_var, x_r, x_i, msgs);
  Ode2 ode2(s2, f, t0, ts, y0_var, theta_var, x_r, x_i, msgs);

  y_a = solver.integrate(ode1);
  y_b = solver.integrate(ode2);

  y1 = stan::math::integrate_ode_adams(f, y0_var, t0, ts, theta_var, x_r, x_i);
  y2 = stan::math::integrate_ode_bdf(f, y0_var, t0, ts, theta_var, x_r, x_i);

  test_cvodes_fwd_sens(y0_var, y_a, y1, 1.E-8, 1.E-5);
  test_cvodes_fwd_sens(y0_var, y_b, y2, 1.E-8, 1.E-5);
  test_cvodes_fwd_sens(theta_var, y_a, y1, 1.E-8, 1.E-5);
  test_cvodes_fwd_sens(theta_var, y_b, y2, 1.E-8, 1.E-5);
}

TEST_F(TorstenOdeTest_sho, fwd_sensitivity_theta) {
  using torsten::dsolve::pk_cvodes_system;
  using torsten::dsolve::pk_cvodes_fwd_system;
  using torsten::dsolve::pk_cvodes_integrator;
  using torsten::dsolve::cvodes_service;

  using stan::math::var;

  pk_cvodes_integrator solver(rtol, atol, 1e8);
  std::vector<var> theta_var = stan::math::to_var(theta);

  std::vector<std::vector<var> > y_a, y_b, y1, y2;
  using Ode1 = pk_cvodes_fwd_system<F, double, double, var, CV_ADAMS>;
  using Ode2 = pk_cvodes_fwd_system<F, double, double, var, CV_BDF>;
  auto f1 = Ode1::rhs();
  auto f2 = Ode2::rhs();
  cvodes_service<typename Ode1::Ode> s1(2, 1, f1);
  cvodes_service<typename Ode2::Ode> s2(2, 1, f2);
  Ode1 ode1(s1, f, t0, ts, y0, theta_var, x_r, x_i, msgs);
  Ode2 ode2(s2, f, t0, ts, y0, theta_var, x_r, x_i, msgs);

  y_a = solver.integrate(ode1);
  y_b = solver.integrate(ode2);

  y1 = stan::math::integrate_ode_adams(f, y0, t0, ts, theta_var, x_r, x_i);
  y2 = stan::math::integrate_ode_bdf(f, y0, t0, ts, theta_var, x_r, x_i);

  test_cvodes_fwd_sens(theta_var, y_a, y1, 1.E-8, 1.E-5);
  test_cvodes_fwd_sens(theta_var, y_b, y2, 1.E-8, 1.E-5);
}

TEST_F(TorstenOdeTest_sho, fwd_sensitivity_y0) {
  using torsten::dsolve::pk_cvodes_system;
  using torsten::dsolve::pk_cvodes_fwd_system;
  using torsten::dsolve::pk_cvodes_integrator;
  using torsten::dsolve::cvodes_service;

  using stan::math::var;

  pk_cvodes_integrator solver(rtol, atol, 1e8);
  std::vector<var> y0_var = stan::math::to_var(y0);

  std::vector<std::vector<var> > y_a, y_b, y1, y2;
  using Ode1 = pk_cvodes_fwd_system<F, double, var, double, CV_ADAMS>;
  using Ode2 = pk_cvodes_fwd_system<F, double, var, double, CV_BDF>;
  auto f1 = Ode1::rhs();
  auto f2 = Ode2::rhs();
  cvodes_service<typename Ode1::Ode> s1(2, 1, f1);
  cvodes_service<typename Ode2::Ode> s2(2, 1, f2);
  Ode1 ode1(s1, f, t0, ts, y0_var, theta, x_r, x_i, msgs);
  Ode2 ode2(s2, f, t0, ts, y0_var, theta, x_r, x_i, msgs);

  y_a = solver.integrate(ode1);
  y_b = solver.integrate(ode2);

  y1 = stan::math::integrate_ode_adams(f, y0_var, t0, ts, theta, x_r, x_i);
  y2 = stan::math::integrate_ode_bdf(f, y0_var, t0, ts, theta, x_r, x_i);

  test_cvodes_fwd_sens(y0_var, y_a, y1, 1.E-8, 1.E-5);
  test_cvodes_fwd_sens(y0_var, y_b, y2, 1.E-8, 1.E-5);
}

TEST_F(TorstenOdeTest_sho, fwd_sensitivity_ts) {
  using torsten::dsolve::pk_cvodes_system;
  using torsten::dsolve::pk_cvodes_fwd_system;
  using torsten::dsolve::pk_cvodes_integrator;
  using torsten::dsolve::cvodes_service;
  using stan::math::value_of;

  using stan::math::var;

  pk_cvodes_integrator solver(rtol, atol, 1e8);
  std::vector<var> ts_var = stan::math::to_var(ts);

  std::vector<std::vector<var> > y, y1, y2;
  using Ode = pk_cvodes_fwd_system<F, var, double, double, CV_ADAMS>;
  cvodes_service<typename Ode::Ode> s(2, 1, Ode::rhs());
  Ode ode(s, f, t0, ts_var, y0, theta, x_r, x_i, msgs);
  y = solver.integrate(ode);

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

TEST_F(TorstenOdeTest_lorenz, fwd_sensitivity_ts) {
  using torsten::dsolve::pk_cvodes_system;
  using torsten::dsolve::pk_cvodes_fwd_system;
  using torsten::dsolve::pk_cvodes_integrator;
  using torsten::dsolve::cvodes_service;
  using stan::math::value_of;

  using stan::math::var;

  pk_cvodes_integrator solver(rtol, atol, 1e8);
  std::vector<var> ts_var = stan::math::to_var(ts);

  std::vector<std::vector<var> > y, y1, y2;
  using Ode = pk_cvodes_fwd_system<F, var, double, double, CV_BDF>;
  cvodes_service<typename Ode::Ode> s(3, 3, Ode::rhs());
  Ode ode(s, f, t0, ts_var, y0, theta, x_r, x_i, msgs);
  y = solver.integrate(ode);

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

TEST_F(TorstenOdeTest_sho, fwd_sensitivity_theta_y0) {
  using torsten::dsolve::pk_cvodes_system;
  using torsten::dsolve::pk_cvodes_fwd_system;
  using torsten::dsolve::pk_cvodes_integrator;
  using torsten::dsolve::cvodes_service;

  using stan::math::var;

  pk_cvodes_integrator solver(rtol, atol, 1e8);
  std::vector<var> theta_var = stan::math::to_var(theta);
  std::vector<var> y0_var = stan::math::to_var(y0);

  std::vector<std::vector<var> > y_a, y_b, y1, y2;
  using Ode1 = pk_cvodes_fwd_system<F, double, var, var, CV_ADAMS>;
  using Ode2 = pk_cvodes_fwd_system<F, double, var, var, CV_BDF>;
  auto f1 = Ode1::rhs();
  auto f2 = Ode2::rhs();
  cvodes_service<typename Ode1::Ode> s1(2, 1, f1);
  cvodes_service<typename Ode2::Ode> s2(2, 1, f2);
  Ode1 ode1(s1, f, t0, ts, y0_var, theta_var, x_r, x_i, msgs);
  Ode2 ode2(s2, f, t0, ts, y0_var, theta_var, x_r, x_i, msgs);

  y_a = solver.integrate(ode1);
  y_b = solver.integrate(ode2);

  y1 = stan::math::integrate_ode_adams(f, y0_var, t0, ts, theta_var, x_r, x_i);
  y2 = stan::math::integrate_ode_bdf(f, y0_var, t0, ts, theta_var, x_r, x_i);

  test_cvodes_fwd_sens(y0_var, y_a, y1, 1.E-8, 1.E-5);
  test_cvodes_fwd_sens(y0_var, y_b, y2, 1.E-8, 1.E-5);
  test_cvodes_fwd_sens(theta_var, y_a, y1, 1.E-8, 1.E-5);
  test_cvodes_fwd_sens(theta_var, y_b, y2, 1.E-8, 1.E-5);
}

TEST_F(TorstenOdeTest_lorenz, fwd_sensitivity_theta_y0_ts) {
  using torsten::dsolve::pk_cvodes_system;
  using torsten::dsolve::pk_cvodes_fwd_system;
  using torsten::dsolve::pk_cvodes_integrator;
  using torsten::dsolve::cvodes_service;
  using stan::math::value_of;
  using stan::math::var;

  pk_cvodes_integrator solver(rtol, atol, 1e8);
  std::vector<var> theta_var = stan::math::to_var(theta);
  std::vector<var> y0_var = stan::math::to_var(y0);
  std::vector<var> ts_var = stan::math::to_var(ts);

  std::vector<std::vector<var> > y, y1, y2;
  using Ode = pk_cvodes_fwd_system<F, var, var, var, CV_ADAMS>;
  cvodes_service<typename Ode::Ode> s(3, 3, Ode::rhs());
  Ode ode(s, f, t0, ts_var, y0_var, theta_var, x_r, x_i, msgs);
  y = solver.integrate(ode);

  y1 = stan::math::integrate_ode_adams(f, y0_var, t0, ts, theta_var, x_r, x_i);

  test_cvodes_fwd_sens(y0_var, y, y1, 1.E-6, 1.E-5);
  test_cvodes_fwd_sens(theta_var, y, y1, 1.E-6, 1.E-4);

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


TEST_F(TorstenOdeTest_chem, fwd_sensitivity_theta_y0_ts) {
  using torsten::dsolve::pk_cvodes_system;
  using torsten::dsolve::pk_cvodes_fwd_system;
  using torsten::dsolve::pk_cvodes_integrator;
  using torsten::dsolve::cvodes_service;
  using stan::math::value_of;
  using stan::math::var;

  pk_cvodes_integrator solver(rtol, atol, 1e8);
  std::vector<var> theta_var = stan::math::to_var(theta);
  std::vector<var> y0_var = stan::math::to_var(y0);
  std::vector<var> ts_var = stan::math::to_var(ts);

  std::vector<std::vector<var> > y, y1, y2;
  using Ode = pk_cvodes_fwd_system<F, var, var, var, CV_BDF>;
  cvodes_service<typename Ode::Ode> s(3, 3, Ode::rhs());
  Ode ode(s, f, t0, ts_var, y0_var, theta_var, x_r, x_i, msgs);
  y = solver.integrate(ode);

  y1 = stan::math::integrate_ode_bdf(f, y0_var, t0, ts, theta_var, x_r, x_i);

  test_cvodes_fwd_sens(y0_var, y, y1, 1.E-6, 1.E-5);
  test_cvodes_fwd_sens(theta_var, y, y1, 1.E-6, 1.E-4);

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

TEST_F(TorstenOdeTest_chem, fwd_sensitivity_y0_ts) {
  using torsten::dsolve::pk_cvodes_system;
  using torsten::dsolve::pk_cvodes_fwd_system;
  using torsten::dsolve::pk_cvodes_integrator;
  using torsten::dsolve::cvodes_service;
  using stan::math::value_of;
  using stan::math::var;

  pk_cvodes_integrator solver(rtol, atol, 1e8);
  std::vector<var> y0_var = stan::math::to_var(y0);
  std::vector<var> ts_var = stan::math::to_var(ts);

  std::vector<std::vector<var> > y, y1, y2;
  using Ode = pk_cvodes_fwd_system<F, var, var, double, CV_BDF>;
  cvodes_service<typename Ode::Ode> s(3, 3, Ode::rhs());
  Ode ode(s, f, t0, ts_var, y0_var, theta, x_r, x_i, msgs);
  y = solver.integrate(ode);

  y1 = stan::math::integrate_ode_bdf(f, y0_var, t0, ts, theta, x_r, x_i);

  test_cvodes_fwd_sens(y0_var, y, y1, 1.E-6, 1.E-5);

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

TEST_F(TorstenOdeTest_chem, fwd_sensitivity_theta_ts) {
  using torsten::dsolve::pk_cvodes_system;
  using torsten::dsolve::pk_cvodes_fwd_system;
  using torsten::dsolve::pk_cvodes_integrator;
  using torsten::dsolve::cvodes_service;
  using stan::math::value_of;
  using stan::math::var;

  pk_cvodes_integrator solver(rtol, atol, 1e8);
  std::vector<var> theta_var = stan::math::to_var(theta);
  std::vector<var> ts_var = stan::math::to_var(ts);

  std::vector<std::vector<var> > y, y1, y2;
  using Ode = pk_cvodes_fwd_system<F, var, double, var, CV_BDF>;
  cvodes_service<typename Ode::Ode> s(3, 3, Ode::rhs());
  Ode ode(s, f, t0, ts_var, y0, theta_var, x_r, x_i, msgs);
  y = solver.integrate(ode);

  y1 = stan::math::integrate_ode_bdf(f, y0, t0, ts, theta_var, x_r, x_i);

  test_cvodes_fwd_sens(theta_var, y, y1, 1.E-6, 1.E-5);

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

TEST_F(TorstenOdeTest_sho, fwd_sensitivity_theta_ts) {
  using torsten::dsolve::pk_cvodes_system;
  using torsten::dsolve::pk_cvodes_fwd_system;
  using torsten::dsolve::pk_cvodes_integrator;
  using torsten::dsolve::cvodes_service;
  using stan::math::value_of;
  using stan::math::var;

  pk_cvodes_integrator solver(rtol, atol, 1e8);
  std::vector<var> theta_var = stan::math::to_var(theta);
  std::vector<var> ts_var = stan::math::to_var(ts);

  std::vector<std::vector<var> > y, y1, y2;
  using Ode = pk_cvodes_fwd_system<F, var, double, var, CV_BDF>;
  cvodes_service<typename Ode::Ode> s(2, 1, Ode::rhs());
  Ode ode(s, f, t0, ts_var, y0, theta_var, x_r, x_i, msgs);
  y = solver.integrate(ode);

  y1 = stan::math::integrate_ode_bdf(f, y0, t0, ts, theta_var, x_r, x_i);

  test_cvodes_fwd_sens(theta_var, y, y1, 1.E-6, 1.E-5);

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

TEST_F(TorstenOdeTest_lorenz, fwd_sens_theta_performance_adams) {
  using torsten::dsolve::pk_integrate_ode_adams;
  using torsten::dsolve::pk_integrate_ode_bdf;
  using stan::math::var;

  std::vector<var> theta_var = stan::math::to_var(theta);
  std::vector<std::vector<var> > y1, y2;

  std::chrono::time_point<std::chrono::system_clock> start, end;
  std::chrono::duration<double> y1_elapsed, y2_elapsed;
 
  start = std::chrono::system_clock::now();
  for (int i = 0; i < 100; ++i) {
    y1 = pk_integrate_ode_adams(f, y0, t0, ts, theta_var, x_r, x_i);
  }
  end = std::chrono::system_clock::now();
  y1_elapsed = end - start;
  std::cout << "torsten solver elapsed time: " << y1_elapsed.count() << "s\n";

  start = std::chrono::system_clock::now();
  for (int i = 0; i < 100; ++i) {
    y2 = stan::math::integrate_ode_adams(f, y0, t0, ts, theta_var, x_r, x_i);
  }
  end = std::chrono::system_clock::now();
  y2_elapsed = end - start;
  std::cout << "stan    solver elapsed time: " << y2_elapsed.count() << "s\n";

  test_cvodes_fwd_sens(theta_var, y1, y2, 1.E-6, 4.E-5);  
}

TEST_F(TorstenOdeTest_chem, fwd_sens_theta_performance_bdf) {
  using torsten::dsolve::pk_integrate_ode_adams;
  using torsten::dsolve::pk_integrate_ode_bdf;
  using stan::math::var;

  std::vector<var> theta_var = stan::math::to_var(theta);
  std::vector<std::vector<var> > y1, y2;

  std::chrono::time_point<std::chrono::system_clock> start, end;
  std::chrono::duration<double> y1_elapsed, y2_elapsed;
 
  start = std::chrono::system_clock::now();
  y1 = pk_integrate_ode_bdf(f, y0, t0, ts, theta_var, x_r, x_i);
  end = std::chrono::system_clock::now();
  y1_elapsed = end - start;
  std::cout << "torsten solver elapsed time: " << y1_elapsed.count() << "s\n";

  start = std::chrono::system_clock::now();
  y2 = stan::math::integrate_ode_bdf(f, y0, t0, ts, theta_var, x_r, x_i);
  end = std::chrono::system_clock::now();
  y2_elapsed = end - start;
  std::cout << "stan    solver elapsed time: " << y2_elapsed.count() << "s\n";

  test_cvodes_fwd_sens(theta_var, y1, y2, 1.E-6, 4.E-5);  
}

TEST_F(TorstenOdeTest_chem, fwd_sens_theta_performance_repeated) {
  using torsten::dsolve::pk_integrate_ode_adams;
  using torsten::dsolve::pk_integrate_ode_bdf;
  using stan::math::var;

  std::vector<var> theta_var = stan::math::to_var(theta);
  std::vector<std::vector<var> > y1, y2;

  std::chrono::time_point<std::chrono::system_clock> start, end;
  std::chrono::duration<double> y1_elapsed, y2_elapsed;
 
  start = std::chrono::system_clock::now();
  for (int i = 0; i < 1000; ++i) {
    y1 = pk_integrate_ode_bdf(f, y0, t0, ts, theta_var, x_r, x_i);    
  }
  end = std::chrono::system_clock::now();
  y1_elapsed = end - start;
  std::cout << "torsten solver elapsed time: " << y1_elapsed.count() << "s\n";

  start = std::chrono::system_clock::now();
  for (int i = 0; i < 1000; ++i) {
    y2 = stan::math::integrate_ode_bdf(f, y0, t0, ts, theta_var, x_r, x_i);
  }
  end = std::chrono::system_clock::now();
  y2_elapsed = end - start;
  std::cout << "stan    solver elapsed time: " << y2_elapsed.count() << "s\n";

  test_cvodes_fwd_sens(theta_var, y1, y2, 1.E-6, 4.E-5);  
}

TEST_F(TorstenOdeTest_chem, fwd_sens_y0_theta_performance_repeated) {
  using torsten::dsolve::pk_integrate_ode_adams;
  using torsten::dsolve::pk_integrate_ode_bdf;
  using stan::math::var;

  std::vector<var> y0_var = stan::math::to_var(y0);
  std::vector<var> theta_var = stan::math::to_var(theta);
  std::vector<std::vector<var> > y1, y2;

  std::chrono::time_point<std::chrono::system_clock> start, end;
  std::chrono::duration<double> y1_elapsed, y2_elapsed;
 
  start = std::chrono::system_clock::now();
  for (int i = 0; i < 1000; ++i) {
    y1 = pk_integrate_ode_bdf(f, y0_var, t0, ts, theta_var, x_r, x_i);
  }
  end = std::chrono::system_clock::now();
  y1_elapsed = end - start;
  std::cout << "torsten solver elapsed time: " << y1_elapsed.count() << "s\n";

  start = std::chrono::system_clock::now();
  for (int i = 0; i < 1000; ++i) {
    y2 = stan::math::integrate_ode_bdf(f, y0_var, t0, ts, theta_var, x_r, x_i);
  }
  end = std::chrono::system_clock::now();
  y2_elapsed = end - start;
  std::cout << "stan    solver elapsed time: " << y2_elapsed.count() << "s\n";

  test_cvodes_fwd_sens(theta_var, y1, y2, 1.E-6, 4.E-5);  
}

TEST_F(TorstenOdeTest_sho, integrate_ode_adams_theta_ts) {
  using torsten::dsolve::pk_integrate_ode_adams;
  using torsten::dsolve::pk_integrate_ode_bdf;
  using stan::math::value_of;
  using stan::math::var;

  std::vector<var> theta_var = stan::math::to_var(theta);
  std::vector<var> ts_var = stan::math::to_var(ts);

  std::vector<std::vector<var> > y, y1, y2;
  for (int i = 0; i < 10; ++i) {
    y = pk_integrate_ode_adams(f, y0, t0, ts_var, theta_var, x_r, x_i);
  }
  y1 = stan::math::integrate_ode_adams(f, y0, t0, ts, theta_var, x_r, x_i);

  test_cvodes_fwd_sens(theta_var, y, y1, 1.E-6, 1.E-5);

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

TEST_F(TorstenOdeTest_lorenz, integrate_ode_bdf_y0_ts) {
  using torsten::dsolve::pk_integrate_ode_bdf;
  using stan::math::value_of;
  using stan::math::var;

  std::vector<var> y0_var = stan::math::to_var(y0);
  std::vector<var> ts_var = stan::math::to_var(ts);

  std::vector<std::vector<var> > y, y1, y2;
  for (int i = 0; i < 10; ++i) {
    y = pk_integrate_ode_bdf(f, y0_var, t0, ts_var, theta, x_r, x_i);
  }
  y1 = stan::math::integrate_ode_bdf(f, y0_var, t0, ts, theta, x_r, x_i);

  test_cvodes_fwd_sens(y0_var, y, y1, 1.E-6, 1.E-5);

  std::vector<double> g(ts.size()), fval(y0.size());
  for (size_t i = 0; i < ts.size(); ++i) {
    fval = f(value_of(ts[i]), value_of(y1[i]), theta, x_r, x_i, msgs);
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

TEST_F(TorstenOdeTest_lorenz, integrate_ode_bdf_theta_ts) {
  using torsten::dsolve::pk_integrate_ode_bdf;
  using stan::math::value_of;
  using stan::math::var;

  std::vector<var> theta_var = stan::math::to_var(theta);
  std::vector<var> ts_var = stan::math::to_var(ts);

  std::vector<std::vector<var> > y, y1, y2;
  for (int i = 0; i < 10; ++i) {
    y = pk_integrate_ode_bdf(f, y0, t0, ts_var, theta_var, x_r, x_i);
  }
  y1 = stan::math::integrate_ode_bdf(f, y0, t0, ts, theta_var, x_r, x_i);

  test_cvodes_fwd_sens(theta_var, y, y1, 1.E-6, 1.E-5);

  std::vector<double> g(ts.size()), fval(y0.size());
  for (size_t i = 0; i < ts.size(); ++i) {
    fval = f(value_of(ts[i]), value_of(y1[i]), theta, x_r, x_i, msgs);
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

// TEST_F(CVODESIntegratorTest, error_handling) {
//   using torsten::dsolve::pk_cvodes_fwd_system;
//   using torsten::dsolve::pk_cvodes_integrator;
//   const double rtol = 1e-4;
//   const double atol = 1e-8;
//   const int n = 600;

//   ASSERT_NO_THROW(pk_cvodes_integrator(rtol, atol, n));

//   EXPECT_THROW_MSG(pk_cvodes_integrator(-1.0E-4, atol, n), std::invalid_argument,
//                    "relative tolerance");

//   EXPECT_THROW_MSG(pk_cvodes_integrator(2.E-3, atol, n), std::invalid_argument,
//                    "relative tolerance");

//   EXPECT_THROW_MSG(pk_cvodes_integrator(rtol, -1.E-9, n), std::invalid_argument,
//                    "absolute tolerance");

//   EXPECT_THROW_MSG(pk_cvodes_integrator(rtol, atol, -100), std::invalid_argument,
//                    "max_num_steps");

//   pk_cvodes_fwd_system<harm_osc_ode_fun, double, double, double> ode{
//       f, y0, theta, x_r, x_i, msgs};
//   pk_cvodes_integrator solver(rtol, atol, n);
//   double bad_t0 = std::numeric_limits<double>::infinity();
//   std::vector<double> bad_ts{std::numeric_limits<double>::infinity()};
//   EXPECT_THROW_MSG(solver.integrate(ode, bad_t0, ts), std::domain_error,
//                    "initial time");
//   EXPECT_THROW_MSG(solver.integrate(ode, t0, bad_ts), std::domain_error,
//                    "times");
//   bad_t0 = 0;
//   bad_ts[0] = -1;
//   EXPECT_THROW_MSG(solver.integrate(ode, bad_t0, bad_ts), std::domain_error,
//                    "initial time");

//   bad_ts[0] = 0.0;
//   bad_ts.push_back(0.0);
//   EXPECT_THROW_MSG(solver.integrate(ode, t0, bad_ts), std::domain_error,
//                    "times");

//   bad_ts.clear();
//   EXPECT_THROW_MSG(solver.integrate(ode, t0, bad_ts), std::invalid_argument,
//                    "times");
// }
