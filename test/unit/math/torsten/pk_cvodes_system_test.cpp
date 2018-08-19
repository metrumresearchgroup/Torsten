#include <stan/math.hpp>
#include <stan/math/rev/core.hpp>
#include <test/unit/math/rev/mat/fun/util.hpp>
#include <test/unit/math/torsten/pk_ode_test_fixture.hpp>
#include <test/unit/math/prim/arr/functor/harmonic_oscillator.hpp>
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

TEST_F(TorstenOdeTest_sho, pk_cvodes_system) {
  using torsten::dsolve::pk_cvodes_system;
  using torsten::dsolve::pk_cvodes_fwd_system;
  using torsten::dsolve::cvodes_service;
  using stan::math::var;
  using stan::math::to_var;

  std::vector<stan::math::var> theta_var{to_var(theta)};
  std::vector<stan::math::var> y0_var{to_var(y0)};

  using Ode1 = pk_cvodes_fwd_system<F, double, double, double, CV_ADAMS>;
  
  cvodes_service<typename Ode1::Ode> s1(2, 1);
  Ode1 ode1(s1, f, t0, ts, y0, theta, x_r, x_i, msgs);
  test_cvodes_system(ode1, y0, theta, ts);

  using Ode2 = pk_cvodes_fwd_system<F, double, var, double, CV_BDF>;
  
  cvodes_service<typename Ode2::Ode> s2(2, 1);
  Ode2 ode2(s2, f, t0, ts, y0_var, theta, x_r, x_i, msgs);
  test_cvodes_system(ode2, y0, theta, ts);

  using Ode3 = pk_cvodes_fwd_system<F, double, var, var, CV_BDF>;
  cvodes_service<typename Ode3::Ode> s3(2, 1);
  Ode3 ode3(s3, f, t0, ts, y0_var, theta_var, x_r, x_i, msgs);
  Ode3 ode4(s3, f, t0, ts, y0_var, theta_var, x_r, x_i, msgs);
  test_cvodes_system(ode3, y0, theta, ts);
  test_cvodes_system(ode4, y0, theta, ts);
}

TEST_F(TorstenOdeTest_chem, pk_cvodes_system) {
  using torsten::dsolve::pk_cvodes_system;
  using torsten::dsolve::pk_cvodes_fwd_system;
  using torsten::dsolve::cvodes_service;
  using stan::math::var;
  using stan::math::to_var;

  std::vector<stan::math::var> theta_var{to_var(theta)};
  std::vector<stan::math::var> y0_var{to_var(y0)};

  using Ode1 = pk_cvodes_fwd_system<F, double, double, double, CV_ADAMS>;
  
  cvodes_service<typename Ode1::Ode> s1(3, 3);
  Ode1 ode1(s1, f, t0, ts, y0, theta, x_r, x_i, msgs);
  test_cvodes_system(ode1, y0, theta, ts);

  using Ode2 = pk_cvodes_fwd_system<F, double, var, double, CV_BDF>;
  
  cvodes_service<typename Ode2::Ode> s2(3, 3);
  Ode2 ode2(s2, f, t0, ts, y0_var, theta, x_r, x_i, msgs);
  test_cvodes_system(ode2, y0, theta, ts);
}

TEST_F(TorstenOdeTest_lorenz, pk_cvodes_system) {
  using torsten::dsolve::pk_cvodes_system;
  using torsten::dsolve::pk_cvodes_fwd_system;
  using torsten::dsolve::cvodes_service;
  using stan::math::var;
  using stan::math::to_var;

  std::vector<stan::math::var> theta_var{to_var(theta)};
  std::vector<stan::math::var> y0_var{to_var(y0)};

  using Ode1 = pk_cvodes_fwd_system<F, double, double, double, CV_ADAMS>;
  
  cvodes_service<typename Ode1::Ode> s1(3, 3);
  Ode1 ode1(s1, f, t0, ts, y0, theta, x_r, x_i, msgs);
  test_cvodes_system(ode1, y0, theta, ts);

  using Ode2 = pk_cvodes_fwd_system<F, double, var, double, CV_BDF>;
  
  cvodes_service<typename Ode2::Ode> s2(3, 3);
  Ode2 ode2(s2, f, t0, ts, y0_var, theta, x_r, x_i, msgs);
  test_cvodes_system(ode2, y0, theta, ts);
}

TEST_F(TorstenOdeTest_sho, cvodes_constructor_errors) {
  using torsten::dsolve::pk_cvodes_fwd_system;
  using torsten::dsolve::cvodes_service;
  using stan::math::var;
  using stan::math::to_var;

  using Ode1 = pk_cvodes_fwd_system<F, double, double, double, CV_ADAMS>;
  using Ode2 = pk_cvodes_fwd_system<F, double, var, var, CV_ADAMS>;
  
  cvodes_service<typename Ode1::Ode> s1(2, 1);
  cvodes_service<typename Ode2::Ode> s2(2, 1);

  std::vector<double> bad_dbl{y0};
  bad_dbl[0] = std::numeric_limits<double>::infinity();

  ASSERT_NO_THROW(Ode1(s1, f, t0, ts, y0, theta, x_r, x_i, msgs));
  EXPECT_THROW_MSG(Ode1(s1, f, t0, ts, bad_dbl, theta, x_r, x_i, msgs),
                   std::domain_error, "initial state");
  EXPECT_THROW_MSG(Ode1(s1, f, t0, ts, y0, bad_dbl, x_r, x_i, msgs),
                   std::domain_error, "parameter vector");
  EXPECT_THROW_MSG(Ode1(s1, f, t0, ts, y0, theta, bad_dbl, x_i, msgs),
                   std::domain_error, "continuous data");

  std::vector<var> bad_var{std::numeric_limits<double>::infinity(), 1.0, 0.1};
  std::vector<var> empty_var;

  std::vector<stan::math::var> theta_var(to_var(theta));
  std::vector<stan::math::var> y0_var(to_var(y0));
  ASSERT_NO_THROW(Ode2(s2, f, t0, ts, y0_var, theta_var, x_r, x_i, msgs));
  EXPECT_THROW_MSG(Ode2(s2, f, t0, ts, bad_var, theta_var, x_r, x_i, msgs),
                   std::domain_error, "initial state");
  EXPECT_THROW_MSG(Ode2(s2, f, t0, ts, y0_var, bad_var, x_r, x_i, msgs),
                   std::domain_error, "parameter vector");
  EXPECT_THROW_MSG(Ode2(s2, f, t0, ts, empty_var, theta_var, x_r, x_i, msgs),
                   std::invalid_argument, "initial state");
}
