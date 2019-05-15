#include <stan/math.hpp>
#include <stan/math/rev/core.hpp>
#include <stan/math/torsten/dsolve/ode_check.hpp>
#include <test/unit/math/rev/mat/fun/util.hpp>
#include <test/unit/math/torsten/pmx_ode_test_fixture.hpp>
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

TEST_F(TorstenOdeTest_sho, PMXCvodesSystem) {
  using torsten::dsolve::PMXCvodesSystem;
  using torsten::dsolve::PMXCvodesFwdSystem;
  using torsten::dsolve::PMXOdeService;
  using stan::math::var;
  using stan::math::to_var;

  std::vector<stan::math::var> theta_var{to_var(theta)};
  std::vector<stan::math::var> y0_var{to_var(y0)};

  using Ode1 = PMXCvodesFwdSystem<F, double, double, double, CV_ADAMS, CSDA>;
  
  PMXOdeService<typename Ode1::Ode> s1(2, 1);
  Ode1 ode1(s1, f, t0, ts, y0, theta, x_r, x_i, msgs);
  test_cvodes_system(ode1, y0, theta, ts);

  using Ode2 = PMXCvodesFwdSystem<F, double, var, double, CV_BDF, CSDA>;
  
  PMXOdeService<typename Ode2::Ode> s2(2, 1);
  Ode2 ode2(s2, f, t0, ts, y0_var, theta, x_r, x_i, msgs);
  test_cvodes_system(ode2, y0, theta, ts);

  using Ode3 = PMXCvodesFwdSystem<F, double, var, var, CV_BDF, CSDA>;
  PMXOdeService<typename Ode3::Ode> s3(2, 1);
  Ode3 ode3(s3, f, t0, ts, y0_var, theta_var, x_r, x_i, msgs);
  Ode3 ode4(s3, f, t0, ts, y0_var, theta_var, x_r, x_i, msgs);
  test_cvodes_system(ode3, y0, theta, ts);
  test_cvodes_system(ode4, y0, theta, ts);
}

TEST_F(TorstenOdeTest_chem, PMXCvodesSystem) {
  using torsten::dsolve::PMXCvodesSystem;
  using torsten::dsolve::PMXCvodesFwdSystem;
  using torsten::dsolve::PMXOdeService;
  using stan::math::var;
  using stan::math::to_var;

  std::vector<stan::math::var> theta_var{to_var(theta)};
  std::vector<stan::math::var> y0_var{to_var(y0)};

  using Ode1 = PMXCvodesFwdSystem<F, double, double, double, CV_ADAMS, CSDA>;
  
  PMXOdeService<typename Ode1::Ode> s1(3, 3);
  Ode1 ode1(s1, f, t0, ts, y0, theta, x_r, x_i, msgs);
  test_cvodes_system(ode1, y0, theta, ts);

  using Ode2 = PMXCvodesFwdSystem<F, double, var, double, CV_BDF, CSDA>;
  
  PMXOdeService<typename Ode2::Ode> s2(3, 3);
  Ode2 ode2(s2, f, t0, ts, y0_var, theta, x_r, x_i, msgs);
  test_cvodes_system(ode2, y0, theta, ts);
}

TEST_F(TorstenOdeTest_lorenz, PMXCvodesSystem) {
  using torsten::dsolve::PMXCvodesSystem;
  using torsten::dsolve::PMXCvodesFwdSystem;
  using torsten::dsolve::PMXOdeService;
  using stan::math::var;
  using stan::math::to_var;

  std::vector<stan::math::var> theta_var{to_var(theta)};
  std::vector<stan::math::var> y0_var{to_var(y0)};

  using Ode1 = PMXCvodesFwdSystem<F, double, double, double, CV_ADAMS, CSDA>;
  
  PMXOdeService<typename Ode1::Ode> s1(3, 3);
  Ode1 ode1(s1, f, t0, ts, y0, theta, x_r, x_i, msgs);
  test_cvodes_system(ode1, y0, theta, ts);

  using Ode2 = PMXCvodesFwdSystem<F, double, var, double, CV_BDF, CSDA>;
  
  PMXOdeService<typename Ode2::Ode> s2(3, 3);
  Ode2 ode2(s2, f, t0, ts, y0_var, theta, x_r, x_i, msgs);
  test_cvodes_system(ode2, y0, theta, ts);
}

TEST_F(TorstenOdeTest_sho, cvodes_constructor_errors) {
  using torsten::dsolve::PMXCvodesFwdSystem;
  using torsten::dsolve::PMXOdeService;
  using torsten::dsolve::ode_check;
  using stan::math::var;
  using stan::math::to_var;

  using Ode1 = PMXCvodesFwdSystem<F, double, double, double, CV_ADAMS, CSDA>;
  using Ode2 = PMXCvodesFwdSystem<F, double, var, var, CV_ADAMS, CSDA>;
  
  PMXOdeService<typename Ode1::Ode> s1(2, 1);
  PMXOdeService<typename Ode2::Ode> s2(2, 1);

  std::vector<double> bad_dbl{y0};
  bad_dbl[0] = std::numeric_limits<double>::infinity();

  static const char* caller = "";

  ASSERT_NO_THROW(ode_check(y0, t0, ts, theta, x_r, x_i, caller));
  EXPECT_THROW_MSG(ode_check(bad_dbl, t0, ts, theta, x_r, x_i, caller),
                   std::domain_error, ": initial state");
  EXPECT_THROW_MSG(ode_check(y0, t0, ts, bad_dbl, x_r, x_i, caller),
                   std::domain_error, ": parameters");
  EXPECT_THROW_MSG(ode_check(y0, t0, ts, theta, bad_dbl, x_i, caller),
                   std::domain_error, ": continuous data");

  std::vector<var> bad_var{std::numeric_limits<double>::infinity(), 1.0, 0.1};
  std::vector<var> empty_var;

  std::vector<stan::math::var> theta_var(to_var(theta));
  std::vector<stan::math::var> y0_var(to_var(y0));
  ASSERT_NO_THROW(ode_check(y0_var, t0, ts, theta_var, x_r, x_i, caller));
  EXPECT_THROW_MSG(ode_check(bad_var, t0, ts, theta_var, x_r, x_i, caller),
                   std::domain_error, ": initial state");
  EXPECT_THROW_MSG(ode_check(y0_var, t0, ts, bad_var, x_r, x_i, caller),
                   std::domain_error, ": parameters");
  EXPECT_THROW_MSG(ode_check(empty_var, t0, ts, theta_var, x_r, x_i, caller),
                   std::invalid_argument, ": initial state");
}
