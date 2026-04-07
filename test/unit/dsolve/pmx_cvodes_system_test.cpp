#include <stan/math.hpp>
#include <stan/math/rev/core.hpp>
#include <stan/math/torsten/dsolve/pmx_ode_system.hpp>
#include <test/unit/math/rev/fun/util.hpp>
#include <stan/math/torsten/test/unit/pmx_ode_test_fixture.hpp>
#include <test/unit/math/prim/functor/harmonic_oscillator.hpp>
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


using stan::math::var;
using stan::math::to_var;
using torsten::dsolve::PMXOdeSystem;

TEST_F(TorstenOdeTest_sho, PMXOdeSystem) {
  std::vector<stan::math::var> theta_var{to_var(theta)};
  std::vector<stan::math::var> y0_var{to_var(y0)};

  using Ode1 = PMXOdeSystem<F, double, double, double>;
  Ode1 ode1(f, t0, ts, y0, theta, x_r, x_i, msgs);
  test_cvodes_system(ode1, y0, theta, ts);

  using Ode2 = PMXOdeSystem<F, double, var, double>;
  Ode2 ode2(f, t0, ts, y0_var, theta, x_r, x_i, msgs);
  test_cvodes_system(ode2, y0, theta, ts);

  using Ode3 = PMXOdeSystem<F, double, var, var>;
  Ode3 ode3(f, t0, ts, y0_var, theta_var, x_r, x_i, msgs);
  Ode3 ode4(f, t0, ts, y0_var, theta_var, x_r, x_i, msgs);
  test_cvodes_system(ode3, y0, theta, ts);
  test_cvodes_system(ode4, y0, theta, ts);
}

TEST_F(TorstenOdeTest_chem, PMXOdeSystem) {
  std::vector<stan::math::var> theta_var{to_var(theta)};
  std::vector<stan::math::var> y0_var{to_var(y0)};

  using Ode1 = PMXOdeSystem<F, double, double, double>;
  Ode1 ode1(f, t0, ts, y0, theta, x_r, x_i, msgs);
  test_cvodes_system(ode1, y0, theta, ts);

  using Ode2 = PMXOdeSystem<F, double, var, double>;
  Ode2 ode2(f, t0, ts, y0_var, theta, x_r, x_i, msgs);
  test_cvodes_system(ode2, y0, theta, ts);

  using Ode3 = PMXOdeSystem<F, double, var, var>;
  Ode3 ode3(f, t0, ts, y0_var, theta_var, x_r, x_i, msgs);
  Ode3 ode4(f, t0, ts, y0_var, theta_var, x_r, x_i, msgs);
  test_cvodes_system(ode3, y0, theta, ts);
  test_cvodes_system(ode4, y0, theta, ts);
}

TEST_F(TorstenOdeTest_lorenz, PMXOdeSystem) {
  std::vector<stan::math::var> theta_var{to_var(theta)};
  std::vector<stan::math::var> y0_var{to_var(y0)};

  using Ode1 = PMXOdeSystem<F, double, double, double>;
  Ode1 ode1(f, t0, ts, y0, theta, x_r, x_i, msgs);
  test_cvodes_system(ode1, y0, theta, ts);

  using Ode2 = PMXOdeSystem<F, double, var, double>;
  Ode2 ode2(f, t0, ts, y0_var, theta, x_r, x_i, msgs);
  test_cvodes_system(ode2, y0, theta, ts);

  using Ode3 = PMXOdeSystem<F, double, var, var>;
  Ode3 ode3(f, t0, ts, y0_var, theta_var, x_r, x_i, msgs);
  Ode3 ode4(f, t0, ts, y0_var, theta_var, x_r, x_i, msgs);
  test_cvodes_system(ode3, y0, theta, ts);
  test_cvodes_system(ode4, y0, theta, ts);
}

TEST_F(TorstenOdeTest_sho, cvodes_constructor_errors) {
  using Ode1 = PMXOdeSystem<F, double, double, double>;
  using Ode2 = PMXOdeSystem<F, double, var, var>;
  
  std::vector<double> bad_dbl{y0};
  bad_dbl[0] = std::numeric_limits<double>::infinity();

  ASSERT_NO_THROW(Ode1(f, t0, ts, y0, theta, x_r, x_i, 0));
  EXPECT_THROW_MSG(Ode1(f, t0, ts, bad_dbl, theta, x_r, x_i, 0),
                   std::domain_error, ": initial state");
  EXPECT_THROW_MSG(Ode1(f, t0, ts, y0, bad_dbl, x_r, x_i, 0),
                   std::domain_error, ": parameters");
  EXPECT_THROW_MSG(Ode1(f, t0, ts, y0, theta, bad_dbl, x_i, 0),
                   std::domain_error, ": continuous data");

  std::vector<var> bad_var{std::numeric_limits<double>::infinity(), 1.0, 0.1};
  std::vector<var> empty_var;

  std::vector<stan::math::var> theta_var(to_var(theta));
  std::vector<stan::math::var> y0_var(to_var(y0));
  ASSERT_NO_THROW(Ode2(f, t0, ts, y0_var, theta_var, x_r, x_i, 0));
  EXPECT_THROW_MSG(Ode2(f, t0, ts, bad_var, theta_var, x_r, x_i, 0),
                   std::domain_error, ": initial state");
  EXPECT_THROW_MSG(Ode2(f, t0, ts, y0_var, bad_var, x_r, x_i, 0),
                   std::domain_error, ": parameters");
  EXPECT_THROW_MSG(Ode2(f, t0, ts, empty_var, theta_var, x_r, x_i, 0),
                   std::invalid_argument, ": initial state");
}
