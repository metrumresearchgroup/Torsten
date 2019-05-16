#include <stan/math.hpp>
#include <stan/math/rev/core.hpp>
#include <test/unit/math/rev/mat/fun/util.hpp>
#include <stan/math/torsten/dsolve/pmx_integrate_dae.hpp>
#include <test/unit/math/torsten/dae_systems.hpp>

#include <nvector/nvector_serial.h>

#include <test/unit/util.hpp>
#include <gtest/gtest.h>

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <iostream>
#include <fstream>
#include <vector>

struct StanIntegrateDAETest : public ::testing::Test {
  chemical_kinetics f;
  std::vector<double> yy0;
  std::vector<double> yp0;
  std::vector<double> theta;
  std::vector<double> x_r;
  std::vector<int> x_i;
  std::ostream* msgs;
  const std::vector<int> eq_id;
  const double t0;
  std::vector<double> ts;

  void SetUp() { stan::math::recover_memory(); }

  StanIntegrateDAETest()
      : yy0{1.0, 0.0, 0.0},
        yp0{-0.04, 0.04, 0.0},
        theta{0.040, 1.0e4, 3.0e7},
        msgs{0},
        eq_id{1, 1, 0},
        t0(0.0) {
    const size_t nout{4};
    const double h{0.4};
    for (size_t i = 0; i < nout; ++i)
      ts.push_back(h * std::pow(10, i));
  }
};

TEST_F(StanIntegrateDAETest, idas_ivp_system_yy0) {
  using torsten::dsolve::pk_idas_fwd_system;
  using torsten::dsolve::pmx_integrate_dae;
  std::vector<std::vector<double> > yy
      = pmx_integrate_dae(f, yy0, yp0, t0, ts, theta, x_r, x_i, 1e-4, 1e-8);
  EXPECT_NEAR(0.985172, yy[0][0], 1e-6);
  EXPECT_NEAR(0.0147939, yy[0][2], 1e-6);
  EXPECT_NEAR(0.905521, yy[1][0], 1e-6);
  EXPECT_NEAR(0.0944571, yy[1][2], 1e-6);
}

TEST_F(StanIntegrateDAETest, forward_sensitivity_theta) {
  using torsten::dsolve::pk_idas_fwd_system;
  using torsten::dsolve::idas_integrator;
  using torsten::dsolve::pmx_integrate_dae;
  using stan::math::to_var;
  using stan::math::value_of;
  using stan::math::var;

  std::vector<var> theta_var = to_var(theta);

  std::vector<std::vector<var> > yy
      = pmx_integrate_dae(f, yy0, yp0, t0, ts, theta_var, x_r, x_i, 1e-5, 1e-12);
  EXPECT_NEAR(0.985172, value_of(yy[0][0]), 1e-6);
  EXPECT_NEAR(0.0147939, value_of(yy[0][2]), 1e-6);
  EXPECT_NEAR(0.905519, value_of(yy[1][0]), 1e-6);
  EXPECT_NEAR(0.0944588, value_of(yy[1][2]), 1e-6);

  // test derivatives against central difference results
  std::vector<double> g;
  const double h = 1.e-2;
  const std::vector<double> theta1{theta[0] - theta[0] * h, theta[1], theta[2]};
  const std::vector<double> theta2{theta[0] + theta[0] * h, theta[1], theta[2]};
  torsten::dsolve::idas_integrator solver(1e-5, 1e-12, 1000);
  pk_idas_fwd_system<chemical_kinetics, double, double, double> dae1(
      f, eq_id, yy0, yp0, theta1, x_r, x_i, msgs),
      dae2(f, eq_id, yy0, yp0, theta2, x_r, x_i, msgs);
  std::vector<std::vector<double> > yy1 = solver.integrate(dae1, t0, ts);
  std::vector<std::vector<double> > yy2 = solver.integrate(dae2, t0, ts);

  double yys_finite_diff = (yy2[3][1] - yy1[3][1]) / (2.0 * theta[0] * h);
  stan::math::set_zero_all_adjoints();
  yy[3][1].grad(theta_var, g);
  EXPECT_NEAR(yys_finite_diff, g[0], 1e-6);
}

TEST_F(StanIntegrateDAETest, inconsistent_ic_error) {
  using torsten::dsolve::pmx_integrate_dae;
  using stan::math::to_var;
  using stan::math::value_of;
  using stan::math::var;

  std::vector<var> theta_var = to_var(theta);

  yy0.back() = -0.1;
  EXPECT_THROW_MSG(
      pmx_integrate_dae(f, yy0, yp0, t0, ts, theta_var, x_r, x_i, 1e-5, 1e-12),
      std::domain_error, "DAE residual at t0");
}
