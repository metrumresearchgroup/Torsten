#ifndef TEST_UNIT_TORSTEN_PK_TWOCPT_MODEL_TEST_FIXTURE
#define TEST_UNIT_TORSTEN_PK_TWOCPT_MODEL_TEST_FIXTURE

#include <stan/math/rev/mat.hpp>
#include <gtest/gtest.h>
#include <boost/numeric/odeint.hpp>
#include <stan/math/torsten/pk_twocpt_model.hpp>
#include <stan/math/torsten/pk_twocpt_solver.hpp>
#include <stan/math/torsten/pk_ode_model.hpp>
#include <stan/math/torsten/pk_rate_adaptor.hpp>
#include <iostream>
#include <sstream>
#include <vector>
#include <string>

struct TorstenTwoCptModelTest : public testing::Test {
  double t0;
  std::vector<double> ts;
  Eigen::Matrix<double, 1, Eigen::Dynamic> y0;
  std::vector<double> rate;
  double CL;
  double Q;
  double V2;
  double V3;
  double ka;
  const std::vector<double> x_r;
  const std::vector<int> x_i;
  std::ostream* msgs;

  void SetUp() {
    // make sure memory's clean before starting each test
    stan::math::recover_memory();
  }
  TorstenTwoCptModelTest() :
    t0(0.0),
    ts{0.1, 0.5, 1.0},
    y0(3),
    rate(3, 0.0),
    CL(10.0),
    Q(28.0),
    V2(80.0),
    V3(70.0),
    ka(1.2),
    msgs{nullptr} {
      y0 << 0.0, 0.0, 0.0;
      SetUp();
  }
};

#endif
