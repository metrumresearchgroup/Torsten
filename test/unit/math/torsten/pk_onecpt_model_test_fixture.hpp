#ifndef TEST_UNIT_TORSTEN_PK_ONECPT_MODEL_TEST_FIXTURE
#define TEST_UNIT_TORSTEN_PK_ONECPT_MODEL_TEST_FIXTURE

#include <stan/math/rev/mat.hpp>
#include <gtest/gtest.h>
#include <boost/numeric/odeint.hpp>
#include <stan/math/torsten/pk_onecpt_model.hpp>
#include <stan/math/torsten/pk_ode_model.hpp>
#include <iostream>
#include <sstream>
#include <vector>
#include <string>

struct TorstenOneCptModelTest : public testing::Test {
  double t0;
  std::vector<double> ts;
  Eigen::Matrix<double, 1, Eigen::Dynamic> y0;
  std::vector<double> rate;
  double CL;
  double V2;
  double ka;
  const std::vector<double> x_r;
  const std::vector<int> x_i;
  std::ostream* msgs;

  void SetUp() {
    // make sure memory's clean before starting each test
    stan::math::recover_memory();
  }
  TorstenOneCptModelTest() :
    t0(0.0),
    ts{0.1, 0.5, 1.0},
    y0(2),
    rate(2, 0.0),
    CL(50.0),
    V2(80.0),
    ka(1.2),
    msgs{nullptr} {
      y0 << 0.0, 0.0;
      SetUp();
  }
};

#endif
