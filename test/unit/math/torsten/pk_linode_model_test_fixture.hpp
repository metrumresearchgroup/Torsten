#ifndef TEST_UNIT_TORSTEN_PK_LINODE_MODEL_TEST_FIXTURE
#define TEST_UNIT_TORSTEN_PK_LINODE_MODEL_TEST_FIXTURE

#include <stan/math/rev/mat.hpp>
#include <gtest/gtest.h>
#include <boost/numeric/odeint.hpp>
#include <stan/math/torsten/pk_linode_model.hpp>
#include <stan/math/torsten/pk_linode_solver.hpp>
#include <stan/math/torsten/pk_ode_model.hpp>
#include <stan/math/torsten/pk_rate_adaptor.hpp>
#include <iostream>
#include <sstream>
#include <vector>
#include <string>

struct TorstenLinOdeModelTest : public testing::Test {
  double t0;
  std::vector<double> ts;
  Eigen::Matrix<double, 1, Eigen::Dynamic> y0;
  std::vector<double> rate;
  double CL;
  double Q;
  double V2;
  double V3;
  double ka;
  std::vector<double> par;
  const std::vector<double> x_r;
  const std::vector<int> x_i;
  std::ostream* msgs;

  void SetUp() {
    // make sure memory's clean before starting each test
    stan::math::recover_memory();
  }

  /* we map two-cpt model's coefs to linear ODE system.
   */
  TorstenLinOdeModelTest() :
    t0(0.0),
    ts{0.1, 0.5, 1.0},
    y0(3),
    rate(3, 0.0),
    CL(10.0),
    Q(28.0),
    V2(80.0),
    V3(70.0),
    ka(1.2),
    par(9),
    msgs{nullptr} {
      y0 << 0.0, 0.0, 0.0;

      double k10 = CL / V2;
      double k12 = Q / V2;
      double k21 = Q / V3;

      /* two-cpt linear system
       * | -ka            0,         0 |
       * |  ka    -(k10 + k12)     k21 |
       * |   0            k12     -k21 |
       */
      std::vector<double> v_par
      {-ka, ka, 0.0, 0.0, -(k10 + k12), k12, 0.0, k21, -k21};
      par = v_par;

      SetUp();
    }
};

#endif
