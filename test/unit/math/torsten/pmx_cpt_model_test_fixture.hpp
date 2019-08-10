#ifndef TEST_UNIT_TORSTEN_PK_CPT_MODEL_TEST_FIXTURE
#define TEST_UNIT_TORSTEN_PK_CPT_MODEL_TEST_FIXTURE

#include <stan/math/rev/mat.hpp>
#include <gtest/gtest.h>
#include <boost/numeric/odeint.hpp>
#include <stan/math/torsten/pmx_onecpt_model.hpp>
#include <stan/math/torsten/pmx_twocpt_model.hpp>
#include <stan/math/torsten/pmx_linode_model.hpp>
#include <stan/math/torsten/pmx_ode_model.hpp>
#include <test/unit/math/torsten/pmx_ode_test_fixture.hpp>
#include <iostream>
#include <sstream>
#include <vector>
#include <string>

struct TorstenTwoCptModelTest : public TorstenOdeTest {

  std::vector<double> ts;
  Eigen::Matrix<double, 1, Eigen::Dynamic> y0;
  std::vector<double> rate;
  double CL;
  double Q;
  double V2;
  double V3;
  double ka;
  std::vector<double> par;
  Eigen::MatrixXd linode_par;      // LinOdeModel parameters

  TorstenTwoCptModelTest() :
    ts{0.1, 0.5, 1.0},
    y0(3),
    rate(3, 0.0),
    CL(10.0),
    Q(28.0),
    V2(80.0),
    V3(70.0),
    ka(1.2),
    par{CL, Q, V2, V3, ka},
    linode_par(3, 3) {
      y0 << 0.0, 0.0, 0.0;

      // to test LinOdeModel, we generate 2-cpt model's
      // linear system.
      double k10 = CL / V2;
      double k12 = Q / V2;
      double k21 = Q / V3;

      /* two-cpt linear system
       * | -ka            0,         0 |
       * |  ka    -(k10 + k12)     k21 |
       * |   0            k12     -k21 |
       */
      // std::vector<double> v_par
      // {-ka, ka, 0.0, 0.0, -(k10 + k12), k12, 0.0, k21, -k21};
      linode_par << -ka, ka, 0.0, 0.0, -(k10 + k12), k12, 0.0, k21, -k21;
  }
};

#endif
