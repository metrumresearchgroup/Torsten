#ifndef TEST_UNIT_TORSTEN_EFFCPT_MODEL_TEST_FIXTURE
#define TEST_UNIT_TORSTEN_EFFCPT_MODEL_TEST_FIXTURE

#include <gtest/gtest.h>
#include <boost/numeric/odeint.hpp>
#include <stan/math/torsten/pmx_onecpt_model.hpp>
#include <stan/math/torsten/pmx_twocpt_model.hpp>
#include <stan/math/torsten/pmx_linode_model.hpp>
#include <stan/math/torsten/pmx_ode_model.hpp>
#include <stan/math/torsten/test/unit/pmx_ode_test_fixture.hpp>
#include <iostream>
#include <sstream>
#include <vector>
#include <string>

struct onecpt_effcpt_model_test : public TorstenOdeTest {
  std::vector<double> ts;
  torsten::PKRec<double> y0;
  std::vector<double> rate;
  double CL;
  double V;
  double ka;
  double ke;
  std::vector<double> par;

  onecpt_effcpt_model_test() :
    ts{0.1, 0.5, 1.0},
    y0(3),
    rate(3, 0.0),
    CL(10.0),
    V(80.0),
    ka(1.2),
    ke(0.5),
    par{CL, V, ka, ke}
  {}
};

#endif
