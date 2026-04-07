#include <stan/math.hpp>
#include <stan/math/rev/core.hpp>
#include <stan/math/torsten/pmx_onecpt_effcpt_model.hpp>
#include <test/unit/math/rev/fun/util.hpp>
#include <stan/math/torsten/test/unit/pmx_effcpt_model_test_fixture.hpp>
#include <test/unit/util.hpp>
#include <gtest/gtest.h>

using stan::math::var;
using stan::math::to_var;
using stan::math::vector_v;
using stan::math::matrix_v;
using torsten::PMXTwoCptODE;
using torsten::PMXOneCptEffCptModel;
using torsten::dsolve::PMXOdeIntegrator;

using torsten::dsolve::PMXCvodesIntegrator;

TEST_F(onecpt_effcpt_model_test, bolus) {
  y0(0) = 745;
  y0(1) = 100;
  y0(2) = 130;  
  PMXOneCptEffCptModel<double> model_1(CL, V, ka, ke);
  auto par = model_1.to_linode_par();
  PMXLinODEModel<double> model_2(par);
  torsten::PKRec<double> y1(y0), y2(y0), y3(y0);
  model_1.solve(y1, t0, ts[0], rate);
  model_2.solve(y2, t0, ts[0], rate);
  EXPECT_FLOAT_EQ(y1(0), y2(0));
  EXPECT_FLOAT_EQ(y1(1), y2(1));
  EXPECT_FLOAT_EQ(y1(2), y2(2));
}

TEST_F(onecpt_effcpt_model_test, infusion) {
  y0(0) = 745;
  y0(1) = 100;
  y0(2) = 130;  
  rate[0] = 1200;
  rate[1] = 200;
  rate[2] = 0;
  PMXOneCptEffCptModel<double> model_1(CL, V, ka, ke);
  auto par = model_1.to_linode_par();
  PMXLinODEModel<double> model_2(par);
  torsten::PKRec<double> y1(y0), y2(y0);
  model_1.solve(y1, t0, ts[0], rate);
  model_2.solve(y2, t0, ts[0], rate);
  EXPECT_FLOAT_EQ(y1(0), y2(0));
  EXPECT_FLOAT_EQ(y1(1), y2(1));
  EXPECT_FLOAT_EQ(y1(2), y2(2));
}

TEST_F(onecpt_effcpt_model_test, bolus_var) {
  y0(0) = 745;
  y0(1) = 100;
  y0(2) = 130;  
  stan::math::var CL_v(CL);
  stan::math::var V_v(V);
  stan::math::var ka_v(ka);
  stan::math::var ke_v(ke);
  PMXOneCptEffCptModel<stan::math::var> model_1(CL_v, V_v, ka_v, ke_v);
  auto par_mat = model_1.to_linode_par();
  PMXLinODEModel<stan::math::var> model_2(par_mat);
  torsten::PKRec<stan::math::var> y1(y0), y2(y0);
  model_1.solve(y1, t0, ts[0], rate);
  model_2.solve(y2, t0, ts[0], rate);

  std::vector<stan::math::var> theta(model_1.par());
  torsten::test::test_grad(theta, y1, y2, 1.e-8, 1.e-8);
}

TEST_F(onecpt_effcpt_model_test, rate_var) {
  rate[0] = 1200;
  rate[1] = 200;
  rate[2] = 300;
  stan::math::var CL_v(CL);
  stan::math::var V_v(V);
  stan::math::var ka_v(ka);
  stan::math::var ke_v(ke);
  PMXOneCptEffCptModel<stan::math::var> model_1(CL_v, V_v, ka_v, ke_v);
  auto par_mat = model_1.to_linode_par();
  PMXLinODEModel<stan::math::var> model_2(par_mat);
  torsten::PKRec<stan::math::var> y1(y0), y2(y0);
  model_1.solve(y1, t0, ts[0], rate);
  model_2.solve(y2, t0, ts[0], rate);

  std::vector<stan::math::var> theta(model_1.par());
  torsten::test::test_grad(theta, y1, y2, 1.e-8, 1.e-8);
}
