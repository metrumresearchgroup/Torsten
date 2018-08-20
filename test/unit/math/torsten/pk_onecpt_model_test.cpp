#include <stan/math.hpp>
#include <stan/math/rev/core.hpp>
#include <test/unit/math/rev/mat/fun/util.hpp>
#include <test/unit/math/torsten/pk_onecpt_model_test_fixture.hpp>
#include <test/unit/util.hpp>
#include <gtest/gtest.h>

TEST_F(TorstenOneCptModelTest, rate_dbl) {
  using stan::math::var;
  using stan::math::to_var;
  using refactor::PKOneCptModel;
  using refactor::PKODERateAdaptor;
  using stan::math::integrate_ode_bdf;
  using torsten::dsolve::pk_integrate_ode_bdf;

  rate[0] = 1200;
  rate[1] = 200;
  using model_t = PKOneCptModel<double, double, double, double>;
  model_t model(t0, y0, rate, CL, V2, ka);
  std::vector<double> yvec(y0.data(), y0.data() + y0.size());
  PKODERateAdaptor<model_t> rate_adaptor(model);

  std::vector<double> y =
    rate_adaptor.model().f()(t0, yvec, model.par(), rate, x_i, msgs);
  EXPECT_FLOAT_EQ(y[0], rate[0]);
  EXPECT_FLOAT_EQ(y[1], rate[1]);
}

TEST_F(TorstenOneCptModelTest, rate_dbl_y0) {
  using stan::math::var;
  using stan::math::to_var;
  using refactor::PKOneCptModel;
  using refactor::PKODERateAdaptor;
  using stan::math::integrate_ode_bdf;
  using torsten::dsolve::pk_integrate_ode_bdf;

  rate[0] = 1200;
  rate[1] = 200;
  y0[0] = 150;
  y0[1] = 50;
  using model_t = PKOneCptModel<double, double, double, double>;
  model_t model(t0, y0, rate, CL, V2, ka);
  std::vector<double> yvec(y0.data(), y0.data() + y0.size());
  PKODERateAdaptor<model_t> rate_adaptor(model);

  std::vector<double> y =
    rate_adaptor.model().f()(t0, yvec, model.par(),
                             rate_adaptor.model().rate(), x_i, msgs);
  EXPECT_FLOAT_EQ(y[0], -ka * y0[0] + rate[0]);
  EXPECT_FLOAT_EQ(y[1], ka * y0[0] - model.k10() * y0[1] + rate[1]);
}

TEST_F(TorstenOneCptModelTest, rate_var) {
  using stan::math::var;
  using stan::math::to_var;
  using refactor::PKOneCptModel;
  using refactor::PKODERateAdaptor;
  using stan::math::integrate_ode_bdf;
  using torsten::dsolve::pk_integrate_ode_bdf;

  rate[0] = 1200;
  rate[1] = 200;
  var CLv = to_var(CL);
  var V2v = to_var(V2);
  var kav = to_var(ka);
  std::vector<stan::math::var> rate_var{to_var(rate)};
  using model_t = PKOneCptModel<double, double, var, var>;
  model_t model(t0, y0, rate_var, CLv, V2v, kav);
  std::vector<double> yvec(y0.data(), y0.data() + y0.size());
  PKODERateAdaptor<model_t> rate_adaptor(model);

  EXPECT_EQ(rate_adaptor.model().par().size(),
            model.par().size() + model.rate().size());
  std::vector<var> y =
    rate_adaptor.model().f()(t0, yvec, rate_adaptor.model().par(),
                             x_r, x_i, msgs);
  EXPECT_FLOAT_EQ(y[0].val(), rate[0]);
  EXPECT_FLOAT_EQ(y[1].val(), rate[1]);
}

TEST_F(TorstenOneCptModelTest, rate_var_y0) {
  using stan::math::var;
  using stan::math::to_var;
  using refactor::PKOneCptModel;
  using refactor::PKODERateAdaptor;
  using stan::math::integrate_ode_bdf;
  using torsten::dsolve::pk_integrate_ode_bdf;

  rate[0] = 1200;
  rate[1] = 200;
  y0[0] = 150;
  y0[1] = 50;
  var CLv = to_var(CL);
  var V2v = to_var(V2);
  var kav = to_var(ka);
  std::vector<stan::math::var> rate_var{to_var(rate)};
  using model_t = PKOneCptModel<double, double, var, var>;
  model_t model(t0, y0, rate_var, CLv, V2v, kav);
  std::vector<double> yvec(y0.data(), y0.data() + y0.size());
  PKODERateAdaptor<model_t> rate_adaptor(model);

  std::vector<var> y =
    rate_adaptor.model().f()(t0, yvec, rate_adaptor.model().par(),
                             x_r, x_i, msgs);
  EXPECT_FLOAT_EQ(y[0].val(), -ka * y0[0] + rate[0]);
  EXPECT_FLOAT_EQ(y[1].val(), ka * y0[0] - model.k10().val() * y0[1] + rate[1]);
}

