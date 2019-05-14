#include <stan/math.hpp>
#include <stan/math/rev/core.hpp>
#include <test/unit/math/rev/mat/fun/util.hpp>
#include <test/unit/math/torsten/pmx_onecpt_model_test_fixture.hpp>
#include <test/unit/util.hpp>
#include <gtest/gtest.h>

TEST_F(TorstenOneCptModelTest, rate_dbl) {
  using stan::math::var;
  using stan::math::to_var;
  using refactor::PMXOneCptModel;
  using torsten::pmx_integrate_ode_bdf;
  using stan::math::integrate_ode_bdf;
  using refactor::PMXOneCptODE;
  using refactor::PMXOdeFunctorRateAdaptor;

  rate[0] = 1200;
  rate[1] = 200;
  using model_t = PMXOneCptModel<double, double, double, double>;
  model_t model(t0, y0, rate, CL, V2, ka);
  std::vector<double> yvec(y0.data(), y0.data() + y0.size());
  PMXOdeFunctorRateAdaptor<PMXOneCptODE, double> f1(model.f());

  std::vector<double> y = f1(t0, yvec, model.par(), rate, x_i, msgs);
  EXPECT_FLOAT_EQ(y[0], rate[0]);
  EXPECT_FLOAT_EQ(y[1], rate[1]);
  EXPECT_FALSE(torsten::has_var_rate<model_t>::value);
}

TEST_F(TorstenOneCptModelTest, rate_var) {
  using stan::math::var;
  using stan::math::to_var;
  using refactor::PMXOneCptModel;
  using stan::math::integrate_ode_bdf;
  using torsten::pmx_integrate_ode_bdf;
  using refactor::PMXOneCptODE;
  using refactor::PMXOdeFunctorRateAdaptor;

  rate[0] = 1200;
  rate[1] = 200;
  var CLv = to_var(CL);
  var V2v = to_var(V2);
  var kav = to_var(ka);
  std::vector<stan::math::var> rate_var{to_var(rate)};
  using model_t = PMXOneCptModel<double, double, var, var>;
  model_t model(t0, y0, rate_var, CLv, V2v, kav);
  std::vector<double> yvec(y0.data(), y0.data() + y0.size());
  std::vector<stan::math::var> theta(model.par());
  PMXOdeFunctorRateAdaptor<PMXOneCptODE, var> f1(model.f(), theta.size());
  theta.insert(theta.end(), rate_var.begin(), rate_var.end());

  std::vector<var> y = f1(t0, yvec, theta, x_r, x_i, msgs);
  EXPECT_FLOAT_EQ(y[0].val(), rate[0]);
  EXPECT_FLOAT_EQ(y[1].val(), rate[1]);

  EXPECT_TRUE(torsten::has_var_rate<model_t>::value);
}

TEST_F(TorstenOneCptModelTest, rate_var_y0) {
  using stan::math::var;
  using stan::math::to_var;
  using refactor::PMXOneCptModel;
  using stan::math::integrate_ode_bdf;
  using torsten::pmx_integrate_ode_bdf;
  using refactor::PMXOneCptODE;
  using refactor::PMXOdeFunctorRateAdaptor;

  rate[0] = 1200;
  rate[1] = 200;
  y0[0] = 150;
  y0[1] = 50;
  var CLv = to_var(CL);
  var V2v = to_var(V2);
  var kav = to_var(ka);
  std::vector<stan::math::var> rate_var{to_var(rate)};
  using model_t = PMXOneCptModel<double, double, var, var>;
  model_t model(t0, y0, rate_var, CLv, V2v, kav);
  std::vector<stan::math::var> theta(model.par());
  std::vector<double> yvec(y0.data(), y0.data() + y0.size());
  PMXOdeFunctorRateAdaptor<PMXOneCptODE, var> f1(model.f(), theta.size());
  theta.insert(theta.end(), rate_var.begin(), rate_var.end());

  std::vector<var> y = f1(t0, yvec, theta, x_r, x_i, msgs);

  EXPECT_TRUE(torsten::has_var_rate<model_t>::value);
  EXPECT_FALSE(torsten::has_var_init<model_t>::value);
  EXPECT_TRUE(torsten::has_var_par<model_t>::value);
  EXPECT_TRUE((std::is_same<torsten::f_t<model_t>, refactor::PMXOneCptODE>::value));
}

TEST_F(TorstenOneCptModelTest, onecpt_solver) {
  using stan::math::var;
  using stan::math::to_var;
  using refactor::PMXOneCptModel;
  using stan::math::integrate_ode_bdf;
  using torsten::pmx_integrate_ode_bdf;
  using refactor::PMXOneCptODE;
  using refactor::PMXOdeFunctorRateAdaptor;

  rate[0] = 1200;
  rate[1] = 200;
  y0[0] = 150;
  y0[1] = 50;
  ts[0] = 10.0;
  ts.resize(1);
  var CLv = to_var(CL);
  var V2v = to_var(V2);
  var kav = to_var(ka);
  std::vector<stan::math::var> rate_var{to_var(rate)};
  using model_t = PMXOneCptModel<double, double, var, var>;
  model_t model(t0, y0, rate_var, CLv, V2v, kav);
  std::vector<stan::math::var> theta(model.par());
  std::vector<double> yvec(y0.data(), y0.data() + y0.size());
  PMXOdeFunctorRateAdaptor<PMXOneCptODE, var> f1(model.f(), theta.size());
  theta.insert(theta.end(), rate_var.begin(), rate_var.end());

  auto y1 = pmx_integrate_ode_bdf(f1, yvec, t0, ts, theta, x_r, x_i, msgs);
  auto y2 = model.solve(ts[0]);
  EXPECT_FLOAT_EQ(y1[0][0].val(), y2(0).val());
  EXPECT_FLOAT_EQ(y1[0][1].val(), y2(1).val());

  std::vector<double> g1, g2;
  for (int i = 0; i < y0.size(); ++i) {
    stan::math::set_zero_all_adjoints();    
    y1[0][i].grad(theta, g1);
    stan::math::set_zero_all_adjoints();    
    y2(i).grad(theta, g2);
    for (size_t j = 0; j < theta.size(); ++j) {
      EXPECT_FLOAT_EQ(g1[j], g2[j]);
    }
  }

  for (int i = 0; i < y0.size(); ++i) {
    stan::math::set_zero_all_adjoints();    
    y1[0][i].grad(rate_var, g1);
    stan::math::set_zero_all_adjoints();    
    y2(i).grad(rate_var, g2);
    for (size_t j = 0; j < rate.size(); ++j) {
      EXPECT_FLOAT_EQ(g1[j], g2[j]);
    }
  }
}

TEST_F(TorstenOneCptModelTest, ss_solver_bolus) {
  using stan::math::var;
  using stan::math::to_var;
  using refactor::PMXOneCptModel;
  using stan::math::integrate_ode_bdf;
  using torsten::pmx_integrate_ode_bdf;
  using refactor::PMXOneCptODE;
  using refactor::PMXOdeFunctorRateAdaptor;

  rate[0] = 0;
  rate[1] = 0;
  y0[0] = 150;
  y0[1] = 50;
  ts[0] = 10.0;
  ts.resize(1);
  var CLv = to_var(CL);
  var V2v = to_var(V2);
  var kav = to_var(ka);
  std::vector<stan::math::var> rate_var{to_var(rate)};
  using model_t = PMXOneCptModel<double, double, var, var>;
  model_t model(t0, y0, rate_var, CLv, V2v, kav);
  std::vector<stan::math::var> theta(model.par());

  std::cout.precision(12);

  double amt = 1800;
  int cmt = 1;
  double ii = 12.0;
  
  auto y1 = model.solve(amt, rate_var[cmt - 1], ii, cmt);
  EXPECT_FLOAT_EQ(y1(0).val(), 1.00330322392E-3);
  EXPECT_FLOAT_EQ(y1(1).val(), 2.07672937446E+0);

  std::vector<double> g1, g2;
  stan::math::set_zero_all_adjoints();
  y1(0).grad(theta, g1);
  EXPECT_FLOAT_EQ(g1[0], 0.0);
  EXPECT_FLOAT_EQ(g1[1], 0.0);
  EXPECT_FLOAT_EQ(g1[2], -0.0120396453978);
  stan::math::set_zero_all_adjoints();
  y1(1).grad(theta, g1);
  EXPECT_FLOAT_EQ(g1[0], -0.266849753086);
  EXPECT_FLOAT_EQ(g1[1], 0.166781095679);
  EXPECT_FLOAT_EQ(g1[2], -1.8559692314);

  cmt = 2;
  y1 = model.solve(amt, rate_var[cmt - 1], ii, cmt);
  EXPECT_FLOAT_EQ(y1(0).val(), 0);
  EXPECT_FLOAT_EQ(y1(1).val(), 0.996102795153);

  stan::math::set_zero_all_adjoints();
  y1(0).grad(theta, g1);
  EXPECT_FLOAT_EQ(g1[0], 0.0);
  EXPECT_FLOAT_EQ(g1[1], 0.0);
  EXPECT_FLOAT_EQ(g1[2], 0.0);
  stan::math::set_zero_all_adjoints();
  y1(1).grad(theta, g1);
  EXPECT_FLOAT_EQ(g1[0], -0.149498104338);
  EXPECT_FLOAT_EQ(g1[1], 0.0934363152112);
  EXPECT_FLOAT_EQ(g1[2], 0.0);
}

TEST_F(TorstenOneCptModelTest, ss_solver_multi_truncated_infusion) {
  using stan::math::var;
  using stan::math::to_var;
  using refactor::PMXOneCptModel;
  using stan::math::integrate_ode_bdf;
  using torsten::pmx_integrate_ode_bdf;
  using refactor::PMXOneCptODE;
  using refactor::PMXOdeFunctorRateAdaptor;

  rate[0] = 1100;
  rate[1] = 770;
  y0[0] = 150;
  y0[1] = 50;
  ts[0] = 10.0;
  ts.resize(1);
  var CLv = to_var(CL);
  var V2v = to_var(V2);
  var kav = to_var(ka);
  std::vector<stan::math::var> rate_var{to_var(rate)};
  using model_t = PMXOneCptModel<double, double, var, var>;
  model_t model(t0, y0, rate_var, CLv, V2v, kav);
  std::vector<var> theta{model.par()};

  double amt = 1800;
  int cmt = 1;
  double ii = 12.0;
  
  auto y1 = model.solve(amt, rate_var[cmt - 1], ii, cmt);
  EXPECT_FLOAT_EQ(y1(0).val(), 0.00312961339574);
  EXPECT_FLOAT_EQ(y1(1).val(), 3.61310672484);

  std::vector<double> g1, g2;
  stan::math::set_zero_all_adjoints();
  y1(0).grad(theta, g1);
  EXPECT_FLOAT_EQ(g1[0], 0.0);
  EXPECT_FLOAT_EQ(g1[1], 0.0);
  EXPECT_FLOAT_EQ(g1[2], -0.0342061212685);
  stan::math::set_zero_all_adjoints();
  y1(1).grad(theta, g1);
  EXPECT_FLOAT_EQ(g1[0], -0.421478621857);
  EXPECT_FLOAT_EQ(g1[1], 0.26342413866);
  EXPECT_FLOAT_EQ(g1[2], -3.20135491073);

  cmt = 2;
  y1 = model.solve(amt, rate_var[cmt - 1], ii, cmt);
  EXPECT_FLOAT_EQ(y1(0).val(), 0.0);
  EXPECT_FLOAT_EQ(y1(1).val(), 2.25697891686);

  stan::math::set_zero_all_adjoints();
  y1(0).grad(theta, g1);
  EXPECT_FLOAT_EQ(g1[0], 0.0);
  EXPECT_FLOAT_EQ(g1[1], 0.0);
  EXPECT_FLOAT_EQ(g1[2], 0.0);
  stan::math::set_zero_all_adjoints();
  y1(1).grad(theta, g1);
  EXPECT_FLOAT_EQ(g1[0], -0.298001025911);
  EXPECT_FLOAT_EQ(g1[1], 0.186250641194);
  EXPECT_FLOAT_EQ(g1[2], 0.0);
}

TEST_F(TorstenOneCptModelTest, ss_solver_const_infusion) {
  using stan::math::var;
  using stan::math::to_var;
  using refactor::PMXOneCptModel;
  using refactor::PMXOdeFunctorRateAdaptor;
  using stan::math::integrate_ode_bdf;
  using torsten::pmx_integrate_ode_bdf;

  rate[0] = 1100;
  rate[1] = 770;
  y0[0] = 150;
  y0[1] = 50;
  ts[0] = 10.0;
  ts.resize(1);
  var CLv = to_var(CL);
  var V2v = to_var(V2);
  var kav = to_var(ka);
  std::vector<var> theta{CLv, V2v, kav};
  std::vector<stan::math::var> rate_var{to_var(rate)};
  using model_t = PMXOneCptModel<double, double, var, var>;
  model_t model(t0, y0, rate_var, CLv, V2v, kav);

  std::cout.precision(12);

  double amt = 1800;
  int cmt = 1;
  double ii = 0.0;
  
  auto y1 = model.solve(amt, rate_var[cmt - 1], ii, cmt);
  EXPECT_FLOAT_EQ(y1(0).val(), 916.666666667);
  EXPECT_FLOAT_EQ(y1(1).val(), 1760);

  std::vector<double> g1, g2;
  stan::math::set_zero_all_adjoints();
  y1(0).grad(theta, g1);
  EXPECT_FLOAT_EQ(g1[0], 0.0);
  EXPECT_FLOAT_EQ(g1[1], 0.0);
  EXPECT_FLOAT_EQ(g1[2], -763.888888889);
  stan::math::set_zero_all_adjoints();
  y1(1).grad(theta, g1);
  EXPECT_FLOAT_EQ(g1[0], -35.2);
  EXPECT_FLOAT_EQ(g1[1], 22);
  EXPECT_FLOAT_EQ(g1[2], 0.0);

  cmt = 2;
  y1 = model.solve(amt, rate_var[cmt - 1], ii, cmt);
  EXPECT_FLOAT_EQ(y1(0).val(), 0.0);
  EXPECT_FLOAT_EQ(y1(1).val(), 1232);

  stan::math::set_zero_all_adjoints();
  y1(0).grad(theta, g1);
  EXPECT_FLOAT_EQ(g1[0], 0.0);
  EXPECT_FLOAT_EQ(g1[1], 0.0);
  EXPECT_FLOAT_EQ(g1[2], 0.0);
  stan::math::set_zero_all_adjoints();
  y1(1).grad(theta, g1);
  EXPECT_FLOAT_EQ(g1[0], -24.64);
  EXPECT_FLOAT_EQ(g1[1], 15.4);
  EXPECT_FLOAT_EQ(g1[2], 0.0);
}
