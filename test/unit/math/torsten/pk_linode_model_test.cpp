#include <stan/math.hpp>
#include <stan/math/rev/core.hpp>
#include <test/unit/math/rev/mat/fun/util.hpp>
#include <test/unit/math/torsten/pk_cpt_model_test_fixture.hpp>
#include <test/unit/util.hpp>
#include <gtest/gtest.h>
#include <stan/math/torsten/pk_ode_model.hpp>

TEST_F(TorstenCptOdeModelTest, linode_rate_dbl) {
  using stan::math::var;
  using stan::math::to_var;
  using refactor::PKLinODEModel;
  using refactor::PKLinODE;
  using refactor::PKOdeFunctorRateAdaptor;
  using stan::math::integrate_ode_bdf;
  using torsten::dsolve::pk_integrate_ode_bdf;

  rate[0] = 1200;
  rate[1] = 200;
  rate[2] = 400;
  using model_t = PKLinODEModel<double, double, double, double>;
  model_t model(t0, y0, rate, linode_par);
  std::vector<double> yvec(y0.data(), y0.data() + y0.size());
  PKOdeFunctorRateAdaptor<PKLinODE, double> f1(model.f());

  std::vector<double> y = f1(t0, yvec, model.par(), rate, x_i, msgs);
  EXPECT_FLOAT_EQ(y[0], rate[0]);
  EXPECT_FLOAT_EQ(y[1], rate[1]);
  EXPECT_FLOAT_EQ(y[2], rate[2]);
}

TEST_F(TorstenCptOdeModelTest, linode_rate_var) {
  using stan::math::var;
  using stan::math::to_var;
  using refactor::PKLinODEModel;
  using refactor::PKLinODE;
  using refactor::PKOdeFunctorRateAdaptor;
  using stan::math::integrate_ode_bdf;
  using torsten::dsolve::pk_integrate_ode_bdf;

  rate[0] = 1200;
  rate[1] = 200;
  rate[2] = 300;
  std::vector<stan::math::var> rate_var{to_var(rate)};
  std::vector<stan::math::var> par_var(to_var(linode_par));
  using model_t = PKLinODEModel<double, double, var, var>;
  model_t model(t0, y0, rate_var, par_var);
  std::vector<double> yvec(y0.data(), y0.data() + y0.size());
  PKOdeFunctorRateAdaptor<PKLinODE, var> f1(model.f(), par_var.size());

  par_var.insert(par_var.end(), rate_var.begin(), rate_var.end());
  std::vector<var> y = f1(t0, yvec, par_var, x_r, x_i, msgs);
  EXPECT_FLOAT_EQ(y[0].val(), rate[0]);
  EXPECT_FLOAT_EQ(y[1].val(), rate[1]);
  EXPECT_FLOAT_EQ(y[2].val(), rate[2]);
}

TEST_F(TorstenCptOdeModelTest, linode_solver) {
  using stan::math::var;
  using stan::math::to_var;
  using refactor::PKLinODEModel;
  using refactor::PKLinODE;
  using refactor::PKOdeFunctorRateAdaptor;
  using stan::math::integrate_ode_bdf;
  using torsten::dsolve::pk_integrate_ode_bdf;
  using Eigen::Matrix;
  using Eigen::Dynamic;
  using stan::math::matrix_exp;

  rate[0] = 1200;
  rate[1] = 200;
  rate[2] = 400;
  y0[0] = 150;
  y0[1] = 500;
  y0[2] = 800;
  ts[0] = 10.0;
  ts.resize(1);
  std::vector<stan::math::var> theta{to_var(linode_par)};  
  std::vector<stan::math::var> rate_var{to_var(rate)};
  using model_t = PKLinODEModel<double, double, var, var>;
  model_t model(t0, y0, rate_var, theta);
  std::vector<double> yvec(y0.data(), y0.data() + y0.size());
  PKOdeFunctorRateAdaptor<PKLinODE, var> f1(model.f(), theta.size());

  theta.insert(theta.end(), rate_var.begin(), rate_var.end());
  auto y1 = pk_integrate_ode_bdf(f1, yvec, t0, ts, theta, x_r, x_i, msgs);
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
      EXPECT_NEAR(g1[j], g2[j], 1.E-5);
    }
  }

  for (int i = 0; i < y0.size(); ++i) {
    stan::math::set_zero_all_adjoints();    
    y1[0][i].grad(rate_var, g1);
    stan::math::set_zero_all_adjoints();    
    y2(i).grad(rate_var, g2);
    for (size_t j = 0; j < rate.size(); ++j) {
      EXPECT_NEAR(g1[j], g2[j], 1.E-5);
    }
  }
}

TEST_F(TorstenCptOdeModelTest, linode_solver_zero_rate) {
  using stan::math::var;
  using stan::math::to_var;
  using refactor::PKLinODEModel;
  using refactor::PKLinODE;
  using refactor::PKOdeFunctorRateAdaptor;
  using stan::math::integrate_ode_bdf;
  using torsten::dsolve::pk_integrate_ode_adams;
  using torsten::dsolve::pk_integrate_ode_bdf;
  using Eigen::Matrix;
  using Eigen::Dynamic;
  using stan::math::matrix_exp;

  y0[0] = 150;
  y0[1] = 500;
  y0[2] = 800;
  ts[0] = 20.0;
  ts.resize(1);
  std::vector<stan::math::var> theta{to_var(linode_par)};  
  using model_t = PKLinODEModel<double, double, double, var>;
  model_t model(t0, y0, rate, theta);
  std::vector<double> yvec(y0.data(), y0.data() + y0.size());
  PKOdeFunctorRateAdaptor<PKLinODE, double> f1(model.f());

  auto y1 = pk_integrate_ode_bdf(f1, yvec, t0, ts, theta, rate, x_i, msgs);
  auto y2 = model.solve(ts[0]);
  EXPECT_NEAR(y1[0][0].val(), y2(0).val(), 1.E-7);
  EXPECT_NEAR(y1[0][1].val(), y2(1).val(), 1.E-7);
  EXPECT_NEAR(y1[0][2].val(), y2(2).val(), 1.E-7);

  std::vector<double> g1, g2;
  for (int i = 0; i < y0.size(); ++i) {
    stan::math::set_zero_all_adjoints();    
    y1[0][i].grad(theta, g1);
    stan::math::set_zero_all_adjoints();    
    y2(i).grad(theta, g2);
    for (size_t j = 0; j < theta.size(); ++j) {
      EXPECT_NEAR(g1[j], g2[j], 1.E-6);
    }
  }
}
