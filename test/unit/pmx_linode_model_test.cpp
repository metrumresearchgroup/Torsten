#include <stan/math.hpp>
#include <test/unit/math/rev/fun/util.hpp>
#include <stan/math/torsten/test/unit/pmx_cpt_model_test_fixture.hpp>
#include <test/unit/util.hpp>
#include <gtest/gtest.h>

using stan::math::var;
using stan::math::to_var;
using torsten::PMXLinODEModel;
using torsten::PMXLinODE;
using torsten::PMXOdeFunctorRateAdaptor;
using stan::math::integrate_ode_bdf;
using torsten::pmx_integrate_ode_bdf;
using Eigen::Matrix;
using Eigen::Dynamic;
using stan::math::matrix_exp;
using stan::math::to_vector;
using torsten::PKODEModel;
using torsten::dsolve::PMXVariadicOdeSystem;
using torsten::dsolve::PMXCvodesIntegrator;

PMXLinODE f0;

TEST_F(TorstenTwoCptModelTest, linode_dbl) {
  y0(0) = 745;
  y0(1) = 100;
  y0(2) = 130;
  rate[0] = 1200;
  rate[1] = 200;
  rate[2] = 400;
  torsten::PMXTwoCptModel<double> model1(CL, Q, V2, V3, ka);
  PMXLinODEModel<double> model2(linode_par, 3);

  torsten::PKRec<double> y1(y0), y2(y0);
  model1.solve(y1, t0, ts[0], rate);
  model2.solve(y2, t0, ts[0], rate);  
}

TEST_F(TorstenTwoCptModelTest, linode_rate_var) {
  rate[0] = 1200;
  rate[1] = 200;
  rate[2] = 300;
  std::vector<var> par_var(to_var(par));
  torsten::PMXTwoCptModel<var> model1(par_var);
  Eigen::Matrix<var,-1,-1> linode_par_var(to_var(linode_par));
  PMXLinODEModel<var> model2(linode_par_var, y0.size());

  std::vector<stan::math::var> rate_var{to_var(rate)};
  std::vector<var> par_var_vec(par_var.data(), par_var.data() + par_var.size());
  par_var_vec.insert(par_var_vec.end(), rate_var.begin(), rate_var.end());

  torsten::PKRec<var> y1(to_var(y0)), y2(to_var(y0));
  model1.solve_analytical(y1, t0, ts[0], rate_var);
  model2.solve(y2, t0, ts[0], rate_var);
  torsten::test::test_grad(rate_var, y1, y2, 1.5e-12, 1e-10);
}

TEST_F(TorstenTwoCptModelTest, linode_par_var) {
  rate[0] = 1200;
  rate[1] = 200;
  rate[2] = 300;
  std::vector<var> par_var(to_var(par));
  torsten::PMXTwoCptModel<var> model1(par_var);
  Eigen::Matrix<var, -1, -1> linode_par_var(model1.to_linode_par());
  PMXLinODEModel<var> model2(linode_par_var, 3);

  torsten::PKRec<var> y1(to_var(y0)), y2(to_var(y0));
  model1.solve_analytical(y1, t0, ts[0], rate);
  model2.solve(y2, t0, ts[0], rate);
  torsten::test::test_grad(par_var, y1, y2, 1.5e-12, 1e-10);
}

TEST_F(TorstenTwoCptModelTest, linode_solver) {
  rate[0] = 1200;
  rate[1] = 200;
  rate[2] = 400;
  y0[0] = 150;
  y0[1] = 500;
  y0[2] = 800;
  ts[0] = 10.0;
  ts.resize(1);
  Eigen::Matrix<var,-1,-1> theta{to_var(linode_par)};  
  std::vector<stan::math::var> rate_var{to_var(rate)};
  using model_t = PMXLinODEModel<var>;
  model_t model(theta, y0.size());
  std::vector<var> theta_vec(theta.data(), theta.data() + theta.size());
  theta_vec.insert(theta_vec.end(), rate_var.begin(), rate_var.end());

  PMXOdeFunctorRateAdaptor<PMXLinODE> f1(f0);
  auto y1 = pmx_ode_bdf(f1, y0, t0, ts, msgs, theta_vec, rate_var, x_r, x_i);
  torsten::PKRec<var> y2(to_var(y0));
  model.solve(y2, t0, ts[0], rate_var);
  EXPECT_FLOAT_EQ(y1[0][0].val(), y2(0).val());
  EXPECT_FLOAT_EQ(y1[0][1].val(), y2(1).val());

  std::vector<double> g1, g2;
  for (int i = 0; i < y0.size(); ++i) {
    stan::math::set_zero_all_adjoints();    
    y1[0][i].grad(theta_vec, g1);
    stan::math::set_zero_all_adjoints();    
    y2(i).grad(theta_vec, g2);
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

TEST_F(TorstenTwoCptModelTest, linode_solver_zero_rate) {
  y0[0] = 150;
  y0[1] = 500;
  y0[2] = 800;
  ts[0] = 20.0;
  ts.resize(1);
  Eigen::Matrix<var,-1,-1> theta{to_var(linode_par)};
  using model_t = PMXLinODEModel<var>;
  model_t model(theta, y0.size());
  std::vector<var> theta_vec(theta.data(), theta.data() + theta.size());

  PMXOdeFunctorRateAdaptor<PMXLinODE> f1(f0);
  auto y1 = pmx_ode_bdf(f1, y0, t0, ts, msgs, theta_vec, rate, x_r, x_i);
  torsten::PKRec<var> y2(to_var(y0));
  model.solve(y2, t0, ts[0], rate);
  EXPECT_NEAR(y1[0][0].val(), y2(0).val(), 1.E-7);
  EXPECT_NEAR(y1[0][1].val(), y2(1).val(), 1.E-7);
  EXPECT_NEAR(y1[0][2].val(), y2(2).val(), 1.E-7);

  std::vector<double> g1, g2;
  for (int i = 0; i < y0.size(); ++i) {
    stan::math::set_zero_all_adjoints();    
    y1[0][i].grad(theta_vec, g1);
    stan::math::set_zero_all_adjoints();    
    y2(i).grad(theta_vec, g2);
    for (size_t j = 0; j < theta.size(); ++j) {
      EXPECT_NEAR(g1[j], g2[j], 1.E-6);
    }
  }
}

TEST_F(TorstenTwoCptModelTest, linode_ss_par_var) {
  double amt = 1300;
  double r = 200;
  double ii = 7.0;
  std::vector<var> par_var(to_var(par));
  torsten::PMXTwoCptModel<var> model1(par_var);
  Eigen::Matrix<var, -1, -1> linode_par_var(model1.to_linode_par());
  PMXLinODEModel<var> model2(linode_par_var, model1.ncmt());

  const dsolve::PMXAnalyiticalIntegrator integ;
  torsten::PKRec<var> y1 = model1.solve(ts[0], amt, r, ii, 1, integ);
  torsten::PKRec<var> y2 = model2.solve(ts[0], amt, r, ii, 1, integ);
  torsten::test::test_grad(par_var, y1, y2, 1e-11, 1e-10);
}

TEST_F(TorstenTwoCptModelTest, linode_ss_input_var) {
  var amt = 1300;
  var r = 200;
  var ii = 7.0;
  std::vector<var> par_var(to_var(par));
  torsten::PMXTwoCptModel<var> model1(par_var);
  Eigen::Matrix<var, -1, -1> linode_par_var(model1.to_linode_par());
  PMXLinODEModel<var> model2(linode_par_var, model1.ncmt());

  const dsolve::PMXAnalyiticalIntegrator integ;
  torsten::PKRec<var> y1 = model1.solve(ts[0], amt, r, ii, 1, integ);
  torsten::PKRec<var> y2 = model2.solve(ts[0], amt, r, ii, 1, integ);
  std::vector<var> params{amt, r, ii};
  torsten::test::test_grad(params, y1, y2, 1e-11, 1e-10);
}

TEST_F(TorstenTwoCptModelTest, linode_ss_vs_ode) {
  double amt = 1300;
  double r = 500;
  double ii = 5.0;
  std::vector<var> par_var(to_var(par));
  torsten::PMXTwoCptModel<var> model1(par_var);
  Eigen::Matrix<var, -1, -1> linode_par_var(model1.to_linode_par());
  PMXLinODEModel<var> model2(linode_par_var, model1.ncmt());
  std::vector<var> ode_par_var(to_array_1d(linode_par_var));
  PKODEModel<var, PMXLinODE> model3(ode_par_var, model2.ncmt(), model2.f());

  const dsolve::PMXAnalyiticalIntegrator integ2;
  const PMXOdeIntegrator<PMXVariadicOdeSystem, PMXCvodesIntegrator<CV_BDF, CV_STAGGERED>> integ3;
  torsten::PKRec<var> y2 = model2.solve(ts[0], amt, r, ii, 1, integ2);
  torsten::PKRec<var> y3 = model3.solve(ts[0], amt, r, ii, 1, integ3);
  torsten::test::test_grad(par_var, y2, y3, 1e-7, 2e-7);
}

TEST_F(TorstenTwoCptModelTest, linode_long_ss_vs_ode) {
  double amt = 1300;
  double r = 500;
  double ii = 2.1;
  std::vector<var> par_var(to_var(par));
  torsten::PMXTwoCptModel<var> model1(par_var);
  Eigen::Matrix<var, -1, -1> linode_par_var(model1.to_linode_par());
  PMXLinODEModel<var> model2(linode_par_var, model1.ncmt());
  std::vector<var> ode_par_var(to_array_1d(linode_par_var));
  PKODEModel<var, PMXLinODE> model3(ode_par_var, model2.ncmt(), model2.f());

  const dsolve::PMXAnalyiticalIntegrator integ2;
  const PMXOdeIntegrator<PMXVariadicOdeSystem, PMXCvodesIntegrator<CV_BDF, CV_STAGGERED>> integ3;
  torsten::PKRec<var> y2 = model2.solve(ts[0], amt, r, ii, 1, integ2);
  torsten::PKRec<var> y3 = model3.solve(ts[0], amt, r, ii, 1, integ3);
  torsten::test::test_grad(par_var, y2, y3, 5e-6, 1e-5);
}

TEST_F(TorstenTwoCptModelTest, linode_long_long_ss_infusion_vs_ode) {
  double amt = 1300;
  double r = 500;
  double ii = 1.2;
  std::vector<var> par_var(to_var(par));
  torsten::PMXTwoCptModel<var> model1(par_var);
  Eigen::Matrix<var, -1, -1> linode_par_var(model1.to_linode_par());
  PMXLinODEModel<var> model2(linode_par_var, model1.ncmt());
  std::vector<var> ode_par_var(to_array_1d(linode_par_var));
  PKODEModel<var, PMXLinODE> model3(ode_par_var, model2.ncmt(), model2.f());

  const dsolve::PMXAnalyiticalIntegrator integ2;
  const PMXOdeIntegrator<PMXVariadicOdeSystem, PMXCvodesIntegrator<CV_BDF, CV_STAGGERED>> integ3;
  torsten::PKRec<var> y2 = model2.solve(ts[0], amt, r, ii, 1, integ2);
  torsten::PKRec<var> y3 = model3.solve(ts[0], amt, r, ii, 1, integ3);
  torsten::test::test_grad(par_var, y2, y3, 1e-6, 1e-5);
}
