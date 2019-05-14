#include <stan/math.hpp>
#include <stan/math/rev/core.hpp>
#include <test/unit/math/rev/mat/fun/util.hpp>
#include <stan/math/torsten/dsolve/pmx_integrate_ode_adams.hpp>
#include <stan/math/torsten/dsolve/pmx_integrate_ode_bdf.hpp>
#include <test/unit/math/torsten/pmx_ode_test_fixture.hpp>
#include <stan/math/rev/mat/functor/integrate_ode_bdf.hpp>
#include <nvector/nvector_serial.h>
#include <test/unit/util.hpp>
#include <gtest/gtest.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <iostream>
#include <fstream>
#include <vector>
#include <limits>
#include <chrono>
#include <ctime>

using torsten::dsolve::PMXCvodesFwdSystem;
using torsten::dsolve::PMXCvodesIntegrator;
using torsten::dsolve::PMXOdeService;
using torsten::PMXCvodesSensMethod;
using stan::math::var;

TEST_F(TorstenOdeTest_sho, t0_var) {
  PMXCvodesIntegrator solver(rtol, atol, 1000);
  
  ts.resize(1); ts[0] = 1.0;
  std::vector<double> t0_vec{t0};
  
  {
    auto f1 = [&] (std::vector<double>& x) {
      auto y = torsten::pmx_integrate_ode_bdf(f, y0, x[0], ts, theta , x_r, x_i);
      Eigen::MatrixXd y1(1, 2);
      y1(0) = y[0][0];
      y1(1) = y[0][1];
      return y1;
    };
    auto f2 = [&] (std::vector<var>& x) {
      double t0 = stan::math::value_of(x[0]);
      std::vector<var> ts_v{t0 + ts[0] - x[0]};
      auto y = torsten::pmx_integrate_ode_bdf(f, y0, t0, ts_v, theta , x_r, x_i);
      Eigen::Matrix<var, Eigen::Dynamic, Eigen::Dynamic> y1(1, 2);
      y1(0) = y[0][0];
      y1(1) = y[0][1];
      return y1;
    };
    torsten::test::test_grad(f1, f2, t0_vec, 2e-5, 1e-6, 1e-3, 1e-5);
  }
}

TEST_F(TorstenOdeTest_chem, t0_var) {
  PMXCvodesIntegrator solver(rtol, atol, 1000);
  
  ts.resize(1); ts[0] = 1.0;
  std::vector<double> t0_vec{t0};
  
  {
    auto f1 = [&] (std::vector<double>& x) {
      auto y = torsten::pmx_integrate_ode_bdf(f, y0, x[0], ts, theta , x_r, x_i);
      Eigen::MatrixXd y1(1, 3);
      y1(0) = y[0][0];
      y1(1) = y[0][1];
      y1(2) = y[0][2];
      return y1;
    };
    auto f2 = [&] (std::vector<var>& x) {
      double t0 = stan::math::value_of(x[0]);
      std::vector<var> ts_v{t0 + ts[0] - x[0]};
      auto y = torsten::pmx_integrate_ode_bdf(f, y0, t0, ts_v, theta , x_r, x_i);
      Eigen::Matrix<var, Eigen::Dynamic, Eigen::Dynamic> y1(1, 3);
      y1(0) = y[0][0];
      y1(1) = y[0][1];
      y1(2) = y[0][2];
      return y1;
    };
    torsten::test::test_grad(f1, f2, t0_vec, 2e-5, 1e-6, 1e-3, 1e-5);
  }
}

TEST_F(TorstenOdeTest_lorenz, t0_var) {
  PMXCvodesIntegrator solver(rtol, atol, 1000);
  
  ts.resize(1); ts[0] = 1.0;
  std::vector<double> t0_vec{t0};
  
  {
    auto f1 = [&] (std::vector<double>& x) {
      auto y = torsten::pmx_integrate_ode_bdf(f, y0, x[0], ts, theta , x_r, x_i);
      Eigen::MatrixXd y1(1, 3);
      y1(0) = y[0][0];
      y1(1) = y[0][1];
      y1(2) = y[0][2];
      return y1;
    };
    auto f2 = [&] (std::vector<var>& x) {
      double t0 = stan::math::value_of(x[0]);
      std::vector<var> ts_v{t0 + ts[0] - x[0]};
      auto y = torsten::pmx_integrate_ode_bdf(f, y0, t0, ts_v, theta , x_r, x_i);
      Eigen::Matrix<var, Eigen::Dynamic, Eigen::Dynamic> y1(1, 3);
      y1(0) = y[0][0];
      y1(1) = y[0][1];
      y1(2) = y[0][2];
      return y1;
    };
    torsten::test::test_grad(f1, f2, t0_vec, 2e-5, 1e-6, 1e-3, 1e-5);
  }
}

TEST_F(TorstenOdeTest_sho, cvodes_ivp_system_csda) {
  PMXCvodesIntegrator solver(rtol, atol, 1000);

  using Ode1 = PMXCvodesFwdSystem<F, double, double, double, CV_BDF, CSDA>;
  
  PMXOdeService<Ode1> s1(2, 1);
  Ode1 ode{s1, f, t0, ts, y0, theta, x_r, x_i, msgs};
  std::vector<std::vector<double> > y = solver.integrate(ode);
  Eigen::MatrixXd y_mat = solver.integrate<Ode1, false>(ode);
  std::vector<std::vector<double> > y1 =
    stan::math::integrate_ode_bdf(f, y0, t0, ts, theta , x_r, x_i);
  torsten::test::test_val(y, y1);
  torsten::test::test_val(y_mat, y1);

  y = torsten::pmx_integrate_ode_bdf(f, y0, t0, ts, theta , x_r, x_i);
  torsten::test::test_val(y, y1);

  using Ode2 = PMXCvodesFwdSystem<F, double, double, double, CV_ADAMS, CSDA>;
  PMXOdeService<Ode2> s2(2, 1);
  Ode2 ode2{s2, f, t0, ts, y0, theta, x_r, x_i, msgs};
  y = solver.integrate(ode2);
  y_mat = solver.integrate<Ode2, false>(ode2);
  y1 = stan::math::integrate_ode_adams(f, y0, t0, ts, theta , x_r, x_i);
  torsten::test::test_val(y, y1);
  torsten::test::test_val(y_mat, y1);

  y = torsten::pmx_integrate_ode_adams(f, y0, t0, ts, theta , x_r, x_i);
  torsten::test::test_val(y, y1);
}

TEST_F(TorstenOdeTest_lorenz, cvodes_ivp_system) {
  PMXCvodesIntegrator solver(rtol, atol, max_num_steps);

  using Ode1 = PMXCvodesFwdSystem<F, double, double, double, CV_BDF, AD>;
  
  PMXOdeService<typename Ode1::Ode> s1(3, 3);
  Ode1 ode{s1, f, t0, ts, y0, theta, x_r, x_i, msgs};
  std::vector<std::vector<double> > y = solver.integrate(ode);
  Eigen::MatrixXd y_mat = solver.integrate<Ode1, false>(ode);
  std::vector<std::vector<double> > y1 =
    stan::math::integrate_ode_bdf(f, y0, t0, ts, theta , x_r, x_i);
  torsten::test::test_val(y, y1);
  torsten::test::test_val(y_mat, y1);

  y = torsten::pmx_integrate_ode_bdf(f, y0, t0, ts, theta , x_r, x_i);
  torsten::test::test_val(y, y1);

  using Ode2 = PMXCvodesFwdSystem<F, double, double, double, CV_ADAMS, AD>;
  
  PMXOdeService<typename Ode2::Ode> s2(3, 3);
  Ode2 ode2{s2, f, t0, ts, y0, theta, x_r, x_i, msgs};
  y = solver.integrate(ode2);
  y_mat = solver.integrate<Ode2, false>(ode2);
  y1 = stan::math::integrate_ode_adams(f, y0, t0, ts, theta , x_r, x_i);
  torsten::test::test_val(y, y1);
  torsten::test::test_val(y_mat, y1);

  y = torsten::pmx_integrate_ode_adams(f, y0, t0, ts, theta , x_r, x_i);
  torsten::test::test_val(y, y1);
}

TEST_F(TorstenOdeTest_chem, fwd_sensitivity_theta_CSDA) {
  PMXCvodesIntegrator solver(rtol, atol, 1e8);
  std::vector<var> theta_var = stan::math::to_var(theta);

  std::vector<std::vector<var> > y_a, y_b, y1, y2;
  using Ode1 = PMXCvodesFwdSystem<F, double, double, var, CV_ADAMS, CSDA>;
  using Ode2 = PMXCvodesFwdSystem<F, double, double, var, CV_BDF, CSDA>;
  
  PMXOdeService<typename Ode1::Ode> s1(3, 3);
  PMXOdeService<typename Ode2::Ode> s2(3, 3);
  Ode1 ode1(s1, f, t0, ts, y0, theta_var, x_r, x_i, msgs);
  Ode2 ode2(s2, f, t0, ts, y0, theta_var, x_r, x_i, msgs);

  y1 = stan::math::integrate_ode_adams(f, y0, t0, ts, theta_var, x_r, x_i);
  y2 = stan::math::integrate_ode_bdf(f, y0, t0, ts, theta_var, x_r, x_i);

  y_a = solver.integrate(ode1);
  y_b = solver.integrate(ode2);
  torsten::test::test_grad(theta_var, y1, y_a, 1.E-8, 1.E-5);
  torsten::test::test_grad(theta_var, y2, y_b, 1.E-8, 1.E-5);

  Eigen::MatrixXd y_a_mat = solver.integrate<Ode1, false>(ode1);
  Eigen::MatrixXd y_b_mat = solver.integrate<Ode2, false>(ode2);
  torsten::test::test_grad(theta_var, y1, y_a_mat, 1.E-8, 1.E-5);
  torsten::test::test_grad(theta_var, y2, y_b_mat, 1.E-8, 1.E-5);
}

TEST_F(TorstenOdeTest_chem, fwd_sensitivity_theta_AD) {
  PMXCvodesIntegrator solver(rtol, atol, 1e8);
  std::vector<var> theta_var = stan::math::to_var(theta);

  std::vector<std::vector<var> > y_a, y_b, y1, y2;
  using Ode1 = PMXCvodesFwdSystem<F, double, double, var, CV_ADAMS, AD>;
  using Ode2 = PMXCvodesFwdSystem<F, double, double, var, CV_BDF, AD>;
  
  PMXOdeService<Ode1> s1(3, 3);
  PMXOdeService<Ode2> s2(3, 3);
  Ode1 ode1(s1, f, t0, ts, y0, theta_var, x_r, x_i, msgs);
  Ode2 ode2(s2, f, t0, ts, y0, theta_var, x_r, x_i, msgs);

  y1 = stan::math::integrate_ode_adams(f, y0, t0, ts, theta_var, x_r, x_i);
  y2 = stan::math::integrate_ode_bdf(f, y0, t0, ts, theta_var, x_r, x_i);

  y_a = solver.integrate(ode1);
  y_b = solver.integrate(ode2);
  torsten::test::test_grad(theta_var, y1, y_a, 1.E-8, 1.E-5);
  torsten::test::test_grad(theta_var, y2, y_b, 1.E-8, 1.E-5);

  y_a = torsten::pmx_integrate_ode_adams(f, y0, t0, ts, theta_var, x_r, x_i);
  y_b = torsten::pmx_integrate_ode_bdf(f, y0, t0, ts, theta_var, x_r, x_i);
  torsten::test::test_grad(theta_var, y1, y_a, 1.E-8, 1.E-5);
  torsten::test::test_grad(theta_var, y2, y_b, 1.E-8, 1.E-5);

  Eigen::MatrixXd y_a_mat = solver.integrate<Ode1, false>(ode1);
  Eigen::MatrixXd y_b_mat = solver.integrate<Ode2, false>(ode2);
  torsten::test::test_grad(theta_var, y1, y_a_mat, 1.E-8, 1.E-5);
  torsten::test::test_grad(theta_var, y2, y_b_mat, 1.E-8, 1.E-5);
}

TEST_F(TorstenOdeTest_chem, fwd_sensitivity_y0_CSDA) {
  PMXCvodesIntegrator solver(rtol, atol, 1e8);
  std::vector<var> y0_var = stan::math::to_var(y0);

  std::vector<std::vector<var> > y_a, y_b, y1, y2;
  using Ode1 = PMXCvodesFwdSystem<F, double, var, double, CV_ADAMS, CSDA>;
  using Ode2 = PMXCvodesFwdSystem<F, double, var, double, CV_BDF, CSDA>;
  
  PMXOdeService<typename Ode1::Ode> s1(3, 3);
  PMXOdeService<typename Ode2::Ode> s2(3, 3);
  Ode1 ode1(s1, f, t0, ts, y0_var, theta, x_r, x_i, msgs);
  Ode2 ode2(s2, f, t0, ts, y0_var, theta, x_r, x_i, msgs);

  y1 = stan::math::integrate_ode_adams(f, y0_var, t0, ts, theta, x_r, x_i);
  y2 = stan::math::integrate_ode_bdf(f, y0_var, t0, ts, theta, x_r, x_i);

  y_a = solver.integrate(ode1);
  y_b = solver.integrate(ode2);
  torsten::test::test_grad(y0_var, y1, y_a, 1.E-8, 1.E-5);
  torsten::test::test_grad(y0_var, y2, y_b, 1.E-8, 1.E-5);

  Eigen::MatrixXd y_a_mat = solver.integrate<Ode1, false>(ode1);
  Eigen::MatrixXd y_b_mat = solver.integrate<Ode2, false>(ode2);
  torsten::test::test_grad(y0_var, y1, y_a_mat, 1.E-8, 1.E-5);
  torsten::test::test_grad(y0_var, y2, y_b_mat, 1.E-8, 1.E-5);
}

TEST_F(TorstenOdeTest_chem, fwd_sensitivity_y0_AD) {
  PMXCvodesIntegrator solver(rtol, atol, 1e8);
  std::vector<var> y0_var = stan::math::to_var(y0);

  std::vector<std::vector<var> > y_a, y_b, y1, y2;
  using Ode1 = PMXCvodesFwdSystem<F, double, var, double, CV_ADAMS, AD>;
  using Ode2 = PMXCvodesFwdSystem<F, double, var, double, CV_BDF, AD>;
  
  PMXOdeService<Ode1> s1(3, 3);
  PMXOdeService<Ode2> s2(3, 3);
  Ode1 ode1(s1, f, t0, ts, y0_var, theta, x_r, x_i, msgs);
  Ode2 ode2(s2, f, t0, ts, y0_var, theta, x_r, x_i, msgs);

  y1 = stan::math::integrate_ode_adams(f, y0_var, t0, ts, theta, x_r, x_i);
  y2 = stan::math::integrate_ode_bdf(f, y0_var, t0, ts, theta, x_r, x_i);

  y_a = solver.integrate(ode1);
  y_b = solver.integrate(ode2);
  torsten::test::test_grad(y0_var, y_a, y1, 1.E-8, 1.E-5);
  torsten::test::test_grad(y0_var, y_b, y2, 1.E-8, 1.E-5);

  y_a = torsten::pmx_integrate_ode_adams(f, y0_var, t0, ts, theta, x_r, x_i);
  y_b = torsten::pmx_integrate_ode_bdf(f, y0_var, t0, ts, theta, x_r, x_i);
  torsten::test::test_grad(y0_var, y_a, y1, 1.E-8, 1.E-5);
  torsten::test::test_grad(y0_var, y_b, y2, 1.E-8, 1.E-5);

  Eigen::MatrixXd y_a_mat = solver.integrate<Ode1, false>(ode1);
  Eigen::MatrixXd y_b_mat = solver.integrate<Ode2, false>(ode2);
  torsten::test::test_grad(y0_var, y1, y_a_mat, 1.E-8, 1.E-5);
  torsten::test::test_grad(y0_var, y2, y_b_mat, 1.E-8, 1.E-5);
}

TEST_F(TorstenOdeTest_chem, fwd_sensitivity_theta_y0) {
  PMXCvodesIntegrator solver(rtol, atol, 1e8);
  std::vector<var> theta_var = stan::math::to_var(theta);
  std::vector<var> y0_var = stan::math::to_var(y0);

  std::vector<std::vector<var> > y_a, y_b, y1, y2;
  using Ode1 = PMXCvodesFwdSystem<F, double, var, var, CV_ADAMS, AD>;
  using Ode2 = PMXCvodesFwdSystem<F, double, var, var, CV_BDF, AD>;
  
  PMXOdeService<typename Ode1::Ode> s1(3, 3);
  PMXOdeService<typename Ode2::Ode> s2(3, 3);
  Ode1 ode1(s1, f, t0, ts, y0_var, theta_var, x_r, x_i, msgs);
  Ode2 ode2(s2, f, t0, ts, y0_var, theta_var, x_r, x_i, msgs);

  y1 = stan::math::integrate_ode_adams(f, y0_var, t0, ts, theta_var, x_r, x_i);
  y2 = stan::math::integrate_ode_bdf(f, y0_var, t0, ts, theta_var, x_r, x_i);

  y_a = solver.integrate(ode1);
  y_b = solver.integrate(ode2);
  torsten::test::test_grad(y0_var, y_a, y1, 1.E-8, 1.E-5);
  torsten::test::test_grad(y0_var, y_b, y2, 1.E-8, 1.E-5);
  torsten::test::test_grad(theta_var, y_a, y1, 1.E-8, 1.E-5);
  torsten::test::test_grad(theta_var, y_b, y2, 1.E-8, 1.E-5);

  y_a = torsten::pmx_integrate_ode_adams(f, y0_var, t0, ts, theta_var, x_r, x_i);
  y_b = torsten::pmx_integrate_ode_bdf(f, y0_var, t0, ts, theta_var, x_r, x_i);
  torsten::test::test_grad(y0_var, y_a, y1, 1.E-8, 1.E-5);
  torsten::test::test_grad(y0_var, y_b, y2, 1.E-8, 1.E-5);
  torsten::test::test_grad(theta_var, y_a, y1, 1.E-8, 1.E-5);
  torsten::test::test_grad(theta_var, y_b, y2, 1.E-8, 1.E-5);

  Eigen::MatrixXd y_a_mat = solver.integrate<Ode1, false>(ode1);
  Eigen::MatrixXd y_b_mat = solver.integrate<Ode2, false>(ode2);
  std::vector<stan::math::var> vars(y0_var);
  vars.insert(vars.end(), theta_var.begin(), theta_var.end());
  torsten::test::test_grad(vars, y1, y_a_mat, 1.E-8, 1.E-5);
  torsten::test::test_grad(vars, y2, y_b_mat, 1.E-8, 1.E-5);
}

TEST_F(TorstenOdeTest_sho, fwd_sensitivity_theta) {
  PMXCvodesIntegrator solver(rtol, atol, 1e8);
  std::vector<var> theta_var = stan::math::to_var(theta);

  std::vector<std::vector<var> > y_a, y_b, y1, y2;
  using Ode1 = PMXCvodesFwdSystem<F, double, double, var, CV_ADAMS, AD>;
  using Ode2 = PMXCvodesFwdSystem<F, double, double, var, CV_BDF, AD>;
  
  PMXOdeService<typename Ode1::Ode> s1(2, 1);
  PMXOdeService<typename Ode2::Ode> s2(2, 1);
  Ode1 ode1(s1, f, t0, ts, y0, theta_var, x_r, x_i, msgs);
  Ode2 ode2(s2, f, t0, ts, y0, theta_var, x_r, x_i, msgs);

  y1 = stan::math::integrate_ode_adams(f, y0, t0, ts, theta_var, x_r, x_i);
  y2 = stan::math::integrate_ode_bdf(f, y0, t0, ts, theta_var, x_r, x_i);

  y_a = solver.integrate(ode1);
  y_b = solver.integrate(ode2);
  torsten::test::test_grad(theta_var, y_a, y1, 1.E-8, 1.E-5);
  torsten::test::test_grad(theta_var, y_b, y2, 1.E-8, 1.E-5);

  y_a = torsten::pmx_integrate_ode_adams(f, y0, t0, ts, theta_var, x_r, x_i);
  y_b = torsten::pmx_integrate_ode_bdf(f, y0, t0, ts, theta_var, x_r, x_i);
  torsten::test::test_grad(theta_var, y_a, y1, 1.E-8, 1.E-5);
  torsten::test::test_grad(theta_var, y_b, y2, 1.E-8, 1.E-5);

  Eigen::MatrixXd y_a_mat = solver.integrate<Ode1, false>(ode1);
  Eigen::MatrixXd y_b_mat = solver.integrate<Ode2, false>(ode2);
  torsten::test::test_grad(theta_var, y1, y_a_mat, 1.E-8, 1.E-5);
  torsten::test::test_grad(theta_var, y2, y_b_mat, 1.E-8, 1.E-5);
}

TEST_F(TorstenOdeTest_sho, fwd_sensitivity_y0) {
  PMXCvodesIntegrator solver(rtol, atol, 1e8);
  std::vector<var> y0_var = stan::math::to_var(y0);

  std::vector<std::vector<var> > y_a, y_b, y1, y2;
  using Ode1 = PMXCvodesFwdSystem<F, double, var, double, CV_ADAMS, AD>;
  using Ode2 = PMXCvodesFwdSystem<F, double, var, double, CV_BDF, AD>;
  
  PMXOdeService<typename Ode1::Ode> s1(2, 1);
  PMXOdeService<typename Ode2::Ode> s2(2, 1);
  Ode1 ode1(s1, f, t0, ts, y0_var, theta, x_r, x_i, msgs);
  Ode2 ode2(s2, f, t0, ts, y0_var, theta, x_r, x_i, msgs);

  y1 = stan::math::integrate_ode_adams(f, y0_var, t0, ts, theta, x_r, x_i);
  y2 = stan::math::integrate_ode_bdf(f, y0_var, t0, ts, theta, x_r, x_i);

  y_a = solver.integrate(ode1);
  y_b = solver.integrate(ode2);
  torsten::test::test_grad(y0_var, y_a, y1, 1.E-8, 1.E-5);
  torsten::test::test_grad(y0_var, y_b, y2, 1.E-8, 1.E-5);

  y_a = torsten::pmx_integrate_ode_adams(f, y0_var, t0, ts, theta, x_r, x_i);
  y_b = torsten::pmx_integrate_ode_bdf(f, y0_var, t0, ts, theta, x_r, x_i);
  torsten::test::test_grad(y0_var, y_a, y1, 1.E-8, 1.E-5);
  torsten::test::test_grad(y0_var, y_b, y2, 1.E-8, 1.E-5);

  Eigen::MatrixXd y_a_mat = solver.integrate<Ode1, false>(ode1);
  Eigen::MatrixXd y_b_mat = solver.integrate<Ode2, false>(ode2);
  torsten::test::test_grad(y0_var, y1, y_a_mat, 1.E-8, 1.E-5);
  torsten::test::test_grad(y0_var, y2, y_b_mat, 1.E-8, 1.E-5);
}

TEST_F(TorstenOdeTest_sho, fwd_sensitivity_ts) {
  using stan::math::value_of;
  PMXCvodesIntegrator solver(rtol, atol, 1e8);
  std::vector<var> ts_var = stan::math::to_var(ts);

  std::vector<std::vector<var> > y, y1, y2;
  using Ode = PMXCvodesFwdSystem<F, var, double, double, CV_ADAMS, AD>;
  PMXOdeService<typename Ode::Ode> s(2, 1);
  Ode ode(s, f, t0, ts_var, y0, theta, x_r, x_i, msgs);
  y = solver.integrate(ode);

  std::vector<double> g(ts.size()), fval(y0.size());
  for (size_t i = 0; i < ts.size(); ++i) {
    fval = f(value_of(ts[i]), value_of(y[i]), theta, x_r, x_i, msgs);
    for (size_t j = 0; j < y0.size(); ++j) {
      stan::math::set_zero_all_adjoints();
      y[i][j].grad(ts_var, g);
      for (size_t k = 0; k < ts.size(); ++k) {
        if (k == i) {
          EXPECT_FLOAT_EQ(g[k], fval[j]);
        } else {
          EXPECT_FLOAT_EQ(g[k], 0.0);
        }
      }
    }
  }
}

TEST_F(TorstenOdeTest_lorenz, fwd_sensitivity_ts) {
  using stan::math::value_of;

  PMXCvodesIntegrator solver(rtol, atol, 1e8);
  std::vector<var> ts_var = stan::math::to_var(ts);

  std::vector<std::vector<var> > y, y1, y2;
  using Ode = PMXCvodesFwdSystem<F, var, double, double, CV_BDF, AD>;
  PMXOdeService<typename Ode::Ode> s(3, 3);
  Ode ode(s, f, t0, ts_var, y0, theta, x_r, x_i, msgs);
  y = solver.integrate(ode);

  std::vector<double> g(ts.size()), fval(y0.size());
  for (size_t i = 0; i < ts.size(); ++i) {
    fval = f(value_of(ts[i]), value_of(y[i]), theta, x_r, x_i, msgs);
    for (size_t j = 0; j < y0.size(); ++j) {
      stan::math::set_zero_all_adjoints();
      y[i][j].grad(ts_var, g);
      for (size_t k = 0; k < ts.size(); ++k) {
        if (k == i) {
          EXPECT_FLOAT_EQ(g[k], fval[j]);
        } else {
          EXPECT_FLOAT_EQ(g[k], 0.0);
        }
      }
    }
  }
}

TEST_F(TorstenOdeTest_sho, fwd_sensitivity_theta_y0) {
  PMXCvodesIntegrator solver(rtol, atol, 1e8);
  std::vector<var> theta_var = stan::math::to_var(theta);
  std::vector<var> y0_var = stan::math::to_var(y0);

  std::vector<std::vector<var> > y_a, y_b, y1, y2;
  using Ode1 = PMXCvodesFwdSystem<F, double, var, var, CV_ADAMS, AD>;
  using Ode2 = PMXCvodesFwdSystem<F, double, var, var, CV_BDF, AD>;
  
  PMXOdeService<typename Ode1::Ode> s1(2, 1);
  PMXOdeService<typename Ode2::Ode> s2(2, 1);
  Ode1 ode1(s1, f, t0, ts, y0_var, theta_var, x_r, x_i, msgs);
  Ode2 ode2(s2, f, t0, ts, y0_var, theta_var, x_r, x_i, msgs);

  y1 = stan::math::integrate_ode_adams(f, y0_var, t0, ts, theta_var, x_r, x_i);
  y2 = stan::math::integrate_ode_bdf(f, y0_var, t0, ts, theta_var, x_r, x_i);

  y_a = solver.integrate(ode1);
  y_b = solver.integrate(ode2);
  torsten::test::test_grad(y0_var, y_a, y1, 1.E-8, 1.E-5);
  torsten::test::test_grad(y0_var, y_b, y2, 1.E-8, 1.E-5);
  torsten::test::test_grad(theta_var, y_a, y1, 1.E-8, 1.E-5);
  torsten::test::test_grad(theta_var, y_b, y2, 1.E-8, 1.E-5);

  y_a = torsten::pmx_integrate_ode_adams(f, y0_var, t0, ts, theta_var, x_r, x_i);
  y_b = torsten::pmx_integrate_ode_bdf(f, y0_var, t0, ts, theta_var, x_r, x_i);
  torsten::test::test_grad(y0_var, y_a, y1, 1.E-8, 1.E-5);
  torsten::test::test_grad(y0_var, y_b, y2, 1.E-8, 1.E-5);
  torsten::test::test_grad(theta_var, y_a, y1, 1.E-8, 1.E-5);
  torsten::test::test_grad(theta_var, y_b, y2, 1.E-8, 1.E-5);

  Eigen::MatrixXd y_a_mat = solver.integrate<Ode1, false>(ode1);
  Eigen::MatrixXd y_b_mat = solver.integrate<Ode2, false>(ode2);
  std::vector<stan::math::var> vars(y0_var);
  vars.insert(vars.end(), theta_var.begin(), theta_var.end());
  torsten::test::test_grad(vars, y1, y_a_mat, 1.E-8, 1.E-5);
  torsten::test::test_grad(vars, y2, y_b_mat, 1.E-8, 1.E-5);
}

TEST_F(TorstenOdeTest_lorenz, fwd_sensitivity_theta_y0_ts) {
  using stan::math::value_of;

  PMXCvodesIntegrator solver(rtol, atol, 1e8);
  std::vector<var> theta_var = stan::math::to_var(theta);
  std::vector<var> y0_var = stan::math::to_var(y0);
  std::vector<var> ts_var = stan::math::to_var(ts);

  std::vector<std::vector<var> > y, y1, y2;
  using Ode = PMXCvodesFwdSystem<F, var, var, var, CV_ADAMS, AD>;
  PMXOdeService<typename Ode::Ode> s(3, 3);
  Ode ode(s, f, t0, ts_var, y0_var, theta_var, x_r, x_i, msgs);
  y = solver.integrate(ode);

  y1 = stan::math::integrate_ode_adams(f, y0_var, t0, ts, theta_var, x_r, x_i);
  torsten::test::test_grad(y0_var, y, y1, 1.E-6, 1.E-5);
  torsten::test::test_grad(theta_var, y, y1, 1.E-6, 1.E-4);

  std::vector<double> g(ts.size()), fval(y0.size());
  for (size_t i = 0; i < ts.size(); ++i) {
    fval = f(value_of(ts[i]), value_of(y[i]), theta, x_r, x_i, msgs);
    for (size_t j = 0; j < y0.size(); ++j) {
      stan::math::set_zero_all_adjoints();
      y[i][j].grad(ts_var, g);
      for (size_t k = 0; k < ts.size(); ++k) {
        if (k == i) {
          EXPECT_FLOAT_EQ(g[k], fval[j]);
        } else {
          EXPECT_FLOAT_EQ(g[k], 0.0);
        }
      }
    }
  }
}

TEST_F(TorstenOdeTest_chem, fwd_sensitivity_theta_y0_ts) {
  using stan::math::value_of;

  PMXCvodesIntegrator solver(rtol, atol, 1e8);
  std::vector<var> theta_var = stan::math::to_var(theta);
  std::vector<var> y0_var = stan::math::to_var(y0);
  std::vector<var> ts_var = stan::math::to_var(ts);

  std::vector<std::vector<var> > y, y1, y2;
  using Ode = PMXCvodesFwdSystem<F, var, var, var, CV_BDF, AD>;
  PMXOdeService<typename Ode::Ode> s(3, 3);
  Ode ode(s, f, t0, ts_var, y0_var, theta_var, x_r, x_i, msgs);
  y = solver.integrate(ode);

  y1 = stan::math::integrate_ode_bdf(f, y0_var, t0, ts, theta_var, x_r, x_i);

  torsten::test::test_grad(y0_var, y, y1, 1.E-6, 1.E-5);
  torsten::test::test_grad(theta_var, y, y1, 1.E-6, 1.E-4);

  std::vector<double> g(ts.size()), fval(y0.size());
  for (size_t i = 0; i < ts.size(); ++i) {
    fval = f(value_of(ts[i]), value_of(y[i]), theta, x_r, x_i, msgs);
    for (size_t j = 0; j < y0.size(); ++j) {
      stan::math::set_zero_all_adjoints();
      y[i][j].grad(ts_var, g);
      for (size_t k = 0; k < ts.size(); ++k) {
        if (k == i) {
          EXPECT_FLOAT_EQ(g[k], fval[j]);
        } else {
          EXPECT_FLOAT_EQ(g[k], 0.0);
        }
      }
    }
  }
}

TEST_F(TorstenOdeTest_chem, fwd_sensitivity_y0_ts) {
  using stan::math::value_of;

  PMXCvodesIntegrator solver(rtol, atol, 1e8);
  std::vector<var> y0_var = stan::math::to_var(y0);
  std::vector<var> ts_var = stan::math::to_var(ts);

  std::vector<std::vector<var> > y, y1, y2;
  using Ode = PMXCvodesFwdSystem<F, var, var, double, CV_BDF, AD>;
  PMXOdeService<typename Ode::Ode> s(3, 3);
  Ode ode(s, f, t0, ts_var, y0_var, theta, x_r, x_i, msgs);
  y = solver.integrate(ode);

  y1 = stan::math::integrate_ode_bdf(f, y0_var, t0, ts, theta, x_r, x_i);
  torsten::test::test_grad(y0_var, y, y1, 1.E-6, 1.E-5);

  std::vector<double> g(ts.size()), fval(y0.size());
  for (size_t i = 0; i < ts.size(); ++i) {
    fval = f(value_of(ts[i]), value_of(y[i]), theta, x_r, x_i, msgs);
    for (size_t j = 0; j < y0.size(); ++j) {
      stan::math::set_zero_all_adjoints();
      y[i][j].grad(ts_var, g);
      for (size_t k = 0; k < ts.size(); ++k) {
        if (k == i) {
          EXPECT_FLOAT_EQ(g[k], fval[j]);
        } else {
          EXPECT_FLOAT_EQ(g[k], 0.0);
        }
      }
    }
  }
}

TEST_F(TorstenOdeTest_chem, fwd_sensitivity_theta_ts) {
  using stan::math::value_of;

  PMXCvodesIntegrator solver(rtol, atol, 1e8);
  std::vector<var> theta_var = stan::math::to_var(theta);
  std::vector<var> ts_var = stan::math::to_var(ts);

  std::vector<std::vector<var> > y, y1, y2;
  using Ode = PMXCvodesFwdSystem<F, var, double, var, CV_BDF, AD>;
  PMXOdeService<typename Ode::Ode> s(3, 3);
  Ode ode(s, f, t0, ts_var, y0, theta_var, x_r, x_i, msgs);
  y = solver.integrate(ode);

  y1 = stan::math::integrate_ode_bdf(f, y0, t0, ts, theta_var, x_r, x_i);
  torsten::test::test_grad(theta_var, y, y1, 1.E-6, 1.E-5);

  std::vector<double> g(ts.size()), fval(y0.size());
  for (size_t i = 0; i < ts.size(); ++i) {
    fval = f(value_of(ts[i]), value_of(y[i]), theta, x_r, x_i, msgs);
    for (size_t j = 0; j < y0.size(); ++j) {
      stan::math::set_zero_all_adjoints();
      y[i][j].grad(ts_var, g);
      for (size_t k = 0; k < ts.size(); ++k) {
        if (k == i) {
          EXPECT_FLOAT_EQ(g[k], fval[j]);
        } else {
          EXPECT_FLOAT_EQ(g[k], 0.0);
        }
      }
    }
  }
}

TEST_F(TorstenOdeTest_sho, fwd_sensitivity_theta_ts) {
  using stan::math::value_of;

  PMXCvodesIntegrator solver(rtol, atol, 1e8);
  std::vector<var> theta_var = stan::math::to_var(theta);
  std::vector<var> ts_var = stan::math::to_var(ts);

  std::vector<std::vector<var> > y, y1, y2;
  using Ode = PMXCvodesFwdSystem<F, var, double, var, CV_BDF, AD>;
  PMXOdeService<typename Ode::Ode> s(2, 1);
  Ode ode(s, f, t0, ts_var, y0, theta_var, x_r, x_i, msgs);
  y = solver.integrate(ode);

  y1 = stan::math::integrate_ode_bdf(f, y0, t0, ts, theta_var, x_r, x_i);
  torsten::test::test_grad(theta_var, y, y1, 1.E-6, 1.E-5);

  std::vector<double> g(ts.size()), fval(y0.size());
  for (size_t i = 0; i < ts.size(); ++i) {
    fval = f(value_of(ts[i]), value_of(y[i]), theta, x_r, x_i, msgs);
    for (size_t j = 0; j < y0.size(); ++j) {
      stan::math::set_zero_all_adjoints();
      y[i][j].grad(ts_var, g);
      for (size_t k = 0; k < ts.size(); ++k) {
        if (k == i) {
          EXPECT_FLOAT_EQ(g[k], fval[j]);
        } else {
          EXPECT_FLOAT_EQ(g[k], 0.0);
        }
      }
    }
  }
}

TEST_F(TorstenOdeTest_lorenz, fwd_sens_theta_performance_adams) {
  using torsten::pmx_integrate_ode_adams;
  using torsten::pmx_integrate_ode_bdf;

  std::vector<var> theta_var = stan::math::to_var(theta);
  std::vector<std::vector<var> > y1, y2;

  std::chrono::time_point<std::chrono::system_clock> start, end;
  std::chrono::duration<double> y1_elapsed, y2_elapsed;
 
  start = std::chrono::system_clock::now();
  for (int i = 0; i < 100; ++i) {
    y1 = pmx_integrate_ode_adams(f, y0, t0, ts, theta_var, x_r, x_i);
  }
  end = std::chrono::system_clock::now();
  y1_elapsed = end - start;
  std::cout << "torsten solver elapsed time: " << y1_elapsed.count() << "s\n";

  start = std::chrono::system_clock::now();
  for (int i = 0; i < 100; ++i) {
    y2 = stan::math::integrate_ode_adams(f, y0, t0, ts, theta_var, x_r, x_i);
  }
  end = std::chrono::system_clock::now();
  y2_elapsed = end - start;
  std::cout << "stan    solver elapsed time: " << y2_elapsed.count() << "s\n";

  torsten::test::test_grad(theta_var, y1, y2, 1.E-6, 4.E-5);  
}

TEST_F(TorstenOdeTest_chem, fwd_sens_theta_performance_bdf) {
  using torsten::pmx_integrate_ode_adams;
  using torsten::pmx_integrate_ode_bdf;

  std::vector<var> theta_var = stan::math::to_var(theta);
  std::vector<std::vector<var> > y1, y2;

  std::chrono::time_point<std::chrono::system_clock> start, end;
  std::chrono::duration<double> y1_elapsed, y2_elapsed;
 
  start = std::chrono::system_clock::now();
  y1 = pmx_integrate_ode_bdf(f, y0, t0, ts, theta_var, x_r, x_i);
  end = std::chrono::system_clock::now();
  y1_elapsed = end - start;
  std::cout << "torsten solver elapsed time: " << y1_elapsed.count() << "s\n";

  start = std::chrono::system_clock::now();
  y2 = stan::math::integrate_ode_bdf(f, y0, t0, ts, theta_var, x_r, x_i);
  end = std::chrono::system_clock::now();
  y2_elapsed = end - start;
  std::cout << "stan    solver elapsed time: " << y2_elapsed.count() << "s\n";

  torsten::test::test_grad(theta_var, y1, y2, 1.E-6, 4.E-5);  
}

TEST_F(TorstenOdeTest_chem, fwd_sens_theta_performance_repeated) {
  using torsten::pmx_integrate_ode_adams;
  using torsten::pmx_integrate_ode_bdf;
  using stan::math::var;

  std::vector<var> theta_var = stan::math::to_var(theta);
  std::vector<std::vector<var> > y1, y2;

  std::chrono::time_point<std::chrono::system_clock> start, end;
  std::chrono::duration<double> y1_elapsed, y2_elapsed;
 
  const int nloop = 100;

  start = std::chrono::system_clock::now();
  for (int i = 0; i < nloop; ++i) {
    y1 = pmx_integrate_ode_bdf(f, y0, t0, ts, theta_var, x_r, x_i);    
  }
  end = std::chrono::system_clock::now();
  y1_elapsed = end - start;
  std::cout << "torsten solver elapsed time: " << y1_elapsed.count() << "s\n";

  start = std::chrono::system_clock::now();
  for (int i = 0; i < nloop; ++i) {
    y2 = stan::math::integrate_ode_bdf(f, y0, t0, ts, theta_var, x_r, x_i);
  }
  end = std::chrono::system_clock::now();
  y2_elapsed = end - start;
  std::cout << "stan    solver elapsed time: " << y2_elapsed.count() << "s\n";

  torsten::test::test_grad(theta_var, y1, y2, 1.E-6, 4.E-5);  
}

TEST_F(TorstenOdeTest_chem, fwd_sens_y0_theta_performance_repeated) {
  using torsten::pmx_integrate_ode_adams;
  using torsten::pmx_integrate_ode_bdf;
  using stan::math::var;

  std::vector<var> y0_var = stan::math::to_var(y0);
  std::vector<var> theta_var = stan::math::to_var(theta);
  std::vector<std::vector<var> > y1, y2;

  std::chrono::time_point<std::chrono::system_clock> start, end;
  std::chrono::duration<double> y1_elapsed, y2_elapsed;
 
  const int nloop = 100;

  start = std::chrono::system_clock::now();
  for (int i = 0; i < nloop; ++i) {
    y1 = pmx_integrate_ode_bdf(f, y0_var, t0, ts, theta_var, x_r, x_i);
  }
  end = std::chrono::system_clock::now();
  y1_elapsed = end - start;
  std::cout << "torsten solver elapsed time: " << y1_elapsed.count() << "s\n";

  start = std::chrono::system_clock::now();
  for (int i = 0; i < nloop; ++i) {
    y2 = stan::math::integrate_ode_bdf(f, y0_var, t0, ts, theta_var, x_r, x_i);
  }
  end = std::chrono::system_clock::now();
  y2_elapsed = end - start;
  std::cout << "stan    solver elapsed time: " << y2_elapsed.count() << "s\n";

  torsten::test::test_grad(theta_var, y1, y2, 1.E-6, 4.E-5);  
}

TEST_F(TorstenOdeTest_sho, integrate_ode_adams_theta_ts) {
  using torsten::pmx_integrate_ode_adams;
  using torsten::pmx_integrate_ode_bdf;
  using stan::math::value_of;

  std::vector<var> theta_var = stan::math::to_var(theta);
  std::vector<var> ts_var = stan::math::to_var(ts);

  std::vector<std::vector<var> > y, y1, y2;
  for (int i = 0; i < 10; ++i) {
    y = pmx_integrate_ode_adams(f, y0, t0, ts_var, theta_var, x_r, x_i);
  }
  y1 = stan::math::integrate_ode_adams(f, y0, t0, ts, theta_var, x_r, x_i);

  torsten::test::test_grad(theta_var, y, y1, 1.E-6, 1.E-5);

  std::vector<double> g(ts.size()), fval(y0.size());
  for (size_t i = 0; i < ts.size(); ++i) {
    fval = f(value_of(ts[i]), value_of(y[i]), theta, x_r, x_i, msgs);
    for (size_t j = 0; j < y0.size(); ++j) {
      stan::math::set_zero_all_adjoints();
      y[i][j].grad(ts_var, g);
      for (size_t k = 0; k < ts.size(); ++k) {
        if (k == i) {
          EXPECT_FLOAT_EQ(g[k], fval[j]);
        } else {
          EXPECT_FLOAT_EQ(g[k], 0.0);
        }
      }
    }
  }
}

TEST_F(TorstenOdeTest_lorenz, integrate_ode_bdf_y0_ts) {
  using torsten::pmx_integrate_ode_bdf;
  using stan::math::value_of;
  using stan::math::var;

  std::vector<var> y0_var = stan::math::to_var(y0);
  std::vector<var> ts_var = stan::math::to_var(ts);

  std::vector<std::vector<var> > y, y1, y2;
  for (int i = 0; i < 10; ++i) {
    y = pmx_integrate_ode_bdf(f, y0_var, t0, ts_var, theta, x_r, x_i);
  }
  y1 = stan::math::integrate_ode_bdf(f, y0_var, t0, ts, theta, x_r, x_i);

  torsten::test::test_grad(y0_var, y, y1, 1.E-6, 1.E-5);

  std::vector<double> g(ts.size()), fval(y0.size());
  for (size_t i = 0; i < ts.size(); ++i) {
    fval = f(value_of(ts[i]), value_of(y1[i]), theta, x_r, x_i, msgs);
    for (size_t j = 0; j < y0.size(); ++j) {
      stan::math::set_zero_all_adjoints();
      y[i][j].grad(ts_var, g);
      for (size_t k = 0; k < ts.size(); ++k) {
        if (k == i) {
          EXPECT_FLOAT_EQ(g[k], fval[j]);
        } else {
          EXPECT_FLOAT_EQ(g[k], 0.0);
        }
      }
    }
  }
}

TEST_F(TorstenOdeTest_lorenz, integrate_ode_bdf_theta_ts) {
  using torsten::pmx_integrate_ode_bdf;
  using stan::math::value_of;
  using stan::math::var;

  std::vector<var> theta_var = stan::math::to_var(theta);
  std::vector<var> ts_var = stan::math::to_var(ts);

  std::vector<std::vector<var> > y, y1, y2;
  for (int i = 0; i < 10; ++i) {
    y = pmx_integrate_ode_bdf(f, y0, t0, ts_var, theta_var, x_r, x_i);
  }
  y1 = stan::math::integrate_ode_bdf(f, y0, t0, ts, theta_var, x_r, x_i);

  torsten::test::test_grad(theta_var, y, y1, 1.E-6, 1.E-5);

  std::vector<double> g(ts.size()), fval(y0.size());
  for (size_t i = 0; i < ts.size(); ++i) {
    fval = f(value_of(ts[i]), value_of(y1[i]), theta, x_r, x_i, msgs);
    for (size_t j = 0; j < y0.size(); ++j) {
      stan::math::set_zero_all_adjoints();
      y[i][j].grad(ts_var, g);
      for (size_t k = 0; k < ts.size(); ++k) {
        if (k == i) {
          EXPECT_FLOAT_EQ(g[k], fval[j]);
        } else {
          EXPECT_FLOAT_EQ(g[k], 0.0);
        }
      }
    }
  }
}

TEST_F(TorstenOdeTest_lorenz, fwd_sensitivity_theta_AD_stan_bdf) {
  using torsten::dsolve::PMXCvodesFwdSystem;
  using torsten::pmx_integrate_ode_bdf;
  using stan::math::var;
  using std::vector;

  // Lorenz system is chaotic in long term.
  std::vector<double> ts0 {ts};

  vector<var> theta_var1 = stan::math::to_var(theta);
  vector<var> theta_var2 = stan::math::to_var(theta);

  vector<vector<var> > y1 = torsten::pmx_integrate_ode_bdf(f, y0, t0, ts0, theta_var1, x_r, x_i);
  vector<vector<var> > y2 = stan::math::integrate_ode_bdf(f, y0, t0, ts0, theta_var2, x_r, x_i);

  torsten::test::test_grad(theta_var1, theta_var2, y1, y2, 1e-7, 1e-6);
}

// TEST_F(CVODESIntegratorTest, error_handling) {
//   using torsten::dsolve::PMXCvodesFwdSystem;
//   using torsten::dsolve::PMXCvodesIntegrator;
//   const double rtol = 1e-4;
//   const double atol = 1e-8;
//   const int n = 600;

//   ASSERT_NO_THROW(PMXCvodesIntegrator(rtol, atol, n));

//   EXPECT_THROW_MSG(PMXCvodesIntegrator(-1.0E-4, atol, n), std::invalid_argument,
//                    "relative tolerance");

//   EXPECT_THROW_MSG(PMXCvodesIntegrator(2.E-3, atol, n), std::invalid_argument,
//                    "relative tolerance");

//   EXPECT_THROW_MSG(PMXCvodesIntegrator(rtol, -1.E-9, n), std::invalid_argument,
//                    "absolute tolerance");

//   EXPECT_THROW_MSG(PMXCvodesIntegrator(rtol, atol, -100), std::invalid_argument,
//                    "max_num_steps");

//   PMXCvodesFwdSystem<harm_osc_ode_fun, double, double, double> ode{
//       f, y0, theta, x_r, x_i, msgs};
//   PMXCvodesIntegrator solver(rtol, atol, n);
//   double bad_t0 = std::numeric_limits<double>::infinity();
//   std::vector<double> bad_ts{std::numeric_limits<double>::infinity()};
//   EXPECT_THROW_MSG(solver.integrate(ode, bad_t0, ts), std::domain_error,
//                    "initial time");
//   EXPECT_THROW_MSG(solver.integrate(ode, t0, bad_ts), std::domain_error,
//                    "times");
//   bad_t0 = 0;
//   bad_ts[0] = -1;
//   EXPECT_THROW_MSG(solver.integrate(ode, bad_t0, bad_ts), std::domain_error,
//                    "initial time");

//   bad_ts[0] = 0.0;
//   bad_ts.push_back(0.0);
//   EXPECT_THROW_MSG(solver.integrate(ode, t0, bad_ts), std::domain_error,
//                    "times");

//   bad_ts.clear();
//   EXPECT_THROW_MSG(solver.integrate(ode, t0, bad_ts), std::invalid_argument,
//                    "times");
// }
