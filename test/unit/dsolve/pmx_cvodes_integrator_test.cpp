#include <stan/math.hpp>
#include <stan/math/rev/core.hpp>
#include <test/unit/math/rev/fun/util.hpp>
#include <stan/math/torsten/dsolve/pmx_integrate_ode_adams.hpp>
#include <stan/math/torsten/dsolve/pmx_integrate_ode_bdf.hpp>
#include <stan/math/torsten/test/unit/pmx_ode_test_fixture.hpp>
#include <stan/math/rev/functor/integrate_ode_bdf.hpp>
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

using torsten::dsolve::PMXOdeSystem;
using torsten::dsolve::PMXCvodesIntegrator;
using torsten::dsolve::PMXOdeService;
using torsten::dsolve::OdeObserver;
using torsten::dsolve::OdeDataObserver;
using torsten::PMXCvodesSensMethod;
using torsten::AD;
using torsten::CSDA;
using torsten::DQ;
using stan::math::var;

TEST_F(TorstenOdeTest_sho, t0_var) {
  ts.resize(1); ts[0] = 1.0;
  std::vector<double> t0_vec{t0};
  
  {
    auto f1 = [&] (const std::vector<double>& x) {
      auto y = torsten::pmx_integrate_ode_bdf(f, y0, x[0], ts, theta , x_r, x_i);
      Eigen::MatrixXd y1(1, 2);
      y1(0) = y[0][0];
      y1(1) = y[0][1];
      return y1;
    };
    auto f2 = [&] (const std::vector<var>& x) {
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
  ts.resize(1); ts[0] = 1.0;
  std::vector<double> t0_vec{t0};
  
  {
    auto f1 = [&] (const std::vector<double>& x) {
      auto y = torsten::pmx_integrate_ode_bdf(f, y0, x[0], ts, theta , x_r, x_i);
      Eigen::MatrixXd y1(1, 3);
      y1(0) = y[0][0];
      y1(1) = y[0][1];
      y1(2) = y[0][2];
      return y1;
    };
    auto f2 = [&] (const std::vector<var>& x) {
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
  ts.resize(1); ts[0] = 1.0;
  std::vector<double> t0_vec{t0};
  
  {
    auto f1 = [&] (const std::vector<double>& x) {
      auto y = torsten::pmx_integrate_ode_bdf(f, y0, x[0], ts, theta , x_r, x_i);
      Eigen::MatrixXd y1(1, 3);
      y1(0) = y[0][0];
      y1(1) = y[0][1];
      y1(2) = y[0][2];
      return y1;
    };
    auto f2 = [&] (const std::vector<var>& x) {
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

TEST_F(TorstenOdeTest_sho, cvodes_ivp_system) {
  using Ode = PMXOdeSystem<F, double, double, double>;
  Ode ode{f, t0, ts, y0, theta, x_r, x_i, msgs};

  {
    PMXCvodesIntegrator<CV_BDF, CV_STAGGERED> solver(rtol, atol, 1000);
    OdeObserver<Ode> observer(ode);
    OdeDataObserver<Ode> observer_mat(ode);
    solver.integrate(ode, observer);
    solver.integrate(ode, observer_mat);
    std::vector<std::vector<double> > y1 =
      stan::math::integrate_ode_bdf(f, y0, t0, ts, theta , x_r, x_i);
    torsten::test::test_val(observer.y, y1);
    torsten::test::test_val(observer_mat.y, y1);

    auto y = torsten::pmx_integrate_ode_bdf(f, y0, t0, ts, theta , x_r, x_i);
    torsten::test::test_val(observer.y, y1);
  }

  {
    PMXCvodesIntegrator<CV_ADAMS, CV_STAGGERED> solver(rtol, atol, 1000);
    OdeObserver<Ode> observer(ode);
    OdeDataObserver<Ode> observer_mat(ode);
    solver.integrate(ode, observer);
    solver.integrate(ode, observer_mat);
    std::vector<std::vector<double> > y1 =
      stan::math::integrate_ode_adams(f, y0, t0, ts, theta , x_r, x_i);
    torsten::test::test_val(observer.y, y1);
    torsten::test::test_val(observer_mat.y, y1);

    auto y = torsten::pmx_integrate_ode_adams(f, y0, t0, ts, theta , x_r, x_i);
    torsten::test::test_val(observer.y, y1);
  }
}

TEST_F(TorstenOdeTest_lorenz, cvodes_ivp_system) {
  using Ode = PMXOdeSystem<F, double, double, double>;
  Ode ode{f, t0, ts, y0, theta, x_r, x_i, msgs};

  {
    PMXCvodesIntegrator<CV_BDF, CV_STAGGERED> solver(rtol, atol, 10000);
    OdeObserver<Ode> observer(ode);
    OdeDataObserver<Ode> observer_mat(ode);
    solver.integrate(ode, observer);
    solver.integrate(ode, observer_mat);
    std::vector<std::vector<double> > y1 =
      stan::math::integrate_ode_bdf(f, y0, t0, ts, theta , x_r, x_i);
    torsten::test::test_val(observer.y, y1);
    torsten::test::test_val(observer_mat.y, y1);

    auto y = torsten::pmx_integrate_ode_bdf(f, y0, t0, ts, theta , x_r, x_i);
    torsten::test::test_val(observer.y, y1);
  }

  {
    PMXCvodesIntegrator<CV_ADAMS, CV_STAGGERED> solver(rtol, atol, 10000);
    OdeObserver<Ode> observer(ode);
    OdeDataObserver<Ode> observer_mat(ode);
    solver.integrate(ode, observer);
    solver.integrate(ode, observer_mat);
    std::vector<std::vector<double> > y1 =
      stan::math::integrate_ode_adams(f, y0, t0, ts, theta , x_r, x_i);
    torsten::test::test_val(observer.y, y1);
    torsten::test::test_val(observer_mat.y, y1);

    auto y = torsten::pmx_integrate_ode_adams(f, y0, t0, ts, theta , x_r, x_i);
    torsten::test::test_val(observer.y, y1);
  }
}

TEST_F(TorstenOdeTest_chem, fwd_sensitivity_theta) {
  using Ode = PMXOdeSystem<F, double, double, var>;
  std::vector<var> theta_var = stan::math::to_var(theta);
  Ode ode(f, t0, ts, y0, theta_var, x_r, x_i, msgs);

  {
    PMXCvodesIntegrator<CV_ADAMS, CV_STAGGERED> solver(rtol, atol, 1e8);
    auto y1 = stan::math::integrate_ode_adams(f, y0, t0, ts, theta_var, x_r, x_i);
    OdeObserver<Ode> observer(ode);
    OdeDataObserver<Ode> observer_mat(ode);
    solver.integrate(ode, observer);
    solver.integrate(ode, observer_mat);
    torsten::test::test_grad(theta_var, y1, observer.y, 1.E-8, 1.E-5);
    torsten::test::test_grad(theta_var, y1, observer_mat.y, 1.E-8, 1.E-5);
  }

  {
    PMXCvodesIntegrator<CV_BDF, CV_STAGGERED> solver(rtol, atol, 1e8);
    auto y1 = stan::math::integrate_ode_bdf(f, y0, t0, ts, theta_var, x_r, x_i);
    OdeObserver<Ode> observer(ode);
    OdeDataObserver<Ode> observer_mat(ode);
    solver.integrate(ode, observer);
    solver.integrate(ode, observer_mat);
    torsten::test::test_grad(theta_var, y1, observer.y, 1.E-8, 1.E-5);
    torsten::test::test_grad(theta_var, y1, observer_mat.y, 1.E-8, 1.E-5);
  }
}

TEST_F(TorstenOdeTest_chem, fwd_sensitivity_y0) {
  using Ode = PMXOdeSystem<F, double, var, double>;
  std::vector<var> y0_var = stan::math::to_var(y0);
  Ode ode(f, t0, ts, y0_var, theta, x_r, x_i, msgs);

  {
    PMXCvodesIntegrator<CV_ADAMS, CV_STAGGERED> solver(rtol, atol, 1e8);
    auto y1 = stan::math::integrate_ode_adams(f, y0_var, t0, ts, theta, x_r, x_i);
    OdeObserver<Ode> observer(ode);
    OdeDataObserver<Ode> observer_mat(ode);
    solver.integrate(ode, observer);
    solver.integrate(ode, observer_mat);
    torsten::test::test_grad(y0_var, y1, observer.y, 1.E-8, 1.E-5);
    torsten::test::test_grad(y0_var, y1, observer_mat.y, 1.E-8, 1.E-5);
  }

  {
    PMXCvodesIntegrator<CV_BDF, CV_STAGGERED> solver(rtol, atol, 1e8);
    auto y1 = stan::math::integrate_ode_bdf(f, y0_var, t0, ts, theta, x_r, x_i);
    OdeObserver<Ode> observer(ode);
    OdeDataObserver<Ode> observer_mat(ode);
    solver.integrate(ode, observer);
    solver.integrate(ode, observer_mat);
    torsten::test::test_grad(y0_var, y1, observer.y, 1.E-8, 1.E-5);
    torsten::test::test_grad(y0_var, y1, observer_mat.y, 1.E-8, 1.E-5);
  }
}

TEST_F(TorstenOdeTest_chem, fwd_sensitivity_theta_y0) {
  using Ode = PMXOdeSystem<F, double, var, var>;
  std::vector<var> theta_var = stan::math::to_var(theta);
  std::vector<var> y0_var = stan::math::to_var(y0);
  Ode ode(f, t0, ts, y0_var, theta_var, x_r, x_i, msgs);

  {
    PMXCvodesIntegrator<CV_ADAMS, CV_STAGGERED> solver(rtol, atol, 1e8);
    auto y1 = stan::math::integrate_ode_adams(f, y0_var, t0, ts, theta_var, x_r, x_i);
    OdeObserver<Ode> observer(ode);
    OdeDataObserver<Ode> observer_mat(ode);
    solver.integrate(ode, observer);
    solver.integrate(ode, observer_mat);    
    torsten::test::test_grad(y0_var, observer.y, y1, 1.E-8, 1.E-5);
    torsten::test::test_grad(theta_var, observer.y, y1, 1.E-8, 1.E-5);

    std::vector<stan::math::var> vars(y0_var);
    vars.insert(vars.end(), theta_var.begin(), theta_var.end());
    torsten::test::test_grad(vars, y1, observer_mat.y, 1.E-8, 1.E-5);
  }

  {
    PMXCvodesIntegrator<CV_BDF, CV_STAGGERED> solver(rtol, atol, 1e8);
    auto y1 = stan::math::integrate_ode_bdf(f, y0_var, t0, ts, theta_var, x_r, x_i);
    OdeObserver<Ode> observer(ode);
    OdeDataObserver<Ode> observer_mat(ode);
    solver.integrate(ode, observer);
    solver.integrate(ode, observer_mat);    
    torsten::test::test_grad(y0_var, observer.y, y1, 1.E-8, 1.E-5);
    torsten::test::test_grad(theta_var, observer.y, y1, 1.E-8, 1.E-5);

    std::vector<stan::math::var> vars(y0_var);
    vars.insert(vars.end(), theta_var.begin(), theta_var.end());
    torsten::test::test_grad(vars, y1, observer_mat.y, 1.E-8, 1.E-5);
  }
}

TEST_F(TorstenOdeTest_sho, fwd_sensitivity_theta) {
  using Ode = PMXOdeSystem<F, double, double, var>;
  std::vector<var> theta_var = stan::math::to_var(theta);
  Ode ode(f, t0, ts, y0, theta_var, x_r, x_i, msgs);

  {
    PMXCvodesIntegrator<CV_ADAMS, CV_STAGGERED> solver(rtol, atol, 1e8);
    auto y1 = stan::math::integrate_ode_adams(f, y0, t0, ts, theta_var, x_r, x_i);
    OdeObserver<Ode> observer(ode);
    OdeDataObserver<Ode> observer_mat(ode);
    solver.integrate(ode, observer);
    solver.integrate(ode, observer_mat);
    torsten::test::test_grad(theta_var, y1, observer.y, 1.E-8, 1.E-5);
    torsten::test::test_grad(theta_var, y1, observer_mat.y, 1.E-8, 1.E-5);
  }

  {
    PMXCvodesIntegrator<CV_BDF, CV_STAGGERED> solver(rtol, atol, 1e8);
    auto y1 = stan::math::integrate_ode_bdf(f, y0, t0, ts, theta_var, x_r, x_i);
    OdeObserver<Ode> observer(ode);
    OdeDataObserver<Ode> observer_mat(ode);
    solver.integrate(ode, observer);
    solver.integrate(ode, observer_mat);
    torsten::test::test_grad(theta_var, y1, observer.y, 1.E-8, 1.E-5);
    torsten::test::test_grad(theta_var, y1, observer_mat.y, 1.E-8, 1.E-5);
  }
}

TEST_F(TorstenOdeTest_sho, fwd_sensitivity_y0) {
  using Ode = PMXOdeSystem<F, double, var, double>;
  std::vector<var> y0_var = stan::math::to_var(y0);
  Ode ode(f, t0, ts, y0_var, theta, x_r, x_i, msgs);

  {
    PMXCvodesIntegrator<CV_ADAMS, CV_STAGGERED> solver(rtol, atol, 1e8);
    auto y1 = stan::math::integrate_ode_adams(f, y0_var, t0, ts, theta, x_r, x_i);
    OdeObserver<Ode> observer(ode);
    OdeDataObserver<Ode> observer_mat(ode);
    solver.integrate(ode, observer);
    solver.integrate(ode, observer_mat);
    torsten::test::test_grad(y0_var, y1, observer.y, 1.E-8, 1.E-5);
    torsten::test::test_grad(y0_var, y1, observer_mat.y, 1.E-8, 1.E-5);
  }

  {
    PMXCvodesIntegrator<CV_BDF, CV_STAGGERED> solver(rtol, atol, 1e8);
    auto y1 = stan::math::integrate_ode_bdf(f, y0_var, t0, ts, theta, x_r, x_i);
    OdeObserver<Ode> observer(ode);
    OdeDataObserver<Ode> observer_mat(ode);
    solver.integrate(ode, observer);
    solver.integrate(ode, observer_mat);
    torsten::test::test_grad(y0_var, y1, observer.y, 1.E-8, 1.E-5);
    torsten::test::test_grad(y0_var, y1, observer_mat.y, 1.E-8, 1.E-5);
  }
}

TEST_F(TorstenOdeTest_sho, fwd_sensitivity_ts) {
  using stan::math::value_of;
  using Ode = PMXOdeSystem<F, var, double, double>;
  PMXCvodesIntegrator<CV_ADAMS, CV_STAGGERED> solver(rtol, atol, 1e8);
  std::vector<var> ts_var = stan::math::to_var(ts);
  Ode ode(f, t0, ts_var, y0, theta, x_r, x_i, msgs);
  OdeObserver<Ode> observer(ode);
  OdeDataObserver<Ode> observer_mat(ode);
  solver.integrate(ode, observer);

  std::vector<std::vector<var>>& y = observer.y;

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
  using Ode = PMXOdeSystem<F, var, double, double>;
  PMXCvodesIntegrator<CV_ADAMS, CV_STAGGERED> solver(rtol, atol, 1e8);
  std::vector<var> ts_var = stan::math::to_var(ts);
  Ode ode(f, t0, ts_var, y0, theta, x_r, x_i, msgs);
  OdeObserver<Ode> observer(ode);
  OdeDataObserver<Ode> observer_mat(ode);
  solver.integrate(ode, observer);

  std::vector<std::vector<var>>& y = observer.y;

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
  using Ode = PMXOdeSystem<F, double, var, var>;
  std::vector<var> theta_var = stan::math::to_var(theta);
  std::vector<var> y0_var = stan::math::to_var(y0);
  Ode ode(f, t0, ts, y0_var, theta_var, x_r, x_i, msgs);

  {
    PMXCvodesIntegrator<CV_ADAMS, CV_STAGGERED> solver(rtol, atol, 1e8);
    auto y1 = stan::math::integrate_ode_adams(f, y0_var, t0, ts, theta_var, x_r, x_i);
    OdeObserver<Ode> observer(ode);
    OdeDataObserver<Ode> observer_mat(ode);
    solver.integrate(ode, observer);
    solver.integrate(ode, observer_mat);    
    torsten::test::test_grad(y0_var, observer.y, y1, 1.E-8, 1.E-5);
    torsten::test::test_grad(theta_var, observer.y, y1, 1.E-8, 1.E-5);

    std::vector<stan::math::var> vars(y0_var);
    vars.insert(vars.end(), theta_var.begin(), theta_var.end());
    torsten::test::test_grad(vars, y1, observer_mat.y, 1.E-8, 1.E-5);
  }

  {
    PMXCvodesIntegrator<CV_BDF, CV_STAGGERED> solver(rtol, atol, 1e8);
    auto y1 = stan::math::integrate_ode_bdf(f, y0_var, t0, ts, theta_var, x_r, x_i);
    OdeObserver<Ode> observer(ode);
    OdeDataObserver<Ode> observer_mat(ode);
    solver.integrate(ode, observer);
    solver.integrate(ode, observer_mat);    
    torsten::test::test_grad(y0_var, observer.y, y1, 1.E-8, 1.E-5);
    torsten::test::test_grad(theta_var, observer.y, y1, 1.E-8, 1.E-5);

    std::vector<stan::math::var> vars(y0_var);
    vars.insert(vars.end(), theta_var.begin(), theta_var.end());
    torsten::test::test_grad(vars, y1, observer_mat.y, 1.E-8, 1.E-5);
  }
}

TEST_F(TorstenOdeTest_lorenz, fwd_sensitivity_theta_y0_ts) {
  using stan::math::value_of;
  using Ode = PMXOdeSystem<F, var, var, var>;
  PMXCvodesIntegrator<CV_ADAMS, CV_STAGGERED> solver(rtol, atol, 1e8);
  std::vector<var> theta_var = stan::math::to_var(theta);
  std::vector<var> y0_var = stan::math::to_var(y0);
  std::vector<var> ts_var = stan::math::to_var(ts);
  Ode ode(f, t0, ts_var, y0_var, theta_var, x_r, x_i, msgs);
  OdeObserver<Ode> observer(ode);
  solver.integrate(ode, observer);
  std::vector<std::vector<var> >& y = observer.y;
  std::vector<std::vector<var> > y1 = stan::math::integrate_ode_adams(f, y0_var, t0, ts, theta_var, x_r, x_i);

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
  using Ode = PMXOdeSystem<F, var, var, var>;
  PMXCvodesIntegrator<CV_BDF, CV_STAGGERED> solver(rtol, atol, 1e8);
  std::vector<var> theta_var = stan::math::to_var(theta);
  std::vector<var> y0_var = stan::math::to_var(y0);
  std::vector<var> ts_var = stan::math::to_var(ts);
  Ode ode(f, t0, ts_var, y0_var, theta_var, x_r, x_i, msgs);
  OdeObserver<Ode> observer(ode);
  solver.integrate(ode, observer);
  std::vector<std::vector<var> >& y = observer.y;
  std::vector<std::vector<var> > y1 = stan::math::integrate_ode_bdf(f, y0_var, t0, ts, theta_var, x_r, x_i);

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
  using Ode = PMXOdeSystem<F, var, var, double>;
  PMXCvodesIntegrator<CV_BDF, CV_STAGGERED> solver(rtol, atol, 1e8);
  std::vector<var> y0_var = stan::math::to_var(y0);
  std::vector<var> ts_var = stan::math::to_var(ts);
  Ode ode(f, t0, ts_var, y0_var, theta, x_r, x_i, msgs);
  OdeObserver<Ode> observer(ode);
  solver.integrate(ode, observer);
  std::vector<std::vector<var> >& y = observer.y;
  std::vector<std::vector<var> > y1 = stan::math::integrate_ode_bdf(f, y0_var, t0, ts_var, theta, x_r, x_i);

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
  using Ode = PMXOdeSystem<F, var, double, var>;
  PMXCvodesIntegrator<CV_BDF, CV_STAGGERED> solver(rtol, atol, 1e8);
  std::vector<var> theta_var = stan::math::to_var(theta);
  std::vector<var> ts_var = stan::math::to_var(ts);
  Ode ode(f, t0, ts_var, y0, theta_var, x_r, x_i, msgs);
  OdeObserver<Ode> observer(ode);
  solver.integrate(ode, observer);
  std::vector<std::vector<var> >& y = observer.y;
  std::vector<std::vector<var> > y1 = stan::math::integrate_ode_bdf(f, y0, t0, ts_var, theta_var, x_r, x_i);

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
  using Ode = PMXOdeSystem<F, var, double, var>;
  PMXCvodesIntegrator<CV_ADAMS, CV_STAGGERED> solver(rtol, atol, 1e8);
  std::vector<var> theta_var = stan::math::to_var(theta);
  std::vector<var> ts_var = stan::math::to_var(ts);
  Ode ode(f, t0, ts_var, y0, theta_var, x_r, x_i, msgs);
  OdeObserver<Ode> observer(ode);
  solver.integrate(ode, observer);
  std::vector<std::vector<var> >& y = observer.y;
  std::vector<std::vector<var> > y1 = stan::math::integrate_ode_adams(f, y0, t0, ts_var, theta_var, x_r, x_i);

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
  for (int i = 0; i < 3; ++i) {
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
  for (int i = 0; i < 3; ++i) {
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
  for (int i = 0; i < 3; ++i) {
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
  using torsten::dsolve::PMXOdeSystem;
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
//   using torsten::dsolve::PMXOdeSystem;
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

//   PMXOdeSystem<harm_osc_ode_fun, double, double, double> ode{
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
