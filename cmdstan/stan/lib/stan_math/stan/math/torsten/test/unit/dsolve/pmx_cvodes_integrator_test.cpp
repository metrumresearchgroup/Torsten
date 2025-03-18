#include <stan/math.hpp>
#include <stan/math/rev/core.hpp>
#include <test/unit/math/rev/fun/util.hpp>
#include <stan/math/torsten/dsolve/pmx_ode_adams.hpp>
#include <stan/math/torsten/dsolve/pmx_ode_bdf.hpp>
#include <stan/math/torsten/dsolve/pmx_ode_bdf.hpp>
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

using torsten::dsolve::PMXVariadicOdeSystem;
using torsten::dsolve::PMXCvodesIntegrator;
using torsten::dsolve::PMXOdeService;
using torsten::dsolve::OdeObserver;
using torsten::dsolve::OdeDataObserver;
using torsten::PMXCvodesSensMethod;
using stan::math::nested_rev_autodiff;
using torsten::AD;
using torsten::CSDA;
using torsten::DQ;
using stan::math::var;

TEST_F(TorstenOdeTest_sho, t0_var) {
  ts.resize(1); ts[0] = 1.0;
  std::vector<double> t0_vec{t0};

  {
    auto f1 = [&] (const std::vector<double>& x) {
      auto y = torsten::pmx_ode_bdf(f_eigen, y0_vec, x[0], ts, msgs, theta, x_r, x_i);
      Eigen::MatrixXd y1(1, 2);
      y1(0) = y[0][0];
      y1(1) = y[0][1];
      return y1;
    };
    auto f2 = [&] (const std::vector<var>& x) {
      double t0 = stan::math::value_of(x[0]);
      std::vector<var> ts_v{t0 + ts[0] - x[0]};
      auto y = torsten::pmx_ode_bdf(f_eigen, y0_vec, t0, ts_v, msgs, theta, x_r, x_i);
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
      auto y = torsten::pmx_ode_bdf(f_eigen, y0_vec, x[0], ts, msgs, theta , x_r, x_i);
      Eigen::MatrixXd y1(1, 3);
      y1(0) = y[0][0];
      y1(1) = y[0][1];
      y1(2) = y[0][2];
      return y1;
    };
    auto f2 = [&] (const std::vector<var>& x) {
      double t0 = stan::math::value_of(x[0]);
      std::vector<var> ts_v{t0 + ts[0] - x[0]};
      auto y = torsten::pmx_ode_bdf(f_eigen, y0_vec, t0, ts_v, msgs, theta , x_r, x_i);
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
      auto y = torsten::pmx_ode_bdf(f_eigen, y0_vec, x[0], ts, msgs, theta , x_r, x_i);
      Eigen::MatrixXd y1(1, 3);
      y1(0) = y[0][0];
      y1(1) = y[0][1];
      y1(2) = y[0][2];
      return y1;
    };
    auto f2 = [&] (const std::vector<var>& x) {
      double t0 = stan::math::value_of(x[0]);
      std::vector<var> ts_v{t0 + ts[0] - x[0]};
      auto y = torsten::pmx_ode_bdf(f_eigen, y0_vec, t0, ts_v, msgs, theta , x_r, x_i);
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
  using Ode = PMXVariadicOdeSystem<F_eigen, double, double, std::vector<double>, std::vector<double>, std::vector<int> >;
  Ode ode{f_eigen, t0, ts, y0_vec, msgs, theta, x_r, x_i};

  {
    PMXCvodesIntegrator<CV_BDF, CV_STAGGERED> solver(rtol, atol, 1000);
    OdeObserver<Ode> observer(ode);
    OdeDataObserver<Ode> observer_mat(ode);
    solver.integrate(ode, observer);
    solver.integrate(ode, observer_mat);
    auto y1 = stan::math::ode_bdf(f_eigen, y0_vec, t0, ts, nullptr, theta , x_r, x_i);
    EXPECT_ARRAY2D_VAL_FLOAT_EQ(observer.y, y1);
    EXPECT_MAT_ARRAY2D_VAL_FLOAT_EQ(observer_mat.y, y1);

    auto y = torsten::pmx_ode_bdf(f_eigen, y0_vec, t0, ts, msgs, theta , x_r, x_i);
    EXPECT_ARRAY2D_VAL_FLOAT_EQ(observer.y, y1);
  }

  {
    PMXCvodesIntegrator<CV_ADAMS, CV_STAGGERED> solver(rtol, atol, 1000);
    OdeObserver<Ode> observer(ode);
    OdeDataObserver<Ode> observer_mat(ode);
    solver.integrate(ode, observer);
    solver.integrate(ode, observer_mat);
    auto y1 = stan::math::ode_adams(f_eigen, y0_vec, t0, ts, nullptr, theta , x_r, x_i);
    EXPECT_ARRAY2D_VAL_FLOAT_EQ(observer.y, y1);
    EXPECT_MAT_ARRAY2D_VAL_FLOAT_EQ(observer_mat.y, y1);

    auto y = torsten::pmx_ode_adams(f_eigen, y0_vec, t0, ts, msgs, theta , x_r, x_i);
    EXPECT_ARRAY2D_VAL_FLOAT_EQ(observer.y, y1);
  }
}

TEST_F(TorstenOdeTest_lorenz, cvodes_ivp_system) {
  using Ode = PMXVariadicOdeSystem<F_eigen, double, double, std::vector<double>, std::vector<double>, std::vector<int> >;
  Ode ode{f_eigen, t0, ts, y0_vec, msgs, theta, x_r, x_i};

  {
    PMXCvodesIntegrator<CV_BDF, CV_STAGGERED> solver(rtol, atol, 10000);
    OdeObserver<Ode> observer(ode);
    OdeDataObserver<Ode> observer_mat(ode);
    solver.integrate(ode, observer);
    solver.integrate(ode, observer_mat);
    auto y1 = stan::math::ode_bdf(f_eigen, y0_vec, t0, ts, nullptr, theta , x_r, x_i);
    EXPECT_ARRAY2D_VAL_FLOAT_EQ(observer.y, y1);
    EXPECT_MAT_ARRAY2D_VAL_FLOAT_EQ(observer_mat.y, y1);

    auto y = torsten::pmx_ode_bdf(f_eigen, y0_vec, t0, ts, msgs, theta , x_r, x_i);
    EXPECT_ARRAY2D_VAL_FLOAT_EQ(observer.y, y1);
  }

  {
    PMXCvodesIntegrator<CV_ADAMS, CV_STAGGERED> solver(rtol, atol, 10000);
    OdeObserver<Ode> observer(ode);
    OdeDataObserver<Ode> observer_mat(ode);
    solver.integrate(ode, observer);
    solver.integrate(ode, observer_mat);
    auto y1 = stan::math::ode_adams(f_eigen, y0_vec, t0, ts, nullptr, theta , x_r, x_i);
    EXPECT_ARRAY2D_VAL_FLOAT_EQ(observer.y, y1);
    EXPECT_MAT_ARRAY2D_VAL_FLOAT_EQ(observer_mat.y, y1);

    auto y = torsten::pmx_ode_adams(f_eigen, y0_vec, t0, ts, msgs, theta , x_r, x_i);
    EXPECT_ARRAY2D_VAL_FLOAT_EQ(observer.y, y1);
  }
}

TEST_F(TorstenOdeTest_chem, fwd_sensitivity_theta) {
  using Ode = PMXVariadicOdeSystem<F_eigen, double, double, std::vector<var>, std::vector<double>, std::vector<int> >;

  nested_rev_autodiff nested;

  std::vector<var> theta_var = stan::math::to_var(theta);
  Ode ode(f_eigen, t0, ts, y0_vec, msgs, theta_var, x_r, x_i);

  {
    PMXCvodesIntegrator<CV_ADAMS, CV_STAGGERED> solver(rtol, atol, 1e8);
    auto y1 = stan::math::ode_adams(f_eigen, y0_vec, t0, ts, nullptr, theta_var, x_r, x_i);
    OdeObserver<Ode> observer(ode);
    OdeDataObserver<Ode> observer_mat(ode);
    solver.integrate(ode, observer);
    solver.integrate(ode, observer_mat);
    EXPECT_ARRAY2D_ADJ_NEAR(y1, observer.y, theta_var, nested, 1.E-5, "THETA");

    auto y2 = torsten::precomputed_gradients(observer_mat.y, theta_var);
    EXPECT_MAT_ARRAY2D_ADJ_NEAR(y2, y1, theta_var, nested, 1.E-5, "THETA");
  }

  {
    PMXCvodesIntegrator<CV_BDF, CV_STAGGERED> solver(rtol, atol, 1e8);
    auto y1 = stan::math::ode_bdf(f_eigen, y0_vec, t0, ts, nullptr, theta_var, x_r, x_i);
    OdeObserver<Ode> observer(ode);
    OdeDataObserver<Ode> observer_mat(ode);
    solver.integrate(ode, observer);
    solver.integrate(ode, observer_mat);
    EXPECT_ARRAY2D_ADJ_NEAR(y1, observer.y, theta_var, nested, 1.E-5, "THETA");

    auto y2 = torsten::precomputed_gradients(observer_mat.y, theta_var);
    EXPECT_MAT_ARRAY2D_ADJ_NEAR(y2, y1, theta_var, nested, 1.E-5, "THETA");
  }
}

TEST_F(TorstenOdeTest_chem, fwd_sensitivity_y0) {
  using Ode = PMXVariadicOdeSystem<F_eigen, double, var, std::vector<double>, std::vector<double>, std::vector<int> >;

  nested_rev_autodiff nested;

  Eigen::Matrix<var, -1, 1> y0_vec_var = stan::math::to_var(y0_vec);
  Ode ode(f_eigen, t0, ts, y0_vec_var, msgs, theta, x_r, x_i);

  {
    PMXCvodesIntegrator<CV_ADAMS, CV_STAGGERED> solver(rtol, atol, 1e8);
    auto y1 = stan::math::ode_adams(f_eigen, y0_vec_var, t0, ts, nullptr, theta, x_r, x_i);
    OdeObserver<Ode> observer(ode);
    OdeDataObserver<Ode> observer_mat(ode);
    solver.integrate(ode, observer);
    solver.integrate(ode, observer_mat);

    auto y0_var = stan::math::to_array_1d(y0_vec_var);
    EXPECT_ARRAY2D_ADJ_NEAR(y1, observer.y, y0_vec_var, nested, 1.E-5, "y0");
    auto y2 = torsten::precomputed_gradients(observer_mat.y, y0_vec_var);
    EXPECT_MAT_ARRAY2D_ADJ_NEAR(y2, y1, y0_vec_var, nested, 1.E-5, "y0");
  }

  {
    PMXCvodesIntegrator<CV_BDF, CV_STAGGERED> solver(rtol, atol, 1e8);
    auto y1 = stan::math::ode_bdf(f_eigen, y0_vec_var, t0, ts, nullptr, theta, x_r, x_i);
    OdeObserver<Ode> observer(ode);
    OdeDataObserver<Ode> observer_mat(ode);
    solver.integrate(ode, observer);
    solver.integrate(ode, observer_mat);

    auto y0_var = stan::math::to_array_1d(y0_vec_var);
    EXPECT_ARRAY2D_ADJ_NEAR(y1, observer.y, y0_vec_var, nested, 1.E-5, "y0");

    auto y2 = torsten::precomputed_gradients(observer_mat.y, y0_vec_var);
    EXPECT_MAT_ARRAY2D_ADJ_NEAR(y2, y1, y0_vec_var, nested, 1.E-5, "y0");
  }
}

// TEST_F(TorstenOdeTest_chem, fwd_sensitivity_theta_y0) {
//   using Ode = PMXVariadicOdeSystem<F_eigen, double, var, std::vector<var>, std::vector<double>, std::vector<int> >;

//   nested_rev_autodiff nested;

//   std::vector<var> theta_var = stan::math::to_var(theta);
//   Eigen::Matrix<var, -1, 1> y0_vec_var = stan::math::to_var(y0_vec);
//   Ode ode(f_eigen, t0, ts, y0_vec_var, msgs, theta_var, x_r, x_i);

//   {
//     PMXCvodesIntegrator<CV_ADAMS, CV_STAGGERED> solver(rtol, atol, 1e8);
//     auto y1 = stan::math::ode_adams(f_eigen, y0_vec_var, t0, ts, nullptr, theta_var, x_r, x_i);
//     OdeObserver<Ode> observer(ode);
//     OdeDataObserver<Ode> observer_mat(ode);
//     solver.integrate(ode, observer);
//     solver.integrate(ode, observer_mat);
//     auto y0_var = stan::math::to_array_1d(y0_vec_var);
//     EXPECT_ARRAY2D_ADJ_NEAR(observer.y, y1, y0_vec_var, nested, 1.E-5, "y0");
//     EXPECT_ARRAY2D_ADJ_NEAR(observer.y, y1, theta_var, nested, 1.E-5, "THETA");

//     std::vector<stan::math::var> vars(stan::math::to_array_1d(y0_vec_var));
//     vars.insert(vars.end(), theta_var.begin(), theta_var.end());

//     auto y2 = torsten::precomputed_gradients(observer_mat.y, vars);
//     EXPECT_MAT_ARRAY2D_ADJ_NEAR(y2, y1, vars, nested, 1.E-5, "y0");
//   }

//   {
//     PMXCvodesIntegrator<CV_BDF, CV_STAGGERED> solver(rtol, atol, 1e8);
//     auto y1 = stan::math::ode_bdf(f_eigen, y0_vec_var, t0, ts, nullptr, theta_var, x_r, x_i);
//     OdeObserver<Ode> observer(ode);
//     OdeDataObserver<Ode> observer_mat(ode);
//     solver.integrate(ode, observer);
//     solver.integrate(ode, observer_mat);
//     auto y0_var = stan::math::to_array_1d(y0_vec_var);
//     EXPECT_ARRAY2D_ADJ_NEAR(observer.y, y1, y0_vec_var, nested, 1.E-5, "y0");
//     EXPECT_ARRAY2D_ADJ_NEAR(observer.y, y1, theta_var, nested, 1.E-5, "THETA");

//     std::vector<stan::math::var> vars(stan::math::to_array_1d(y0_vec_var));
//     vars.insert(vars.end(), theta_var.begin(), theta_var.end());

//     auto y2 = torsten::precomputed_gradients(observer_mat.y, vars);
//     EXPECT_MAT_ARRAY2D_ADJ_NEAR(y2, y1, vars, nested, 1.E-5, "y0");
//   }
// }

// TEST_F(TorstenOdeTest_sho, fwd_sensitivity_theta) {
//   using Ode = PMXVariadicOdeSystem<F_eigen, double, double, std::vector<var>, std::vector<double>, std::vector<int> >;

//   nested_rev_autodiff nested;

//   std::vector<var> theta_var = stan::math::to_var(theta);
//   Ode ode(f_eigen, t0, ts, y0_vec, msgs, theta_var, x_r, x_i);

//   {
//     PMXCvodesIntegrator<CV_ADAMS, CV_STAGGERED> solver(rtol, atol, 1e8);
//     auto y1 = stan::math::ode_adams(f_eigen, y0_vec, t0, ts, nullptr, theta_var, x_r, x_i);
//     OdeObserver<Ode> observer(ode);
//     OdeDataObserver<Ode> observer_mat(ode);
//     solver.integrate(ode, observer);
//     solver.integrate(ode, observer_mat);
//     EXPECT_ARRAY2D_ADJ_NEAR(y1, observer.y, theta_var, nested, 1.E-5, "THETA");

//     auto y2 = torsten::precomputed_gradients(observer_mat.y, theta_var);
//     EXPECT_MAT_ARRAY2D_ADJ_NEAR(y2, y1, theta_var, nested, 1.E-5, "y0");
//   }

//   {
//     PMXCvodesIntegrator<CV_BDF, CV_STAGGERED> solver(rtol, atol, 1e8);
//     auto y1 = stan::math::ode_bdf(f_eigen, y0_vec, t0, ts, nullptr, theta_var, x_r, x_i);
//     OdeObserver<Ode> observer(ode);
//     OdeDataObserver<Ode> observer_mat(ode);
//     solver.integrate(ode, observer);
//     solver.integrate(ode, observer_mat);
//     EXPECT_ARRAY2D_ADJ_NEAR(y1, observer.y, theta_var, nested, 1.E-5, "THETA");

//     auto y2 = torsten::precomputed_gradients(observer_mat.y, theta_var);
//     EXPECT_MAT_ARRAY2D_ADJ_NEAR(y2, y1, theta_var, nested, 1.E-5, "y0");
//   }
// }

// TEST_F(TorstenOdeTest_sho, fwd_sensitivity_y0) {
//   using Ode = PMXVariadicOdeSystem<F_eigen, double, var, std::vector<double>, std::vector<double>, std::vector<int> >;

//   nested_rev_autodiff nested;

//   Eigen::Matrix<var, -1, 1> y0_vec_var = stan::math::to_var(y0_vec);
//   Ode ode(f_eigen, t0, ts, y0_vec_var, msgs, theta, x_r, x_i);

//   {
//     PMXCvodesIntegrator<CV_ADAMS, CV_STAGGERED> solver(rtol, atol, 1e8);
//     auto y1 = stan::math::ode_adams(f_eigen, y0_vec_var, t0, ts, nullptr, theta, x_r, x_i);
//     OdeObserver<Ode> observer(ode);
//     OdeDataObserver<Ode> observer_mat(ode);
//     solver.integrate(ode, observer);
//     solver.integrate(ode, observer_mat);
//     auto y0_var = stan::math::to_array_1d(y0_vec_var);
//     EXPECT_ARRAY2D_ADJ_NEAR(y1, observer.y, y0_vec_var, nested, 1.E-5, "y0");

//     auto y2 = torsten::precomputed_gradients(observer_mat.y, y0_vec_var);
//     EXPECT_MAT_ARRAY2D_ADJ_NEAR(y2, y1, y0_vec_var, nested, 1.E-5, "y0");
//   }

//   {
//     PMXCvodesIntegrator<CV_BDF, CV_STAGGERED> solver(rtol, atol, 1e8);
//     auto y1 = stan::math::ode_bdf(f_eigen, y0_vec_var, t0, ts, nullptr, theta, x_r, x_i);
//     OdeObserver<Ode> observer(ode);
//     OdeDataObserver<Ode> observer_mat(ode);
//     solver.integrate(ode, observer);
//     solver.integrate(ode, observer_mat);
//     auto y0_var = stan::math::to_array_1d(y0_vec_var);
//     EXPECT_ARRAY2D_ADJ_NEAR(y1, observer.y, y0_vec_var, nested, 1.E-5, "y0");

//     auto y2 = torsten::precomputed_gradients(observer_mat.y, y0_vec_var);
//     EXPECT_MAT_ARRAY2D_ADJ_NEAR(y2, y1, y0_vec_var, nested, 1.E-5, "y0");
//   }
// }

TEST_F(TorstenOdeTest_sho, fwd_sensitivity_ts) {
  using stan::math::value_of;
  using Ode = PMXVariadicOdeSystem<F_eigen, var, double, std::vector<double>, std::vector<double>, std::vector<int> >;
  PMXCvodesIntegrator<CV_ADAMS, CV_STAGGERED> solver(rtol, atol, 1e8);
  std::vector<var> ts_var = stan::math::to_var(ts);
  Ode ode(f_eigen, t0, ts_var, y0_vec, msgs, theta, x_r, x_i);
  OdeObserver<Ode> observer(ode);
  OdeDataObserver<Ode> observer_mat(ode);
  solver.integrate(ode, observer);

  auto& y = observer.y;
  std::vector<double> y_i(y[0].size());

  std::vector<double> g(ts.size()), fval(y0.size());
  for (size_t i = 0; i < ts.size(); ++i) {
    for (size_t m = 0; m < y_i.size(); ++m) {
      y_i[m] = y[i](m).val();
    }
    fval = f(value_of(ts[i]), y_i, theta, x_r, x_i, msgs);
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
  using Ode = PMXVariadicOdeSystem<F_eigen, var, double, std::vector<double>, std::vector<double>, std::vector<int> >;
  PMXCvodesIntegrator<CV_ADAMS, CV_STAGGERED> solver(rtol, atol, 1e8);
  std::vector<var> ts_var = stan::math::to_var(ts);
  Ode ode(f_eigen, t0, ts_var, y0_vec, msgs, theta, x_r, x_i);
  OdeObserver<Ode> observer(ode);
  OdeDataObserver<Ode> observer_mat(ode);
  solver.integrate(ode, observer);

  auto& y = observer.y;
  std::vector<double> y_i(y[0].size());

  std::vector<double> g(ts.size()), fval(y0.size());
  for (size_t i = 0; i < ts.size(); ++i) {
    for (size_t m = 0; m < y_i.size(); ++m) {
      y_i[m] = y[i](m).val();
    }
    fval = f(value_of(ts[i]), y_i, theta, x_r, x_i, msgs);
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
  using Ode = PMXVariadicOdeSystem<F_eigen, double, var, std::vector<var>, std::vector<double>, std::vector<int> >;

  nested_rev_autodiff nested;

  std::vector<var> theta_var = stan::math::to_var(theta);
  Eigen::Matrix<var, -1, 1> y0_vec_var = stan::math::to_var(y0_vec);
  Ode ode(f_eigen, t0, ts, y0_vec_var, msgs, theta_var, x_r, x_i);

  {
    PMXCvodesIntegrator<CV_ADAMS, CV_STAGGERED> solver(rtol, atol, 1e8);
    auto y1 = stan::math::ode_adams(f_eigen, y0_vec_var, t0, ts, nullptr, theta_var, x_r, x_i);
    OdeObserver<Ode> observer(ode);
    OdeDataObserver<Ode> observer_mat(ode);
    solver.integrate(ode, observer);
    solver.integrate(ode, observer_mat);
    auto y0_var = stan::math::to_array_1d(y0_vec_var);
    EXPECT_ARRAY2D_ADJ_NEAR(observer.y, y1, y0_vec_var, nested, 1.E-5, "y0");
    EXPECT_ARRAY2D_ADJ_NEAR(observer.y, y1, theta_var, nested, 1.E-5, "THETA");

    std::vector<stan::math::var> vars(stan::math::to_array_1d(y0_vec_var));
    vars.insert(vars.end(), theta_var.begin(), theta_var.end());

    auto y2 = torsten::precomputed_gradients(observer_mat.y, vars);
    EXPECT_MAT_ARRAY2D_ADJ_NEAR(y2, y1, vars, nested, 1.E-5, "vars");
  }

  {
    PMXCvodesIntegrator<CV_BDF, CV_STAGGERED> solver(rtol, atol, 1e8);
    auto y1 = stan::math::ode_bdf(f_eigen, y0_vec_var, t0, ts, nullptr, theta_var, x_r, x_i);
    OdeObserver<Ode> observer(ode);
    OdeDataObserver<Ode> observer_mat(ode);
    solver.integrate(ode, observer);
    solver.integrate(ode, observer_mat);
    auto y0_var = stan::math::to_array_1d(y0_vec_var);
    EXPECT_ARRAY2D_ADJ_NEAR(observer.y, y1, y0_vec_var, nested, 1.E-5, "y0");
    EXPECT_ARRAY2D_ADJ_NEAR(observer.y, y1, theta_var, nested, 1.E-5, "THETA");

    std::vector<stan::math::var> vars(stan::math::to_array_1d(y0_vec_var));
    vars.insert(vars.end(), theta_var.begin(), theta_var.end());

    auto y2 = torsten::precomputed_gradients(observer_mat.y, vars);
    EXPECT_MAT_ARRAY2D_ADJ_NEAR(y2, y1, vars, nested, 1.E-5, "vars");
  }
}

// TEST_F(TorstenOdeTest_lorenz, fwd_sensitivity_theta_y0_ts) {
//   using stan::math::value_of;

//   nested_rev_autodiff nested;

//   using Ode = PMXVariadicOdeSystem<F_eigen, var, var, std::vector<var>, std::vector<double>, std::vector<int> >;
//   PMXCvodesIntegrator<CV_ADAMS, CV_STAGGERED> solver(rtol, atol, 1e8);
//   std::vector<var> theta_var = stan::math::to_var(theta);
//   Eigen::Matrix<var, -1, 1> y0_vec_var = stan::math::to_var(y0_vec);
//   std::vector<var> ts_var = stan::math::to_var(ts);
//   Ode ode(f_eigen, t0, ts_var, y0_vec_var, msgs, theta_var, x_r, x_i);
//   OdeObserver<Ode> observer(ode);
//   solver.integrate(ode, observer);
//   auto& y = observer.y;
//   auto y1 = stan::math::ode_adams(f_eigen, y0_vec_var, t0, ts, nullptr, theta_var, x_r, x_i);

//   EXPECT_ARRAY2D_ADJ_NEAR(y , y1, y0_vec_var   , nested, 1.E-5, "y0");
//   EXPECT_ARRAY2D_ADJ_NEAR(y , y1, theta_var, nested, 1.E-4, "THETA");

//   std::vector<double> y_i(y[0].size());

//   std::vector<double> g(ts.size()), fval(y0.size());
//   for (size_t i = 0; i < ts.size(); ++i) {
//     for (size_t m = 0; m < y_i.size(); ++m) {
//       y_i[m] = y[i](m).val();
//     }
//     fval = f(value_of(ts[i]), y_i, theta, x_r, x_i, msgs);
//     for (size_t j = 0; j < y0.size(); ++j) {
//       stan::math::set_zero_all_adjoints();
//       y[i][j].grad(ts_var, g);
//       for (size_t k = 0; k < ts.size(); ++k) {
//         if (k == i) {
//           EXPECT_FLOAT_EQ(g[k], fval[j]);
//         } else {
//           EXPECT_FLOAT_EQ(g[k], 0.0);
//         }
//       }
//     }
//   }
// }

// TEST_F(TorstenOdeTest_chem, fwd_sensitivity_theta_y0_ts) {
//   using stan::math::value_of;
//   using Ode = PMXVariadicOdeSystem<F_eigen, var, var, std::vector<var>, std::vector<double>, std::vector<int> >;

//   nested_rev_autodiff nested;

//   PMXCvodesIntegrator<CV_BDF, CV_STAGGERED> solver(rtol, atol, 1e8);
//   std::vector<var> theta_var = stan::math::to_var(theta);
//   Eigen::Matrix<var, -1, 1> y0_vec_var = stan::math::to_var(y0_vec);
//   std::vector<var> ts_var = stan::math::to_var(ts);
//   Ode ode(f_eigen, t0, ts_var, y0_vec_var, msgs, theta_var, x_r, x_i);
//   OdeObserver<Ode> observer(ode);
//   solver.integrate(ode, observer);
//   auto& y = observer.y;
//   auto y1 = stan::math::ode_bdf(f_eigen, y0_vec_var, t0, ts, nullptr, theta_var, x_r, x_i);

//   EXPECT_ARRAY2D_ADJ_NEAR(y , y1, y0_vec_var   , nested, 1.E-5, "y0");
//   EXPECT_ARRAY2D_ADJ_NEAR(y , y1, theta_var, nested, 1.E-4, "THETA");

//   std::vector<double> y_i(y[0].size());

//   std::vector<double> g(ts.size()), fval(y0.size());
//   for (size_t i = 0; i < ts.size(); ++i) {
//     for (size_t m = 0; m < y_i.size(); ++m) {
//       y_i[m] = y[i](m).val();
//     }
//     fval = f(value_of(ts[i]), y_i, theta, x_r, x_i, msgs);
//     for (size_t j = 0; j < y0.size(); ++j) {
//       stan::math::set_zero_all_adjoints();
//       y[i][j].grad(ts_var, g);
//       for (size_t k = 0; k < ts.size(); ++k) {
//         if (k == i) {
//           EXPECT_FLOAT_EQ(g[k], fval[j]);
//         } else {
//           EXPECT_FLOAT_EQ(g[k], 0.0);
//         }
//       }
//     }
//   }
// }

// TEST_F(TorstenOdeTest_chem, fwd_sensitivity_y0_ts) {
//   using stan::math::value_of;
//   using Ode = PMXVariadicOdeSystem<F_eigen, var, var, std::vector<double>, std::vector<double>, std::vector<int> >;

//   nested_rev_autodiff nested;

//   PMXCvodesIntegrator<CV_BDF, CV_STAGGERED> solver(rtol, atol, 1e8);
//   Eigen::Matrix<var, -1, 1> y0_vec_var = stan::math::to_var(y0_vec);
//   std::vector<var> ts_var = stan::math::to_var(ts);
//   Ode ode(f_eigen, t0, ts_var, y0_vec_var, msgs, theta, x_r, x_i);
//   OdeObserver<Ode> observer(ode);
//   solver.integrate(ode, observer);
//   auto& y = observer.y;
//   auto y1 = stan::math::ode_bdf(f_eigen, y0_vec_var, t0, ts_var, nullptr, theta, x_r, x_i);

//   EXPECT_ARRAY2D_ADJ_NEAR(y , y1, y0_vec_var   , nested, 1.E-5, "y0");

//   std::vector<double> y_i(y[0].size());

//   std::vector<double> g(ts.size()), fval(y0.size());
//   for (size_t i = 0; i < ts.size(); ++i) {
//     for (size_t m = 0; m < y_i.size(); ++m) {
//       y_i[m] = y[i](m).val();
//     }
//     fval = f(value_of(ts[i]), y_i, theta, x_r, x_i, msgs);
//     for (size_t j = 0; j < y0.size(); ++j) {
//       stan::math::set_zero_all_adjoints();
//       y[i][j].grad(ts_var, g);
//       for (size_t k = 0; k < ts.size(); ++k) {
//         if (k == i) {
//           EXPECT_FLOAT_EQ(g[k], fval[j]);
//         } else {
//           EXPECT_FLOAT_EQ(g[k], 0.0);
//         }
//       }
//     }
//   }
// }

// TEST_F(TorstenOdeTest_chem, fwd_sensitivity_theta_ts) {
//   using stan::math::value_of;
//   using Ode = PMXVariadicOdeSystem<F_eigen, var, double, std::vector<var>, std::vector<double>, std::vector<int> >;

//   nested_rev_autodiff nested;

//   PMXCvodesIntegrator<CV_BDF, CV_STAGGERED> solver(rtol, atol, 1e8);
//   std::vector<var> theta_var = stan::math::to_var(theta);
//   std::vector<var> ts_var = stan::math::to_var(ts);
//   Ode ode(f_eigen, t0, ts_var, y0_vec, msgs, theta_var, x_r, x_i);
//   OdeObserver<Ode> observer(ode);
//   solver.integrate(ode, observer);
//   auto& y = observer.y;
//   auto y1 = stan::math::ode_bdf(f_eigen, y0_vec, t0, ts_var, nullptr, theta_var, x_r, x_i);

//   EXPECT_ARRAY2D_ADJ_NEAR(y , y1, theta_var, nested, 1.E-5, "THETA");

//   std::vector<double> y_i(y[0].size());

//   std::vector<double> g(ts.size()), fval(y0.size());
//   for (size_t i = 0; i < ts.size(); ++i) {
//     for (size_t m = 0; m < y_i.size(); ++m) {
//       y_i[m] = y[i](m).val();
//     }
//     fval = f(value_of(ts[i]), y_i, theta, x_r, x_i, msgs);
//     for (size_t j = 0; j < y0.size(); ++j) {
//       stan::math::set_zero_all_adjoints();
//       y[i][j].grad(ts_var, g);
//       for (size_t k = 0; k < ts.size(); ++k) {
//         if (k == i) {
//           EXPECT_FLOAT_EQ(g[k], fval[j]);
//         } else {
//           EXPECT_FLOAT_EQ(g[k], 0.0);
//         }
//       }
//     }
//   }
// }

TEST_F(TorstenOdeTest_sho, fwd_sensitivity_theta_ts) {
  using stan::math::value_of;
  using Ode = PMXVariadicOdeSystem<F_eigen, var, double, std::vector<var>, std::vector<double>, std::vector<int> >;

  nested_rev_autodiff nested;

  PMXCvodesIntegrator<CV_ADAMS, CV_STAGGERED> solver(rtol, atol, 1e8);
  std::vector<var> theta_var = stan::math::to_var(theta);
  std::vector<var> ts_var = stan::math::to_var(ts);
  Ode ode(f_eigen, t0, ts_var, y0_vec, msgs, theta_var, x_r, x_i);
  OdeObserver<Ode> observer(ode);
  solver.integrate(ode, observer);
  auto& y = observer.y;
  auto y1 = stan::math::ode_adams(f_eigen, y0_vec, t0, ts_var, nullptr,theta_var, x_r, x_i);

  EXPECT_ARRAY2D_ADJ_NEAR(y , y1, theta_var, nested, 1.E-5, "THETA");

  std::vector<double> y_i(y[0].size());

  std::vector<double> g(ts.size()), fval(y0.size());
  for (size_t i = 0; i < ts.size(); ++i) {
    for (size_t m = 0; m < y_i.size(); ++m) {
      y_i[m] = y[i](m).val();
    }
    fval = f(value_of(ts[i]), y_i, theta, x_r, x_i, msgs);
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


TEST_F(TorstenOdeTest_sho, integrate_ode_adams_theta_ts) {
  using torsten::pmx_ode_adams;
  using torsten::pmx_ode_bdf;
  using stan::math::value_of;

  nested_rev_autodiff nested;

  std::vector<var> theta_var = stan::math::to_var(theta);
  std::vector<var> ts_var = stan::math::to_var(ts);

  std::vector<Eigen::Matrix<stan::math::var, -1, 1> > y, y1;

  for (int i = 0; i < 3; ++i) {
    y = pmx_ode_adams(f_eigen, y0_vec, t0, ts_var, msgs, theta_var, x_r, x_i);
  }
  y1 = stan::math::ode_adams(f_eigen, y0_vec, t0, ts, nullptr, theta_var, x_r, x_i);

  EXPECT_ARRAY2D_ADJ_NEAR(y , y1, theta_var, nested, 1.E-5, "THETA");

  std::vector<double> y_i(y1[0].size());

  std::vector<double> g(ts.size()), fval(y0.size());
  for (size_t i = 0; i < ts.size(); ++i) {
    for (size_t m = 0; m < y_i.size(); ++m) {
      y_i[m] = y[i][m].val();
    }
    fval = f(value_of(ts[i]), y_i, theta, x_r, x_i, msgs);
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
