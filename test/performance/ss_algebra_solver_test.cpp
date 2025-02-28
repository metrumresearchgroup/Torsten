#include <stan/math/prim/fun/Eigen.hpp>
#include <stan/math/torsten/test/unit/pmx_ode_test_fixture.hpp>
#include <stan/math/torsten/test/unit/pmx_onecpt_test_fixture.hpp>
#include <stan/math/torsten/test/unit/pmx_twocpt_test_fixture.hpp>
#include <stan/math/torsten/test/unit/pmx_cpt_model_test_fixture.hpp>
#include <stan/math/torsten/test/unit/pmx_friberg_karlsson_test_fixture.hpp>
#include <stan/math/torsten/test/unit/expect_near_matrix_eq.hpp>
#include <stan/math/torsten/test/unit/expect_matrix_eq.hpp>
#include <stan/math/torsten/pmx_solve_rk45.hpp>
#include <stan/math/torsten/pmx_solve_bdf.hpp>
#include <stan/math/torsten/pmx_solve_adams.hpp>
#include <stan/math/torsten/pmx_twocpt_model.hpp>
#include <stan/math/torsten/pmx_onecpt_model.hpp>
#include <stan/math/torsten/pmx_ode_model.hpp>
#include <stan/math/torsten/test/unit/test_util.hpp>
#include <gtest/gtest.h>

auto f_onecpt = torsten::PMXOneCptModel<double>::f_;
auto f_twocpt = torsten::PMXTwoCptModel<double>::f_;

using stan::math::var;
using std::vector;
using Eigen::Matrix;
using Eigen::MatrixXd;
using Eigen::Dynamic;

using torsten::pmx_solve_rk45;
using torsten::pmx_solve_bdf;
using torsten::pmx_solve_adams;
using torsten::NONMENEventsRecord;
using torsten::PKODEModel;
using torsten::dsolve::PMXOdeIntegrator;
using torsten::dsolve::PMXOdeintIntegrator;


using torsten::dsolve::PMXCvodesIntegrator;
using torsten::PKRec;

TEST_F(TorstenTwoCptModelTest, algebra_solver_perf_bolus) {
  using torsten::PMXTwoCptODE;

  const int ncmt = 3;
  PMXTwoCptODE f;
  const std::vector<double> par{CL, Q, V2, V3, ka};

  PKODEModel<double, PMXTwoCptODE> model(par, ncmt, f);
  double ss_amt = 80 * 1000.0;
  double ss_rate = 0.0;
  double ss_ii = 8.5;
  const int n_dose = 50;

  const torsten::dsolve::PMXOdeIntegrator<PMXOdeSystem,  PMXCvodesIntegrator<CV_ADAMS, CV_STAGGERED>> integ(1.e-8, 1.e-8, 100000, 5.e-8, 5.e-8, 100000, 0);

  for (int ss_cmt = 1; ss_cmt < ncmt + 1; ++ss_cmt) {
    // steady state solution
    auto y = model.solve(0, ss_amt, ss_rate, ss_ii, ss_cmt, integ);

    // long-integration solution
    std::cout << "long term integration" << "\n";
    PKRec<double> y1 = PKRec<double>::Zero(ncmt);
    const std::vector<double> rate(ncmt, 0.0);
    double t0 = 0;
    double t1 = t0 + ss_ii;
    for (int i = 0; i < n_dose; ++i) {
      y1(ss_cmt - 1) += ss_amt;
      model.solve(y1, t0, t1, rate, integ); 
      t0 = t1;
      t1 += ss_ii;
    }

    for (int i = 1; i < ncmt; ++i) {
      EXPECT_NEAR(y[i], y1[i], 1.e-6);
    }
  }
}

TEST_F(TorstenTwoCptModelTest, algebra_solver_perf_bolus_par_sens) {
  using torsten::PMXTwoCptODE;

  const int ncmt = 3;
  PMXTwoCptODE f;
  std::vector<stan::math::var> par{CL, Q, V2, V3, ka};

  PKODEModel<stan::math::var, PMXTwoCptODE> model(par, ncmt, f);
  double ss_amt = 80 * 1000.0;
  double ss_rate = 0.0;
  double ss_ii = 8.5;
  const int n_dose = 50;

  const torsten::dsolve::PMXOdeIntegrator<PMXOdeSystem,  PMXCvodesIntegrator<CV_ADAMS, CV_STAGGERED>> integ(1.e-8, 1.e-8, 100000, 5.e-8, 5.e-8, 100000, 0);

  for (int ss_cmt = 1; ss_cmt < ncmt + 1; ++ss_cmt) {
    // steady state solution
    auto y = model.solve(0, ss_amt, ss_rate, ss_ii, ss_cmt, integ);

    // long-integration solution
    std::cout << "long term integration" << "\n";
    PKRec<stan::math::var> y1 = PKRec<stan::math::var>::Zero(ncmt);
    const std::vector<double> rate(ncmt, 0.0);
    double t0 = 0;
    double t1 = t0 + ss_ii;
    for (int i = 0; i < n_dose; ++i) {
      y1(ss_cmt - 1) += ss_amt;
      model.solve(y1, t0, t1, rate, integ); 
      t0 = t1;
      t1 += ss_ii;
    }

    torsten::test::test_grad(par, y, y1, 1.e-6, 1.e-6);
  }
}

TEST_F(TorstenTwoCptModelTest, algebra_solver_perf_truncated_infusion) {
  using torsten::PMXTwoCptODE;

  const int ncmt = 3;
  PMXTwoCptODE f;
  const std::vector<stan::math::var> par{CL, Q, V2, V3, ka};

  PKODEModel<stan::math::var, PMXTwoCptODE> model(par, ncmt, f);
  double ss_amt = 80 * 1000.0;
  double ss_rate = 15000.0;
  double ss_ii = 4.5;
  const int n_dose = 50;

  using scheme_t = boost::numeric::odeint::runge_kutta_dopri5<std::vector<double>, double, std::vector<double>, double>;
  const torsten::dsolve::PMXOdeIntegrator<PMXOdeSystem, PMXOdeintIntegrator<scheme_t>> integ(1.e-8, 1.e-8, 100000, 5.e-8, 5.e-8, 100000, 0);

  const int ss_cmt = 2;
  // steady state solution
  auto y = model.solve(0, ss_amt, ss_rate, ss_ii, ss_cmt, integ);
  EXPECT_FLOAT_EQ(y(ss_cmt - 1).val(), 138380.44);
}

TEST_F(FribergKarlssonTest, algebra_solver_perf_bolus) {
  const int npar = 9;
  const PMXOdeIntegrator<PMXOdeSystem,  PMXCvodesIntegrator<CV_ADAMS, CV_STAGGERED>> integ(1.e-8, 1.e-8, 100000, 5.e-8, 5.e-8, 10000, 0); 
  std::vector<double> par(npar);
  std::vector<double> y0(nCmt);
  par[0] = 7.25654;             // CL
  par[1] = 17.0481;           // Q
  par[2] = 107.166;            // VC
  par[3] = 124.171;           // VP
  par[4] = 1.47284;           // ka
  par[5] = 117.039;           // MTT
  par[6] = 5.49146;           // circ0
  par[7] = 0.172795;          // gamma
  par[8] = 0.00014864;        // alpha
  PKODEModel<double, FribergKarlsson> model(par, nCmt, f);
  double ss_amt = 8000.0;
  double ss_rate = 0.0;
  double ss_ii = 504;
  int ss_cmt = 2;

  // steady state solution
  auto y = model.solve(0, ss_amt, ss_rate, ss_ii, ss_cmt, integ);

  // long-integration solution
  int n_dose = 30;
  std::cout << "torsten info: " << "long integ" << "\n";
  PKRec<double> y1 = PKRec<double>::Zero(nCmt);
  std::vector<double> rate(nCmt, 0.0);
    double t0 = 0;
    double t1 = t0 + ss_ii;
    for (int i = 0; i < n_dose; ++i) {
      y1(ss_cmt - 1) += ss_amt;
      model.solve(y1, t0, t1, rate, integ); 
      t0 = t1;
      t1 += ss_ii;
    }

  for (int i = 1; i < nCmt; ++i) {
    EXPECT_NEAR(y[i], y1[i], 1.e-6);
  }
}

TEST_F(FribergKarlssonTest, algebra_solver_perf_truncated_infusion) {
  const int npar = 9;
  const PMXOdeIntegrator<PMXOdeSystem,  PMXCvodesIntegrator<CV_ADAMS, CV_STAGGERED>> integ(1.e-8, 1.e-8, 100000, 5.e-8, 5.e-8, 10000, 0); 
  std::vector<double> par(npar);
  std::vector<double> y0(nCmt);
  par[0] = 7.25654;             // CL
  par[1] = 17.0481;           // Q
  par[2] = 107.166;            // VC
  par[3] = 124.171;           // VP
  par[4] = 1.47284;           // ka
  par[5] = 117.039;           // MTT
  par[6] = 5.49146;           // circ0
  par[7] = 0.172795;          // gamma
  par[8] = 0.00014864;        // alpha
  PKODEModel<double, FribergKarlsson> model(par, nCmt, f);
  double ss_amt = 8000.0;
  double ss_rate = 4000.0;
  double ss_ii = 504;
  int ss_cmt = 2;

  // steady state solution
  auto y = model.solve(0, ss_amt, ss_rate, ss_ii, ss_cmt, integ);

  // long-integration solution
  std::cout << "torsten info: " << "long integ" << "\n";
  PKRec<double> y1 = PKRec<double>::Zero(nCmt);
  std::vector<double> rate(nCmt, 0.0);
  double t0 = 0, t1 = 0;
  for (int i = 0; i < 30; ++i) {
    t0 = t1;
    t1 = t0 + ss_amt/ss_rate;
    rate[ss_cmt - 1] = ss_rate;
    model.solve(y1, t0, t1, rate, integ); 
    t0 = t1;
    t1 += ss_ii - ss_amt/ss_rate;
    rate[ss_cmt - 1] = 0.0;
    model.solve(y1, t0, t1, rate, integ); 
  }

  for (int i = 1; i < nCmt; ++i) {
    EXPECT_NEAR(y[i], y1[i], 1.e-6);
  }
}

TEST_F(FribergKarlssonTest, algebra_solver_perf_stationary_solution) {
  const int npar = 9;
  using scheme_t = boost::numeric::odeint::runge_kutta_dopri5<std::vector<double>, double, std::vector<double>, double>;
  const PMXOdeIntegrator<PMXOdeSystem, PMXOdeintIntegrator<scheme_t>> integ(1.e-8, 1.e-8, 100000, 5.e-8, 5.e-8, 10000, 0); 
  std::vector<double> par(npar);
  std::vector<double> y0(nCmt);
  par[0] = 7.25654;             // CL
  par[1] = 17.0481;           // Q
  par[2] = 107.166;            // VC
  par[3] = 124.171;           // VP
  par[4] = 1.47284;           // ka
  par[5] = 117.039;           // MTT
  par[6] = 5.49146;           // circ0
  par[7] = 0.001;          // gamma
  par[8] = 0.1;        // alpha
  PKODEModel<double, FribergKarlsson> model(par, nCmt, f);
  double ss_amt = 8000.0;
  double ss_rate = 0.0;
  double ss_ii = 504;
  int ss_cmt = 2;

  // steady state solution
  auto y = model.solve(0, ss_amt, ss_rate, ss_ii, ss_cmt, integ);

  // long-integration solution
  int n_dose = 30;
  std::cout << "torsten info: " << "long integ" << "\n";
  PKRec<double> y1 = PKRec<double>::Zero(nCmt);
  std::vector<double> rate(nCmt, 0.0);
    double t0 = 0;
    double t1 = t0 + ss_ii;
    for (int i = 0; i < n_dose; ++i) {
      y1(ss_cmt - 1) += ss_amt;
      model.solve(y1, t0, t1, rate, integ); 
      t0 = t1;
      t1 += ss_ii;
    }

  for (int i = 1; i < nCmt; ++i) {
    EXPECT_NEAR(y[i], y1[i], 1.e-6);
  }
}
