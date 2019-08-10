#include <stan/math/prim/mat/fun/Eigen.hpp>
#include <test/unit/math/torsten/pmx_ode_test_fixture.hpp>
#include <test/unit/math/torsten/pmx_onecpt_test_fixture.hpp>
#include <test/unit/math/torsten/pmx_twocpt_test_fixture.hpp>
#include <test/unit/math/torsten/pmx_friberg_karlsson_test_fixture.hpp>
#include <test/unit/math/torsten/expect_near_matrix_eq.hpp>
#include <test/unit/math/torsten/expect_matrix_eq.hpp>
#include <stan/math/torsten/pmx_solve_rk45.hpp>
#include <stan/math/torsten/pmx_solve_bdf.hpp>
#include <stan/math/torsten/pmx_twocpt_model.hpp>
#include <stan/math/torsten/pmx_onecpt_model.hpp>
#include <stan/math/torsten/pmx_ode_model.hpp>
#include <test/unit/math/torsten/util_generalOdeModel.hpp>
#include <gtest/gtest.h>

auto f  = refactor::PMXOneCptModel<double,double,double,double>::f_;
auto f2 = refactor::PMXTwoCptModel<double,double,double,double>::f_;

using stan::math::var;
using std::vector;
using Eigen::Matrix;
using Eigen::MatrixXd;
using Eigen::Dynamic;

using torsten::pmx_solve_rk45;
using torsten::pmx_solve_bdf;
using torsten::pmx_solve_adams;
using torsten::NONMENEventsRecord;

TEST_F(TorstenOneCptTest, ss_zero_rate) {
  // Steady state induced by multiple bolus doses (SS = 1, rate = 0)
  using std::vector;
  using Eigen::Matrix;
  using Eigen::MatrixXd;
  using Eigen::Dynamic;

  time[0] = 0.0;
  time[1] = 0.0;
  for(int i = 2; i < 10; i++) time[i] = time[i - 1] + 5;

  amt[0] = 1200;
  addl[0] = 10;
  ss[0] = 1;

  double rel_tol = 1e-8, abs_tol = 1e-8;
  long int max_num_steps = 1e8;
  MatrixXd x_rk45 = torsten::pmx_solve_rk45(f, nCmt,
                                time, amt, rate, ii, evid, cmt, addl, ss,
                                pMatrix, biovar, tlag,
                                0,
                                rel_tol, abs_tol, max_num_steps);

  MatrixXd x_bdf = torsten::pmx_solve_bdf(f, nCmt,
                              time, amt, rate, ii, evid, cmt, addl, ss,
                              pMatrix, biovar, tlag,
                              0,
                              rel_tol, abs_tol, max_num_steps);

  Eigen::MatrixXd x(10, 2);
  x << 1200.0      , 384.7363,
    1200.0      , 384.7363,
    2.974504    , 919.6159,
    7.373062e-3 , 494.0040,
    3.278849e+1 , 1148.4725,
    8.127454e-2 , 634.2335,
    3.614333e+2 , 1118.2043,
    8.959035e-1 , 813.4883,
    2.220724e-3 , 435.9617,
    9.875702    , 1034.7998;
  MatrixXd xt = x.transpose();

  torsten::test::test_val(x_rk45, xt, 1e-6, 1e-4);
  torsten::test::test_val(x_bdf, xt, 1e-4, 1e-4);

  TORSTEN_ODE_GRAD_THETA_TEST(pmx_solve_rk45, f, nCmt,
                              time, amt, rate, ii, evid, cmt, addl, ss, pMatrix, biovar, tlag,
                              rel_tol, abs_tol, max_num_steps,
                              2e-5, 1e-6, 1e-5, 1e-5);

  TORSTEN_ODE_GRAD_THETA_TEST(pmx_solve_bdf, f, nCmt,
                              time, amt, rate, ii, evid, cmt, addl, ss, pMatrix, biovar, tlag,
                              rel_tol, abs_tol, max_num_steps,
                              2e-5, 1e-6, 5e-5, 1e-5);

  TORSTEN_ODE_GRAD_THETA_TEST(pmx_solve_adams, f, nCmt,
                              time, amt, rate, ii, evid, cmt, addl, ss, pMatrix, biovar, tlag,
                              rel_tol, abs_tol, max_num_steps,
                              2e-5, 1e-4, 1e-3, 1e-4);

  TORSTEN_ODE_GRAD_BIOVAR_TEST(pmx_solve_rk45, f, nCmt,
                              time, amt, rate, ii, evid, cmt, addl, ss, pMatrix, biovar, tlag,
                              rel_tol, abs_tol, max_num_steps,
                              2e-5, 1e-6, 1e-6, 1e-6);

  TORSTEN_ODE_GRAD_BIOVAR_TEST(pmx_solve_bdf, f, nCmt,
                              time, amt, rate, ii, evid, cmt, addl, ss, pMatrix, biovar, tlag,
                              rel_tol, abs_tol, max_num_steps,
                              2e-5, 1e-6, 1e-4, 1e-6);

  TORSTEN_ODE_GRAD_BIOVAR_TEST(pmx_solve_adams, f, nCmt,
                              time, amt, rate, ii, evid, cmt, addl, ss, pMatrix, biovar, tlag,
                              rel_tol, abs_tol, max_num_steps,
                              2e-5, 1e-6, 1e-5, 1e-6);
}

TEST_F(TorstenOneCptTest, single_bolus_tlag) {
  nt = 2;
  time.resize(nt);
  amt.resize(nt);
  rate.resize(nt);
  cmt.resize(nt);
  evid.resize(nt);
  ii.resize(nt);
  addl.resize(nt);
  ss.resize(nt);
  evid[0] = 1;
  cmt[0] = 1;
  ii[0] = 0;
  addl[0] = 0;
  time[0] = 0.0;
  tlag[0][0] = 1.5;

  time[1] = 2.5;

  double rel_tol = 1e-8, abs_tol = 1e-8;
  long int max_num_steps = 1e8;

  TORSTEN_ODE_GRAD_TLAG_TEST(pmx_solve_bdf, f, nCmt,
                             time, amt, rate, ii, evid, cmt, addl, ss, pMatrix, biovar, tlag,
                             rel_tol, abs_tol, max_num_steps,
                             2e-5, 1e-10, 1e-4, 1e-5);

  TORSTEN_ODE_GRAD_TLAG_TEST(pmx_solve_adams, f, nCmt,
                             time, amt, rate, ii, evid, cmt, addl, ss, pMatrix, biovar, tlag,
                             rel_tol, abs_tol, max_num_steps,
                             2e-5, 1e-10, 1e-4, 1e-5);
}

TEST_F(TorstenOneCptTest, multiple_bolus_tlag) {
  time[0] = 0.0;
  time[1] = 0.0;
  for(int i = 2; i < 10; i++) time[i] = time[i - 1] + 5;
  time[3] = 7.0;

  amt[0] = 0.0;
  evid[0] = 0;
  evid[1] = 1;
  evid[3] = 1;
  amt[1] = 1200;
  amt[3] = 1000;
  cmt[1] = 1;
  cmt[3] = 1;
  ss[1] = 0;
  tlag.resize(nt);
  for (int i = 0; i < nt; ++i) {
    tlag[i].resize(nCmt);
    tlag[i][0] = 0.0;
    tlag[i][1] = 0.0;
  }
  tlag[1][0] = 1.0;
  tlag[3][0] = 2.5;
  ii[0] = 0.0;
  addl[0] = 0;

  double rel_tol = 1e-10, abs_tol = 1e-10;
  long int max_num_steps = 1e8;

  TORSTEN_ODE_GRAD_TLAG_TEST(pmx_solve_bdf, f, nCmt,
                             time, amt, rate, ii, evid, cmt, addl, ss, pMatrix, biovar, tlag,
                             rel_tol, abs_tol, max_num_steps,
                             2e-5, 1e-10, 1e-4, 1e-5);

  TORSTEN_ODE_GRAD_TLAG_TEST(pmx_solve_adams, f, nCmt,
                             time, amt, rate, ii, evid, cmt, addl, ss, pMatrix, biovar, tlag,
                             rel_tol, abs_tol, max_num_steps,
                             2e-5, 1e-7, 1e-4, 1e-5);
}

TEST_F(TorstenOneCptTest, single_iv_tlag) {
  nt = 2;
  time.resize(nt);
  amt.resize(nt);
  rate.resize(nt);
  cmt.resize(nt);
  evid.resize(nt);
  ii.resize(nt);
  addl.resize(nt);
  ss.resize(nt);
  evid[0] = 1;
  cmt[0] = 2;
  ii[0] = 0;
  addl[0] = 0;
  time[0] = 0.0;
  tlag[0][0] = 1.5;
  tlag[0][1] = 1.5;
  rate[0] = 300;

  time[1] = 2.5;

  double rel_tol = 1e-10, abs_tol = 1e-10;
  long int max_num_steps = 1e8;

  TORSTEN_ODE_GRAD_TLAG_TEST(pmx_solve_bdf, f, nCmt,
                             time, amt, rate, ii, evid, cmt, addl, ss, pMatrix, biovar, tlag,
                             rel_tol, abs_tol, max_num_steps,
                             2e-5, 1e-9, 1e-6, 1e-6);

  TORSTEN_ODE_GRAD_TLAG_TEST(pmx_solve_adams, f, nCmt,
                             time, amt, rate, ii, evid, cmt, addl, ss, pMatrix, biovar, tlag,
                             rel_tol, abs_tol, max_num_steps,
                             2e-5, 1e-8, 1e-6, 1e-6);
}

TEST_F(TorstenOneCptTest, multiple_iv_tlag) {
  nt = 4;
  time.resize(nt);
  amt.resize(nt);
  rate.resize(nt);
  cmt.resize(nt);
  evid.resize(nt);
  ii.resize(nt);
  addl.resize(nt);
  ss.resize(nt);
  tlag[0] = std::vector<double>{1.3, 0.8};

  time[0] = 0.0;
  time[1] = 0.0;
  for(int i = 2; i < nt; i++) time[i] = time[i - 1] + 5;

  amt[0] = 0.0;
  evid[0] = 0;
  evid[1] = 1;
  evid[2] = 1;
  amt[1] = 1200;
  amt[2] = 1000;
  rate[1] = 230;
  rate[2] = 280;
  cmt[1] = 2;
  cmt[2] = 2;
  ss[1] = 0;
  ii[0] = 0.0;
  addl[0] = 0;

  double rel_tol = 1e-10, abs_tol = 1e-10;
  long int max_num_steps = 1e8;

  TORSTEN_ODE_GRAD_TLAG_TEST(pmx_solve_bdf, f, nCmt,
                             time, amt, rate, ii, evid, cmt, addl, ss, pMatrix, biovar, tlag,
                             rel_tol, abs_tol, max_num_steps,
                             2e-5, 1e-5, 1e-6, 1e-6);

  TORSTEN_ODE_GRAD_TLAG_TEST(pmx_solve_adams, f, nCmt,
                             time, amt, rate, ii, evid, cmt, addl, ss, pMatrix, biovar, tlag,
                             rel_tol, abs_tol, max_num_steps,
                             2e-5, 1e-7, 1e-6, 1e-6);
}

TEST_F(TorstenOneCptTest, ss_multiple_infusion_tlag) {
  time[0] = 0.0;
  time[1] = 0.0;
  for(int i = 2; i < nt; i++) time[i] = time[i - 1] + 5;
  time[3] = 6.0;

  amt[0] = 0.0;
  evid[0] = 0;
  evid[1] = 1;
  evid[3] = 1;
  amt[1] = 1200;
  amt[3] = 1000;
  rate[1] = 200;
  rate[3] = 200;
  cmt[1] = 1;
  cmt[3] = 1;
  ss[1] = 1;
  tlag.resize(nt);
  for (int i = 0; i < nt; ++i) {
    tlag[i].resize(nCmt);
    tlag[i][0] = 0.0;
    tlag[i][1] = 0.0;
  }
  tlag[1][0] = 1.0;
  tlag[3][0] = 0.5;
  ii[1] = 6.5;
  addl[0] = 0;

  double rel_tol = 1e-8, abs_tol = 1e-8;
  long int max_num_steps = 1e8;

  {
    auto f1 = [&] (std::vector<double>& x) {
      std::vector<std::vector<double> > tlag1(nt);
      for (int i = 0; i < nt; ++i) tlag1[i] = tlag[i];
      tlag1[3] = x;
      NONMENEventsRecord<double, double, double, double, std::vector<double>, double, double>
      events_rec(nCmt, time, amt, rate, ii, evid, cmt, addl, ss, pMatrix, biovar, tlag1);
      torsten::EventsManager<NONMENEventsRecord<double, double, double, double, std::vector<double>, double, double>>
      em(events_rec);
      return torsten::pmx_solve_bdf(f, nCmt, time, amt, rate, ii, evid, cmt, addl, ss, pMatrix, biovar, tlag1,
                                          0, rel_tol, abs_tol, max_num_steps);
    };
    auto f2 = [&] (std::vector<stan::math::var>& x) {
      std::vector<std::vector<var> > tlag1(nt);
      for (int i = 0; i < nt; ++i) tlag1[i] = stan::math::to_var(tlag[i]);
      tlag1[3] = x;
      NONMENEventsRecord<double, double, double, double, std::vector<double>, double, var>
      events_rec(nCmt, time, amt, rate, ii, evid, cmt, addl, ss, pMatrix, biovar, tlag1);
      torsten::EventsManager<NONMENEventsRecord<double, double, double, double, std::vector<double>, double, var>>
      em(events_rec);
      return torsten::pmx_solve_bdf(f, nCmt, time, amt, rate, ii, evid, cmt, addl, ss, pMatrix, biovar, tlag1,
                                          0, rel_tol, abs_tol, max_num_steps);
    };
    std::vector<double> tlag_test(tlag[3]);
    torsten::test::test_grad(f1, f2, tlag_test, 2e-5, 1e-5, 1e-3, 1e-5);
  }
}

TEST_F(TorstenOneCptTest, ss_bolus) {
  resize(3);
  time[0] = 0.0;
  time[1] = 0.0;
  for(int i = 2; i < nt; i++) time[i] = time[i - 1] + 5;

  amt[0] = 1200;
  addl[0] = 2;
  ss[0] = 1;

  double rel_tol = 1e-12, abs_tol = 1e-12;
  long int max_num_steps = 1e8;

  biovar[0] = std::vector<double>{0.8, 0.9};
  tlag[0] = std::vector<double>{2.8, 3.9};

  TORSTEN_ODE_GRAD_THETA_TEST(pmx_solve_rk45, f, nCmt,
                              time, amt, rate, ii, evid, cmt, addl, ss, pMatrix, biovar, tlag,
                              rel_tol, abs_tol, max_num_steps,
                              2e-5, 1e-6, 1e-5, 1e-5);

  TORSTEN_ODE_GRAD_THETA_TEST(pmx_solve_bdf, f, nCmt,
                              time, amt, rate, ii, evid, cmt, addl, ss, pMatrix, biovar, tlag,
                              rel_tol, abs_tol, max_num_steps,
                              2e-5, 1e-6, 5e-5, 1e-5);

  TORSTEN_ODE_GRAD_THETA_TEST(pmx_solve_adams, f, nCmt,
                              time, amt, rate, ii, evid, cmt, addl, ss, pMatrix, biovar, tlag,
                              rel_tol, abs_tol, max_num_steps,
                              2e-5, 1e-4, 1e-3, 1e-4);

  TORSTEN_ODE_GRAD_BIOVAR_TEST(pmx_solve_rk45, f, nCmt,
                              time, amt, rate, ii, evid, cmt, addl, ss, pMatrix, biovar, tlag,
                              rel_tol, abs_tol, max_num_steps,
                              2e-5, 1e-6, 1e-6, 1e-6);

  TORSTEN_ODE_GRAD_BIOVAR_TEST(pmx_solve_bdf, f, nCmt,
                              time, amt, rate, ii, evid, cmt, addl, ss, pMatrix, biovar, tlag,
                              rel_tol, abs_tol, max_num_steps,
                              2e-5, 1e-6, 1e-6, 1e-6);

  TORSTEN_ODE_GRAD_BIOVAR_TEST(pmx_solve_adams, f, nCmt,
                              time, amt, rate, ii, evid, cmt, addl, ss, pMatrix, biovar, tlag,
                              rel_tol, abs_tol, max_num_steps,
                              2e-5, 1e-6, 1e-5, 1e-6);
}

TEST_F(TorstenOneCptTest, ss_constant_infusion) {
  time[0] = 0.0;
  for(int i = 1; i < 10; i++) time[i] = time[i - 1] + 2.5;

  amt[0] = 0.0;
  rate[0] = 150;
  ii[0] = 0.0;
  addl[0] = 0;
  ss[0] = 1;

  double rel_tol = 1e-12, abs_tol = 1e-12;
  long int max_num_steps = 1e8;

  TORSTEN_ODE_GRAD_THETA_TEST(pmx_solve_rk45, f, nCmt,
                              time, amt, rate, ii, evid, cmt, addl, ss, pMatrix, biovar, tlag,
                              rel_tol, abs_tol, max_num_steps,
                              2e-5, 1e-6, 1e-5, 1e-5);

  TORSTEN_ODE_GRAD_THETA_TEST(pmx_solve_bdf, f, nCmt,
                              time, amt, rate, ii, evid, cmt, addl, ss, pMatrix, biovar, tlag,
                              rel_tol, abs_tol, max_num_steps,
                              2e-5, 1e-6, 5e-5, 1e-5);

  TORSTEN_ODE_GRAD_THETA_TEST(pmx_solve_adams, f, nCmt,
                              time, amt, rate, ii, evid, cmt, addl, ss, pMatrix, biovar, tlag,
                              rel_tol, abs_tol, max_num_steps,
                              2e-5, 1e-4, 1e-3, 1e-4);

  TORSTEN_ODE_GRAD_BIOVAR_TEST(pmx_solve_rk45, f, nCmt,
                              time, amt, rate, ii, evid, cmt, addl, ss, pMatrix, biovar, tlag,
                              rel_tol, abs_tol, max_num_steps,
                              2e-5, 1e-6, 1e-6, 1e-6);

  TORSTEN_ODE_GRAD_BIOVAR_TEST(pmx_solve_bdf, f, nCmt,
                              time, amt, rate, ii, evid, cmt, addl, ss, pMatrix, biovar, tlag,
                              rel_tol, abs_tol, max_num_steps,
                              2e-5, 1e-6, 1e-6, 1e-6);

  TORSTEN_ODE_GRAD_BIOVAR_TEST(pmx_solve_adams, f, nCmt,
                              time, amt, rate, ii, evid, cmt, addl, ss, pMatrix, biovar, tlag,
                              rel_tol, abs_tol, max_num_steps,
                              2e-5, 1e-6, 1e-5, 1e-6);
}

TEST_F(TorstenOneCptTest, ss_multiple_infusion) {
  time[0] = 0.0;
  for(int i = 1; i < 10; i++) time[i] = time[i - 1] + 2.5;

  amt[0] = 1200;
  rate[0] = 150;
  ii[0] = 16;
  addl[0] = 0;
  ss[0] = 1;

  double rel_tol = 1e-8, abs_tol = 1e-8;
  long int max_num_steps = 1e8;
  MatrixXd x_rk45 = torsten::pmx_solve_rk45(f, nCmt,
                                time, amt, rate, ii, evid, cmt, addl, ss,
                                pMatrix, biovar, tlag,
                                0,
                                rel_tol, abs_tol, max_num_steps);

  MatrixXd x_bdf = torsten::pmx_solve_bdf(f, nCmt,
                              time, amt, rate, ii, evid, cmt, addl, ss,
                              pMatrix, biovar, tlag,
                              0,
                              rel_tol, abs_tol, max_num_steps);

  MatrixXd x(10, 2);
  x << 8.465519e-03, 360.2470,
             1.187770e+02, 490.4911,
             1.246902e+02, 676.1759,
             1.249846e+02, 816.5263,
             1.133898e+01, 750.0054,
             5.645344e-01, 557.3459,
             2.810651e-02, 408.1926,
             1.399341e-03, 298.6615,
             6.966908e-05, 218.5065,
             3.468619e-06, 159.8628;
  MatrixXd xt = x.transpose();

  torsten::test::test_val(x_rk45, xt, 1e-4, 1e-4);
  torsten::test::test_val(x_bdf, xt, 1e-3, 1e-4);

  rel_tol = 1e-12;
  abs_tol = 1e-12;

  TORSTEN_ODE_GRAD_THETA_TEST(pmx_solve_rk45, f, nCmt,
                              time, amt, rate, ii, evid, cmt, addl, ss, pMatrix, biovar, tlag,
                              rel_tol, abs_tol, max_num_steps,
                              2e-5, 1e-6, 1e-5, 1e-5);

  TORSTEN_ODE_GRAD_THETA_TEST(pmx_solve_bdf, f, nCmt,
                              time, amt, rate, ii, evid, cmt, addl, ss, pMatrix, biovar, tlag,
                              rel_tol, abs_tol, max_num_steps,
                              2e-5, 1e-6, 5e-5, 1e-5);

  TORSTEN_ODE_GRAD_THETA_TEST(pmx_solve_adams, f, nCmt,
                              time, amt, rate, ii, evid, cmt, addl, ss, pMatrix, biovar, tlag,
                              rel_tol, abs_tol, max_num_steps,
                              2e-5, 1e-4, 1e-3, 1e-4);

  biovar[0] = std::vector<double>{0.8, 0.7};

  TORSTEN_ODE_GRAD_BIOVAR_TEST(pmx_solve_rk45, f, nCmt,
                              time, amt, rate, ii, evid, cmt, addl, ss, pMatrix, biovar, tlag,
                              rel_tol, abs_tol, max_num_steps,
                              2e-5, 1e-6, 1e-6, 1e-6);

  TORSTEN_ODE_GRAD_BIOVAR_TEST(pmx_solve_bdf, f, nCmt,
                               time, amt, rate, ii, evid, cmt, addl, ss, pMatrix, biovar, tlag,
                               rel_tol, abs_tol, max_num_steps,
                               2e-5, 1e-10, 1.5e-6, 2e-8);

  TORSTEN_ODE_GRAD_BIOVAR_TEST(pmx_solve_adams, f, nCmt,
                               time, amt, rate, ii, evid, cmt, addl, ss, pMatrix, biovar, tlag,
                               rel_tol, abs_tol, max_num_steps,
                               2e-5, 1e-10, 1e-6, 1e-8);
}

TEST_F(TorstenOneCptTest, multiple_bolus) {
  ii[0] = 12;
  addl[0] = 14;

  double rel_tol = 1e-8, abs_tol = 1e-8;
  long int max_num_steps = 1e8;
  MatrixXd x_rk45 = torsten::pmx_solve_rk45(f, nCmt,
                                time, amt, rate, ii, evid, cmt, addl, ss,
                                pMatrix, biovar, tlag,
                                0,
                                rel_tol, abs_tol, max_num_steps);

  MatrixXd x_bdf = torsten::pmx_solve_bdf(f, nCmt,
                              time, amt, rate, ii, evid, cmt, addl, ss,
                              pMatrix, biovar, tlag,
                              0,
                              rel_tol, abs_tol, max_num_steps);

  MatrixXd x(10, nCmt);
  x << 1000.0, 0.0,
    740.8182, 254.97490,
    548.8116, 436.02020,
    406.5697, 562.53846,
    301.1942, 648.89603,
    223.1302, 705.72856,
    165.2989, 740.90816,
    122.4564, 760.25988,
    90.71795, 768.09246,
    8.229747, 667.87079;
  MatrixXd xt = x.transpose();

  torsten::test::test_val(x_rk45, xt, 1e-5, 1e-5);
  torsten::test::test_val(x_bdf, xt, 1e-5, 1e-5);

  rel_tol = 1e-12;
  abs_tol = 1e-12;
  ii[0] = 3.5;
  addl[0] = 2;

  TORSTEN_ODE_GRAD_THETA_TEST(pmx_solve_rk45, f, nCmt,
                              time, amt, rate, ii, evid, cmt, addl, ss, pMatrix, biovar, tlag,
                              rel_tol, abs_tol, max_num_steps,
                              2e-5, 1e-6, 1e-5, 1e-5);

  TORSTEN_ODE_GRAD_THETA_TEST(pmx_solve_bdf, f, nCmt,
                              time, amt, rate, ii, evid, cmt, addl, ss, pMatrix, biovar, tlag,
                              rel_tol, abs_tol, max_num_steps,
                              2e-5, 1e-6, 1e-3, 1e-4);

  TORSTEN_ODE_GRAD_THETA_TEST(pmx_solve_adams, f, nCmt,
                              time, amt, rate, ii, evid, cmt, addl, ss, pMatrix, biovar, tlag,
                              rel_tol, abs_tol, max_num_steps,
                              2e-5, 1e-4, 1e-3, 1e-4);

  TORSTEN_ODE_GRAD_BIOVAR_TEST(pmx_solve_rk45, f, nCmt,
                              time, amt, rate, ii, evid, cmt, addl, ss, pMatrix, biovar, tlag,
                              rel_tol, abs_tol, max_num_steps,
                              2e-5, 1e-6, 1e-6, 1e-6);

  TORSTEN_ODE_GRAD_BIOVAR_TEST(pmx_solve_bdf, f, nCmt,
                              time, amt, rate, ii, evid, cmt, addl, ss, pMatrix, biovar, tlag,
                              rel_tol, abs_tol, max_num_steps,
                              2e-5, 1e-6, 1e-6, 1e-6);

  TORSTEN_ODE_GRAD_BIOVAR_TEST(pmx_solve_adams, f, nCmt,
                              time, amt, rate, ii, evid, cmt, addl, ss, pMatrix, biovar, tlag,
                              rel_tol, abs_tol, max_num_steps,
                              2e-5, 1e-6, 1e-5, 1e-6);
}

TEST_F(TorstenOneCptTest, multiple_bolus_tlag_overload) {
  resize(3);
  addl[0] = 1;
  tlag[0] = std::vector<double>{2.8, 3.9};

  TORSTEN_ODE_PARAM_OVERLOAD_TEST(torsten::pmx_solve_bdf, f, nCmt,
                                  time, amt, rate, ii, evid, cmt, addl, ss, pMatrix, biovar, tlag,
                                  1e-10, 1e-10);
}

template <typename T0, typename T1, typename T2, typename T3>
inline
std::vector<typename boost::math::tools::promote_args<T0, T1, T2, T3>::type>
oneCptModelODE_abstime(const T0& t,
                       const std::vector<T1>& x,
	                   const std::vector<T2>& parms,
	                   const std::vector<T3>& rate,
	                   const std::vector<int>& dummy, std::ostream* pstream__) {
  typedef typename boost::math::tools::promote_args<T0, T1, T2, T3>::type scalar;

  scalar CL0 = parms[0], V1 = parms[1], ka = parms[2], CLSS = parms[3],
    K = parms[4];

  scalar CL = CL0 + (CLSS - CL0) * (1 - stan::math::exp(-K * t));
  scalar k10 = CL / V1;

  std::vector<scalar> y(2, 0);

  y[0] = -ka * x[0];
  y[1] = ka * x[0] - k10 * x[1];

  return y;
}

struct oneCptModelODE_abstime_functor {
  template <typename T0, typename T1, typename T2, typename T3>
  inline
  std::vector<typename boost::math::tools::promote_args<T0, T1, T2, T3>::type>
  operator()(const T0& t,
             const std::vector<T1>& x,
             const std::vector<T2>& parms,
             const std::vector<T3>& rate,
             const std::vector<int>& dummy, std::ostream* pstream__) const {
        return oneCptModelODE_abstime(t, x, parms, rate, dummy, pstream__);
    }
};

TEST(Torsten, genCpt_One_abstime_SingleDose) {
  using std::vector;
  using Eigen::Matrix;
  using Eigen::MatrixXd;
  using Eigen::Dynamic;

  double rel_err = 1e-6;

  vector<vector<double> > pMatrix(1);
  pMatrix[0].resize(5);
  pMatrix[0][0] = 10; // CL0
  pMatrix[0][1] = 80; // Vc
  pMatrix[0][2] = 1.2; // ka
  pMatrix[0][3] = 2; // CLSS
  pMatrix[0][4] = 1; // K

  int nCmt = 2;
  vector<vector<double> > biovar(1);
  biovar[0].resize(nCmt);
  biovar[0][0] = 1;  // F1
  biovar[0][1] = 1;  // F2

  vector<vector<double> > tlag(1);
  tlag[0].resize(nCmt);
  tlag[0][0] = 0;  // tlag1
  tlag[0][1] = 0;  // tlag2


  vector<double> time(10);
  time[0] = 0.0;
  for(int i = 1; i < 9; i++) time[i] = time[i - 1] + 0.25;
  time[9] = 4.0;

  vector<double> amt(10, 0);
  amt[0] = 1000;

  vector<double> rate(10, 0);

  vector<int> cmt(10, 2);
  cmt[0] = 1;

  vector<int> evid(10, 0);
  evid[0] = 1;

  vector<double> ii(10, 0);
  ii[0] = 12;

  vector<int> addl(10, 0);
  addl[0] = 14;

  vector<int> ss(10, 0);

  double rel_tol = 1e-8, abs_tol = 1e-8;
  long int max_num_steps = 1e8;

  oneCptModelODE_abstime_functor f;

  Matrix<double, Eigen::Dynamic, Eigen::Dynamic> x_rk45;
  x_rk45 = torsten::pmx_solve_rk45(f, 2,
                                time, amt, rate, ii, evid, cmt, addl, ss,
                                pMatrix, biovar, tlag,
                                0,
                                rel_tol, abs_tol, max_num_steps);

  Matrix<double, Eigen::Dynamic, Eigen::Dynamic> x_bdf;
  x_bdf = torsten::pmx_solve_bdf(f, 2,
                              time, amt, rate, ii, evid, cmt, addl, ss,
                              pMatrix, biovar, tlag,
                              0,
                              rel_tol, abs_tol, max_num_steps);

  MatrixXd amounts(10, 2);
  amounts << 1000.0, 0.0,
    740.8182, 255.4765,
    548.8116, 439.2755,
    406.5697, 571.5435,
    301.1942, 666.5584,
    223.1302, 734.5274,
    165.2989, 782.7979,
    122.4564, 816.6868,
    90.71795, 840.0581,
    8.229747, 869.0283;
  MatrixXd xt = amounts.transpose();

  expect_near_matrix_eq(xt, x_rk45, rel_err);
  expect_near_matrix_eq(xt, x_bdf, rel_err);

  TORSTEN_ODE_GRAD_THETA_TEST(pmx_solve_rk45, f, 2,
                              time, amt, rate, ii, evid, cmt, addl, ss, pMatrix, biovar, tlag,
                              rel_tol, abs_tol, max_num_steps,
                              2e-5, 1e-6, 1e-5, 1e-5);

  TORSTEN_ODE_GRAD_THETA_TEST(pmx_solve_bdf, f, 2,
                              time, amt, rate, ii, evid, cmt, addl, ss, pMatrix, biovar, tlag,
                              rel_tol, abs_tol, max_num_steps,
                              2e-5, 1e-6, 2e-4, 1e-4);

  TORSTEN_ODE_GRAD_THETA_TEST(pmx_solve_adams, f, 2,
                              time, amt, rate, ii, evid, cmt, addl, ss, pMatrix, biovar, tlag,
                              rel_tol, abs_tol, max_num_steps,
                              2e-5, 1e-4, 1e-4, 1e-5);

  TORSTEN_ODE_GRAD_BIOVAR_TEST(pmx_solve_rk45, f, 2,
                              time, amt, rate, ii, evid, cmt, addl, ss, pMatrix, biovar, tlag,
                              rel_tol, abs_tol, max_num_steps,
                              2e-5, 1e-6, 1e-6, 1e-6);

  TORSTEN_ODE_GRAD_BIOVAR_TEST(pmx_solve_bdf, f, 2,
                              time, amt, rate, ii, evid, cmt, addl, ss, pMatrix, biovar, tlag,
                              rel_tol, abs_tol, max_num_steps,
                              2e-5, 1e-6, 1e-4, 1e-6);

  TORSTEN_ODE_GRAD_BIOVAR_TEST(pmx_solve_adams, f, 2,
                              time, amt, rate, ii, evid, cmt, addl, ss, pMatrix, biovar, tlag,
                              rel_tol, abs_tol, max_num_steps,
                              2e-5, 1e-6, 1e-5, 1e-6);
}

TEST_F(TorstenOneCptTest, multiple_bolus_time_dependent_param) {
  pMatrix.resize(nt);
  for (int i = 0; i < nt; i++) {
    pMatrix[i].resize(3);
    if (i < 6) pMatrix[i][0] = 10; // CL
    else pMatrix[i][0] = 50;
    pMatrix[i][1] = 80; // Vc
    pMatrix[i][2] = 1.2; // ka
  }

  time[0] = 0.0;
  for(int i = 1; i < nt; i++) time[i] = time[i - 1] + 2.5;
  addl[0] = 1;

  double rel_tol = 1e-8, abs_tol = 1e-8;
  long int max_num_steps = 1e8;
  MatrixXd x_rk45 = torsten::pmx_solve_rk45(f, nCmt,
                                         time, amt, rate, ii, evid, cmt, addl, ss,
                                         pMatrix, biovar, tlag,
                                         0,
                                         rel_tol, abs_tol, max_num_steps);

  MatrixXd x_bdf = torsten::pmx_solve_bdf(f, nCmt,
                                       time, amt, rate, ii, evid, cmt, addl, ss,
                                       pMatrix, biovar, tlag,
                                       0,
                                       rel_tol, abs_tol, max_num_steps);

  MatrixXd x(nt, 2);
  x << 1000.0, 0.0,
    4.978707e+01, 761.1109513,
    2.478752e+00, 594.7341503,
    1.234098e-01, 437.0034049,
    6.144212e-03, 319.8124495,
    5.488119e+02, 670.0046601,
    2.732374e+01, 323.4948561,
    1.360369e+00, 76.9219400,
    6.772877e-02, 16.5774607,
    3.372017e-03, 3.4974152;
  MatrixXd xt = x.transpose();

  torsten::test::test_val(xt, x_rk45, 1e-6, 1e-5);
  torsten::test::test_val(xt, x_bdf, 1e-4, 1e-5);

  rel_tol = 1e-12;
  abs_tol = 1e-12;

  TORSTEN_ODE_GRAD_THETA_TEST(pmx_solve_rk45, f, nCmt,
                              time, amt, rate, ii, evid, cmt, addl, ss, pMatrix, biovar, tlag,
                              rel_tol, abs_tol, max_num_steps,
                              2e-5, 1e-6, 1e-3, 1e-5);

  TORSTEN_ODE_GRAD_THETA_TEST(pmx_solve_bdf, f, nCmt,
                              time, amt, rate, ii, evid, cmt, addl, ss, pMatrix, biovar, tlag,
                              rel_tol, abs_tol, max_num_steps,
                              1e-5, 1e-6, 5e-3, 1e-4);

  TORSTEN_ODE_GRAD_THETA_TEST(pmx_solve_adams, f, nCmt,
                              time, amt, rate, ii, evid, cmt, addl, ss, pMatrix, biovar, tlag,
                              rel_tol, abs_tol, max_num_steps,
                              2e-5, 1e-4, 1e-3, 1e-4);

  biovar[0] = std::vector<double>{0.8, 0.8};
  TORSTEN_ODE_GRAD_BIOVAR_TEST(pmx_solve_rk45, f, nCmt,
                              time, amt, rate, ii, evid, cmt, addl, ss, pMatrix, biovar, tlag,
                              rel_tol, abs_tol, max_num_steps,
                              2e-5, 1e-6, 1e-4, 1e-6);

  TORSTEN_ODE_GRAD_BIOVAR_TEST(pmx_solve_bdf, f, nCmt,
                              time, amt, rate, ii, evid, cmt, addl, ss, pMatrix, biovar, tlag,
                              rel_tol, abs_tol, max_num_steps,
                              1e-5, 1e-6, 1e-4, 1e-6);

  TORSTEN_ODE_GRAD_BIOVAR_TEST(pmx_solve_adams, f, nCmt,
                              time, amt, rate, ii, evid, cmt, addl, ss, pMatrix, biovar, tlag,
                              rel_tol, abs_tol, max_num_steps,
                              2e-5, 1e-4, 1e-4, 1e-6);
}

TEST_F(TorstenOneCptTest, rate_par) {
  using std::vector;
  using Eigen::Matrix;
  using Eigen::MatrixXd;
  using Eigen::Dynamic;

  amt[0] = 1200;
  rate[0] = 1200;

  double rel_tol = 1e-6, abs_tol = 1e-6;
  long int max_num_steps = 1e6;

  MatrixXd x_rk45 = torsten::pmx_solve_rk45(f, nCmt,
                                  time, amt, rate, ii, evid, cmt, addl, ss,
                                  pMatrix, biovar, tlag,
                                  0,
                                  rel_tol, abs_tol, max_num_steps);

  MatrixXd x_bdf = torsten::pmx_solve_bdf(f, nCmt,
                                time, amt, rate, ii, evid, cmt, addl, ss,
                                pMatrix, biovar, tlag,
                                0,
                                rel_tol, abs_tol, max_num_steps);

  MatrixXd x(nt, 2);
  x << 0.00000,   0.00000,
    259.18178,  40.38605,
    451.18836, 145.61440,
    593.43034, 296.56207,
    698.80579, 479.13371,
    517.68806, 642.57025,
    383.51275, 754.79790,
    284.11323, 829.36134,
    210.47626, 876.28631,
    19.09398, 844.11769;
  MatrixXd xt = x.transpose();
  
  torsten::test::test_val(xt, x_rk45, 1e-6, 1e-5);
  torsten::test::test_val(xt, x_bdf, 1e-5, 1e-5);

  TORSTEN_ODE_GRAD_THETA_TEST(pmx_solve_rk45, f, nCmt,
                              time, amt, rate, ii, evid, cmt, addl, ss, pMatrix, biovar, tlag,
                              rel_tol, abs_tol, max_num_steps,
                              2e-5, 5e-5, 1e-5, 1e-5);

  TORSTEN_ODE_GRAD_THETA_TEST(pmx_solve_bdf, f, nCmt,
                              time, amt, rate, ii, evid, cmt, addl, ss, pMatrix, biovar, tlag,
                              rel_tol, abs_tol, max_num_steps,
                              2e-5, 1e-6, 1e-3, 1e-4);

  TORSTEN_ODE_GRAD_THETA_TEST(pmx_solve_adams, f, nCmt,
                              time, amt, rate, ii, evid, cmt, addl, ss, pMatrix, biovar, tlag,
                              rel_tol, abs_tol, max_num_steps,
                              2e-5, 1e-4, 1e-4, 1e-5);

  TORSTEN_ODE_GRAD_BIOVAR_TEST(pmx_solve_rk45, f, nCmt,
                              time, amt, rate, ii, evid, cmt, addl, ss, pMatrix, biovar, tlag,
                              rel_tol, abs_tol, max_num_steps,
                              2e-5, 5e-5, 1e-6, 1e-6);

  TORSTEN_ODE_GRAD_BIOVAR_TEST(pmx_solve_bdf, f, nCmt,
                              time, amt, rate, ii, evid, cmt, addl, ss, pMatrix, biovar, tlag,
                              rel_tol, abs_tol, max_num_steps,
                              2e-5, 1e-6, 1e-4, 1e-6);

  TORSTEN_ODE_GRAD_BIOVAR_TEST(pmx_solve_adams, f, nCmt,
                              time, amt, rate, ii, evid, cmt, addl, ss, pMatrix, biovar, tlag,
                              rel_tol, abs_tol, max_num_steps,
                              2e-5, 1e-6, 1e-5, 1e-6);
}

TEST_F(TorstenTwoCptTest, rate_par) {
  auto f = refactor::PMXTwoCptModel<double,double,double,double>::f_;

  pMatrix[0][0] = 5;  // CL
  pMatrix[0][1] = 8;  // Q
  pMatrix[0][2] = 35;  // Vc
  pMatrix[0][3] = 105;  // Vp
  pMatrix[0][4] = 1.2;  // ka

  amt[0] = 1200;
  rate[0] = 1200;

  double rel_tol = 1e-6, abs_tol = 1e-6;
  long int max_num_steps = 1e6;

  MatrixXd x_rk45, x_bdf;
  x_rk45 = pmx_solve_rk45(f2, nCmt,
                                         time, amt, rate, ii, evid, cmt, addl, ss,
                                         pMatrix, biovar, tlag,
                                         0,
                                         rel_tol, abs_tol, max_num_steps);
  x_bdf = pmx_solve_bdf(f2, nCmt,
                                       time, amt, rate, ii, evid, cmt, addl, ss,
                                       pMatrix, biovar, tlag,
                                       0,
                                       rel_tol, abs_tol, max_num_steps);

  MatrixXd amounts(10, 3);
  amounts << 0.00000,   0.00000,   0.0000000,
    259.18178,  39.55748,   0.7743944,
    451.18836, 139.65573,   5.6130073,
    593.43034, 278.43884,  17.2109885,
    698.80579, 440.32663,  37.1629388,
    517.68806, 574.76950,  65.5141658,
    383.51275, 653.13596,  99.2568509,
    284.11323, 692.06145, 135.6122367,
    210.47626, 703.65965, 172.6607082,
    19.09398, 486.11014, 406.6342765;
  MatrixXd xt = amounts.transpose();

  // relative error determined empirically
  double rel_err_rk45 = 1e-6, rel_err_bdf = 1e-4;
  expect_near_matrix_eq(xt, x_rk45, rel_err_rk45);
  expect_near_matrix_eq(xt, x_bdf, rel_err_bdf);

  TORSTEN_ODE_GRAD_THETA_TEST(pmx_solve_rk45, f, nCmt,
                              time, amt, rate, ii, evid, cmt, addl, ss, pMatrix, biovar, tlag,
                              rel_tol, abs_tol, max_num_steps,
                              2e-5, 4e-5, 1e-3, 1e-4);

  TORSTEN_ODE_GRAD_THETA_TEST(pmx_solve_bdf, f, nCmt,
                              time, amt, rate, ii, evid, cmt, addl, ss, pMatrix, biovar, tlag,
                              rel_tol, abs_tol, max_num_steps,
                              2e-5, 1e-6, 1e-2, 5e-4);

  TORSTEN_ODE_GRAD_THETA_TEST(pmx_solve_adams, f, nCmt,
                              time, amt, rate, ii, evid, cmt, addl, ss, pMatrix, biovar, tlag,
                              rel_tol, abs_tol, max_num_steps,
                              2e-5, 1e-4, 1e-4, 1e-5);

  TORSTEN_ODE_GRAD_BIOVAR_TEST(pmx_solve_rk45, f, nCmt,
                              time, amt, rate, ii, evid, cmt, addl, ss, pMatrix, biovar, tlag,
                              rel_tol, abs_tol, max_num_steps,
                              2e-5, 5e-5, 1e-6, 1e-6);

  TORSTEN_ODE_GRAD_BIOVAR_TEST(pmx_solve_bdf, f, nCmt,
                              time, amt, rate, ii, evid, cmt, addl, ss, pMatrix, biovar, tlag,
                              rel_tol, abs_tol, max_num_steps,
                              2e-5, 1e-6, 1e-4, 1e-6);

  TORSTEN_ODE_GRAD_BIOVAR_TEST(pmx_solve_adams, f, nCmt,
                              time, amt, rate, ii, evid, cmt, addl, ss, pMatrix, biovar, tlag,
                              rel_tol, abs_tol, max_num_steps,
                              2e-5, 1e-6, 1e-5, 1e-6);
}

TEST_F(FribergKarlssonTest, steady_state) {
  ss[0] = 1;

  double rel_tol = 1e-6, abs_tol = 1e-6;
  long int max_num_steps = 1e6;

  MatrixXd x_rk45, x_bdf;
  x_rk45 = torsten::pmx_solve_rk45(f, nCmt,
                                time, amt, rate, ii, evid, cmt, addl, ss,
                                theta, biovar, tlag,
                                0,
                                rel_tol, abs_tol, max_num_steps);

  x_bdf = torsten::pmx_solve_bdf(f, nCmt,
                              time, amt, rate, ii, evid, cmt, addl, ss,
                              theta, biovar, tlag,
                              0,
                              rel_tol, abs_tol, max_num_steps);

  MatrixXd x(10, 8);
  x << 8.000000e+04, 11996.63, 55694.35, -3.636308, -3.653620,  -3.653933, -3.653748, -3.653622,
    6.566800e+03, 53123.67, 70649.28, -3.650990, -3.653172, -3.653910, -3.653755, -3.653627,
    5.390358e+02, 34202.00, 80161.15, -3.662446, -3.653349, -3.653883, -3.653761, -3.653632,
    4.424675e+01, 23849.69, 80884.40, -3.665321, -3.653782, -3.653870, -3.653765, -3.653637,
    3.631995e+00, 19166.83, 78031.24, -3.664114, -3.654219, -3.653876, -3.653769, -3.653642,
    2.981323e-01, 16799.55, 74020.00, -3.660988, -3.654550, -3.653896, -3.653774, -3.653647,
    2.447219e-02, 15333.26, 69764.65, -3.656791, -3.654722, -3.653926, -3.653779, -3.653653,
    2.008801e-03, 14233.96, 65591.05, -3.651854, -3.654708, -3.653957, -3.653786, -3.653658,
    1.648918e-04, 13303.26, 61607.92, -3.646317, -3.654488, -3.653983, -3.653793, -3.653663,
    1.353552e-05, 12466.56, 57845.10, -3.640244, -3.654050, -3.653995, -3.653801, -3.653668;
  MatrixXd xt = x.transpose();
  torsten::test::test_val(xt, x_rk45, 1.5e-2, 1e-5);
  torsten::test::test_val(xt, x_bdf, 1.5e-2, 1e-5);

  amt[0] = 700;
  ii[0] = 5.0;
  resize(2);
  time[0] = 0;
  for(int i = 1; i < nt; i++) time[i] = time[i - 1] + 5.0;
  addl[0] = 2;
  rel_tol = 1e-12;
  abs_tol = 1e-12;

  TORSTEN_ODE_GRAD_THETA_TEST(pmx_solve_rk45, f, nCmt,
                              time, amt, rate, ii, evid, cmt, addl, ss, theta, biovar, tlag,
                              rel_tol, abs_tol, max_num_steps,
                              2e-5, 1e-6, 1e-3, 1e-4);

  TORSTEN_ODE_GRAD_THETA_TEST(pmx_solve_bdf, f, nCmt,
                              time, amt, rate, ii, evid, cmt, addl, ss, theta, biovar, tlag,
                              rel_tol, abs_tol, max_num_steps,
                              1e-5, 1e-6, 1e-2, 2e-3);

  TORSTEN_ODE_GRAD_THETA_TEST(pmx_solve_adams, f, nCmt,
                              time, amt, rate, ii, evid, cmt, addl, ss, theta, biovar, tlag,
                              rel_tol, abs_tol, max_num_steps,
                              2e-5, 1e-6, 3e-2, 2e-4);
}

TEST_F(FribergKarlssonTest, multiple_bolus) {
  using std::vector;
  using Eigen::Matrix;
  using Eigen::MatrixXd;
  using Eigen::Dynamic;

  double rel_tol = 1e-6, abs_tol = 1e-6;
  long int max_num_steps = 1e6;

  rel_tol = 1e-12;
  abs_tol = 1e-12;
  addl[0] = 2;

  TORSTEN_ODE_GRAD_THETA_TEST(pmx_solve_rk45, f, nCmt,
                              time, amt, rate, ii, evid, cmt, addl, ss, theta, biovar, tlag,
                              rel_tol, abs_tol, max_num_steps,
                              2e-5, 1e-6, 1e-3, 1e-4);

  TORSTEN_ODE_GRAD_THETA_TEST(pmx_solve_bdf, f, nCmt,
                              time, amt, rate, ii, evid, cmt, addl, ss, theta, biovar, tlag,
                              rel_tol, abs_tol, max_num_steps,
                              1e-5, 1e-6, 1e-2, 2e-3);

  TORSTEN_ODE_GRAD_THETA_TEST(pmx_solve_adams, f, nCmt,
                              time, amt, rate, ii, evid, cmt, addl, ss, theta, biovar, tlag,
                              rel_tol, abs_tol, max_num_steps,
                              2e-5, 1e-6, 1e-2, 2e-3);
}

TEST_F(TorstenOdeTest, exception) {
  pMatrix[0][0] = 1.0E-10;
  pMatrix[0][1] = 1.0E-10;
  pMatrix[0][2] = 1.0E-20;
  pMatrix[0][3] = 1.0E+80;
  pMatrix[0][4] = 1.0E+70;

  auto& f = refactor::PMXTwoCptModel<double, double, double, double>::f_;
  int ncmt = refactor::PMXTwoCptModel<double, double, double, double>::Ncmt;

  EXPECT_NO_THROW(torsten::pmx_solve_bdf(f, ncmt, time, amt, rate, ii,
                                               evid, cmt, addl, ss, pMatrix,
                                               biovar, tlag));
}

TEST_F(TorstenOdeTest_neutropenia, max_cvodes_fails) {
  t0 = 12.0;
  ts.resize(1);
  ts[0] = 12.10;
  std::vector<double> y0_bad{80000, 0.004713515595, 0.0002254631631, -1.956259535, -1.956258419, -1.956257852, -1.956257737, -1.956346997};
  std::vector<double> y0_good{80000,0.004717011055, 0.0002256303627, -1.956258109, -1.956258109, -1.956258109, -1.956258109, -1.957624238};
  theta = std::vector<double> {1.508795366, 0.7444735025, 0.9989168227, 0.04373528316, 2.85770127, 0.2803519852, 1.956258109, 0.9692560566, 0.3025591511};

  x_r.resize(y0.size());
  std::fill(x_r.begin(), x_r.end(), 0.0);

  double rtol = 1e-6;
  double atol = 1e-6;
  long int max_num_steps = 1e3;

  y0 = y0_bad;
  EXPECT_THROW(stan::math::integrate_ode_bdf(f, y0, t0, ts, theta, x_r, x_i, 0, rtol, atol, max_num_steps), // NOLINT
               std::runtime_error);
  EXPECT_THROW(torsten::pmx_integrate_ode_bdf(f, y0, t0, ts, theta, x_r, x_i, 0, rtol, atol, max_num_steps), // NOLINT
               std::runtime_error);

  y0 = y0_good;
  EXPECT_NO_THROW(stan::math::integrate_ode_bdf(f, y0, t0, ts, theta, x_r, x_i, 0, rtol, atol, max_num_steps)); // NOLINT
  EXPECT_NO_THROW(torsten::pmx_integrate_ode_bdf(f, y0, t0, ts, theta, x_r, x_i, 0, rtol, atol, max_num_steps)); // NOLINT
}

TEST_F(TorstenOneCptTest, ss_multiple_infusion_rate) {
  time[0] = 0.0;
  for(int i = 1; i < 10; i++) time[i] = time[i - 1] + 8.0;

  amt[0] = 800;
  std::fill(rate.begin(), rate.end(), 150.0);
  ii[0] = 5;
  ii[5] = 2;
  addl[0] = 0;
  ss[0] = 1;
  ss[5] = 1;

  double rel_tol = 1e-8, abs_tol = 1e-8;
  long int max_num_steps = 1e8;
  biovar[0] = std::vector<double>{0.8, 0.7};

  TORSTEN_ODE_GRAD_RATE_TEST(pmx_solve_rk45, f, nCmt,
                             time, amt, rate, ii, evid, cmt, addl, ss, pMatrix, biovar, tlag,
                             rel_tol, abs_tol, max_num_steps,
                             1e-3, 1e-6, 1e-6, 1e-6);

  TORSTEN_ODE_GRAD_RATE_TEST(pmx_solve_bdf, f, nCmt,
                             time, amt, rate, ii, evid, cmt, addl, ss, pMatrix, biovar, tlag,
                             rel_tol, abs_tol, max_num_steps,
                             1e-3, 1e-10, 6e-5, 2e-8);

  TORSTEN_ODE_GRAD_RATE_TEST(pmx_solve_adams, f, nCmt,
                             time, amt, rate, ii, evid, cmt, addl, ss, pMatrix, biovar, tlag,
                             rel_tol, abs_tol, max_num_steps,
                             1e-3, 1e-8, 1e-6, 1e-8);
}
