#include <test/unit/math/torsten/test_util.hpp>
#include <gtest/gtest.h>
#include <test/unit/math/torsten/pmx_onecpt_test_fixture.hpp>
#include <stan/math/torsten/pmx_solve_onecpt.hpp>
#include <stan/math/torsten/pmx_onecpt_model.hpp>
#include <stan/math/torsten/pmx_twocpt_model.hpp>
#include <gtest/gtest.h>
#include <vector>

using std::vector;
using Eigen::Matrix;
using Eigen::MatrixXd;
using Eigen::Dynamic;
using torsten::pmx_solve_onecpt;

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

  TORSTEN_CPT_GRAD_TLAG_TEST(pmx_solve_onecpt, time, amt, rate, ii, evid, cmt, addl, ss, pMatrix, biovar, tlag, 2e-5, 1e-6, 1e-6, 1e-6);
}

TEST_F(TorstenOneCptTest, multiple_bolus) {
  MatrixXd x = pmx_solve_onecpt(time, amt, rate, ii, evid, cmt, addl, ss, pMatrix, biovar, tlag);
  MatrixXd amounts(10, 2);
  amounts << 1000.0, 0.0,
    740.8182, 254.97490,
    548.8116, 436.02020,
    406.5697, 562.53846,
    301.1942, 648.89603,
    223.1302, 705.72856,
    165.2989, 740.90816,
    122.4564, 760.25988,
    90.71795, 768.09246,
    8.229747, 667.87079;
  MatrixXd xt = amounts.transpose();
  torsten::test::test_val(xt, x);

  std::vector<std::vector<double>> biovar_test(1, {0.8, 0.9});
  std::vector<std::vector<double>> tlag_test(1, {0.4, 0.5});
  TORSTEN_CPT_GRAD_THETA_TEST(pmx_solve_onecpt, time, amt, rate, ii, evid, cmt, addl, ss, pMatrix, biovar, tlag, 2e-5, 1e-6, 1e-4, 1e-5);
  TORSTEN_CPT_GRAD_BIOVAR_TEST(pmx_solve_onecpt, time, amt, rate, ii, evid, cmt, addl, ss, pMatrix, biovar_test, tlag, 2e-5, 1e-6, 1e-4, 1e-5);
  TORSTEN_CPT_GRAD_TLAG_TEST(pmx_solve_onecpt, time, amt, rate, ii, evid, cmt, addl, ss, pMatrix, biovar, tlag_test, 2e-5, 1e-6, 1e-4, 1e-5);
}

TEST_F(TorstenOneCptTest, multiple_bolus_overload) {
  resize(4);
  time[0] = 0;
  for(int i = 1; i < nt; i++) time[i] = time[i - 1] + 0.9;
  addl[0] = 1;

  biovar[0] = std::vector<double>{0.8, 0.9};
  tlag[0] = std::vector<double>{0.5, 0.8};
  TORSTEN_CPT_PARAM_OVERLOAD_TEST(pmx_solve_onecpt, time, amt, rate, ii, evid, cmt, addl, ss,
                                  pMatrix, biovar, tlag, 1e-6, 1e-6);
}

TEST_F(TorstenOneCptTest, multiple_bolus_central_cmt) {
  cmt[0] = 2;
  MatrixXd x;
  x = pmx_solve_onecpt(time, amt, rate, ii, evid, cmt, addl, ss, pMatrix, biovar, tlag);

  MatrixXd amounts(10, 2);
  amounts << 0, 1000.0000,
    0,  969.2332,
    0,  939.4131,
    0 , 910.5104,
    0,  882.4969,
    0,  855.3453,
    0,  829.0291,
    0,  803.5226,
    0,  778.8008,
    0,  606.5307;
  MatrixXd xt = amounts.transpose();
  torsten::test::test_val(xt, x);

  std::vector<std::vector<double> > biovar_test(1, {0.8, 0.9});
  std::vector<std::vector<double> > tlag_test(1, {0.4, 0.8});
  TORSTEN_CPT_GRAD_THETA_TEST(pmx_solve_onecpt, time, amt, rate, ii, evid, cmt, addl, ss, pMatrix, biovar, tlag, 2e-5, 1e-6, 1e-4, 1e-5);
  TORSTEN_CPT_GRAD_BIOVAR_TEST(pmx_solve_onecpt, time, amt, rate, ii, evid, cmt, addl, ss, pMatrix, biovar_test, tlag, 2e-5, 1e-6, 1e-4, 1e-5);
  TORSTEN_CPT_GRAD_TLAG_TEST(pmx_solve_onecpt, time, amt, rate, ii, evid, cmt, addl, ss, pMatrix, biovar, tlag_test, 2e-5, 1e-6, 1e-4, 1e-5);

  using model_t = refactor::PMXOneCptModel<double, double, double, double>;
  TORSTEN_CPT_ODE_GRAD_TEST(pmx_solve_onecpt, torsten::pmx_solve_bdf, model_t::f_, model_t::Ncmt, 
                            time, amt, rate, ii, evid, cmt, addl, ss, pMatrix, biovar, tlag, 1.E-2, 1.E-2);
}

TEST_F(TorstenOneCptTest, multiple_iv) {
  cmt[0] = 2;
  rate[0] = 350;
  addl[0] = 2;

  std::vector<std::vector<double> > biovar_test(1, {0.8, 0.9});
  std::vector<std::vector<double> > tlag_test(1, {0.4, 0.8});
  TORSTEN_CPT_GRAD_THETA_TEST(pmx_solve_onecpt, time, amt, rate, ii, evid, cmt, addl, ss, pMatrix, biovar, tlag, 2e-5, 1e-6, 1e-4, 1e-5);
  TORSTEN_CPT_GRAD_BIOVAR_TEST(pmx_solve_onecpt, time, amt, rate, ii, evid, cmt, addl, ss, pMatrix, biovar_test, tlag, 2e-5, 1e-6, 1e-4, 1e-5);
  TORSTEN_CPT_GRAD_TLAG_TEST(pmx_solve_onecpt, time, amt, rate, ii, evid, cmt, addl, ss, pMatrix, biovar, tlag_test, 2e-5, 1e-6, 1e-4, 1e-5);

  using model_t = refactor::PMXOneCptModel<double, double, double, double>;
  TORSTEN_CPT_ODE_GRAD_TEST(pmx_solve_onecpt, torsten::pmx_solve_bdf, model_t::f_, model_t::Ncmt, 
                            time, amt, rate, ii, evid, cmt, addl, ss, pMatrix, biovar_test, tlag_test, 5.E-3, 1.E-2);
}

TEST_F(TorstenOneCptTest, multiple_bolus_tlag) {
  tlag[0].resize(nCmt);
  tlag[0][0] = 5;  // tlag1
  tlag[0][1] = 0;  // tlag2

  time[0] = 0;
  for(int i = 1; i < 10; i++) time[i] = time[i - 1] + 1;

  amt[0] = 1200;
  
  MatrixXd x = pmx_solve_onecpt(time, amt, rate, ii, evid, cmt, addl, ss, pMatrix, biovar, tlag);
  
  MatrixXd amounts(10, 2);
  amounts << 0, 0,
             0, 0,
             0, 0,
             0, 0,
             0, 0,
             0, 0,
             361.433054, 778.6752,
             108.861544, 921.7110,
             32.788467, 884.0469,
             9.875696, 801.4449;
  MatrixXd xt = amounts.transpose();
  torsten::test::test_val(xt, x);

  TORSTEN_CPT_GRAD_THETA_TEST(pmx_solve_onecpt, time, amt, rate, ii, evid, cmt, addl, ss, pMatrix, biovar, tlag, 2e-5, 1e-6, 1e-4, 1e-5);
  TORSTEN_CPT_GRAD_BIOVAR_TEST(pmx_solve_onecpt, time, amt, rate, ii, evid, cmt, addl, ss, pMatrix, biovar, tlag, 2e-5, 1e-6, 1e-4, 1e-5);

  tlag[0][0] = 1.7;  // tlag1
  tlag[0][1] = 0;  // tlag2
  TORSTEN_CPT_GRAD_TLAG_TEST(pmx_solve_onecpt, time, amt, rate, ii, evid, cmt, addl, ss, pMatrix, biovar, tlag, 2e-5, 1e-6, 1e-4, 1e-5);
}

TEST_F(TorstenOneCptTest, steady_state) {
  time[0] = 0.0;
  time[1] = 0.0;
  for(int i = 2; i < 10; i++) time[i] = time[i - 1] + 5;
  amt[0] = 1200;
  addl[0] = 10;
  ss[0] = 1;

  MatrixXd x;
  x = pmx_solve_onecpt(time, amt, rate, ii, evid, cmt, addl, ss,
                             pMatrix, biovar, tlag);

  MatrixXd amounts(10, 2);
  amounts << 1200.0, 384.7363,
    1200.0, 384.7363,
    2.974504, 919.6159,
    7.373062e-3, 494.0040,
    3.278849e+1, 1148.4725,
    8.127454e-2, 634.2335,
    3.614333e+2, 1118.2043,
    8.959035e-1, 813.4883,
    2.220724e-3, 435.9617,
    9.875702, 1034.7998;
  MatrixXd xt = amounts.transpose();

  for(int i = 0; i < amounts.rows(); i++) {
    EXPECT_NEAR(xt(0, i), x(0, i), std::max(xt(0, i), x(0, i)) * 1e-6);
    EXPECT_NEAR(xt(1, i), x(1, i), std::max(xt(1, i), x(1, i)) * 1e-6);
  }

  TORSTEN_CPT_GRAD_THETA_TEST(pmx_solve_onecpt, time, amt, rate, ii, evid, cmt, addl, ss, pMatrix, biovar, tlag, 2e-5, 1e-6, 1e-4, 1e-5);
  TORSTEN_CPT_GRAD_BIOVAR_TEST(pmx_solve_onecpt, time, amt, rate, ii, evid, cmt, addl, ss, pMatrix, biovar, tlag, 2e-5, 1e-6, 1e-4, 1e-5);

  tlag[0][0] = 1.7;  // tlag1
  tlag[0][1] = 0;  // tlag2
  TORSTEN_CPT_GRAD_TLAG_TEST(pmx_solve_onecpt, time, amt, rate, ii, evid, cmt, addl, ss, pMatrix, biovar, tlag, 2e-5, 1e-6, 1e-4, 1e-5);
}

TEST_F(TorstenOneCptTest, steady_state_overload) {
  resize(3);
  time[0] = 0.0;
  time[1] = 0.0;
  for(int i = 2; i < nt; i++) time[i] = time[i - 1] + 5;
  amt[0] = 1200;
  addl[0] = 10;
  ss[0] = 1;

  tlag[0][0] = 1.7;  // tlag1
  tlag[0][1] = 0;  // tlag2
  TORSTEN_CPT_PARAM_OVERLOAD_TEST(pmx_solve_onecpt, time, amt, rate, ii, evid, cmt, addl, ss,
                                  pMatrix, biovar, tlag, 1e-6, 1e-6);
}

TEST_F(TorstenOneCptTest, steady_state_multiple_infusion) {
  time[0] = 0.0;
  time[1] = 0.0;
  for(int i = 2; i < 10; i++) time[i] = time[i - 1] + 5;
  amt[0] = 1200;
  rate[0] = 150;
  addl[0] = 10;
  ss[0] = 1;
  MatrixXd x = pmx_solve_onecpt(time, amt, rate, ii, evid, cmt, addl, ss, pMatrix, biovar, tlag);

  MatrixXd amounts(10, 2);
  amounts << 1.028649, 659.9385,
    1.028649, 659.9385,
    124.692706, 837.1959,
    11.338982, 836.1947,
    121.612641, 737.4911,
    124.991604, 950.4222,
    87.660547, 642.9529,
    124.907445, 879.6271,
    3.415236, 745.2971,
    123.979747, 789.6393;
  MatrixXd xt = amounts.transpose();

  torsten::test::test_val(x, xt, 1e-6, 1e-8);

  std::vector<std::vector<double> > biovar_test(1, {0.8, 0.9});
  std::vector<std::vector<double> > tlag_test(1, {2.4, 1.7});
  TORSTEN_CPT_GRAD_THETA_TEST(pmx_solve_onecpt, time, amt, rate, ii, evid, cmt, addl, ss, pMatrix, biovar, tlag, 2e-5, 1e-6, 1e-4, 1e-5);
  TORSTEN_CPT_GRAD_BIOVAR_TEST(pmx_solve_onecpt, time, amt, rate, ii, evid, cmt, addl, ss, pMatrix, biovar_test, tlag, 2e-5, 1e-6, 1e-4, 1e-5);
  TORSTEN_CPT_GRAD_TLAG_TEST(pmx_solve_onecpt, time, amt, rate, ii, evid, cmt, addl, ss, pMatrix, biovar, tlag_test, 2e-5, 1e-6, 1e-4, 1e-5);
}

TEST_F(TorstenOneCptTest, steady_state_multiple_infusion_vs_ode) {
  resize(3);
  amt[0] = 1200;
  rate[0] = 100;
  addl[0] = 0;
  ii[0] = amt[0]/rate[0] + 5.0;
  ss[0] = 1;
  time[0] = 0.0;
  time[1] = ii[0] * 0.5;
  time[2] = ii[0];

  MatrixXd x = pmx_solve_onecpt(time, amt, rate, ii, evid, cmt, addl, ss, pMatrix, biovar, tlag);

  // compare with long-run numerical integration solution
  resize(2);
  amt[0] = 1200;
  rate[0] = 100;
  addl[0] = 20;
  ii[0] = amt[0]/rate[0] + 5.0;
  ss[0] = 0;
  time[0] = 0.0;
  time[1] = addl[0] * ii[0];

  auto f = refactor::PMXOneCptModel<double, double, double, double>::f_;
  MatrixXd x_ode = torsten::pmx_solve_rk45(f, nCmt, time, amt, rate, ii, evid, cmt, addl, ss, pMatrix, biovar, tlag);

  // The last columns from the two runs are the steady state results
  for (int i = 0; i < nCmt; ++i) {
    EXPECT_NEAR(x.rightCols(1)(i), x_ode.rightCols(1)(i), 1.e-5);
  }
}

TEST_F(TorstenOneCptTest, multiple_steady_state_iv_overload) {
  resize(3);
  time[0] = 0.0;
  time[1] = 0.0;
  for(int i = 2; i < nt; i++) time[i] = time[i - 1] + 5;
  amt[0] = 1200;
  rate[0] = 150;
  addl[0] = 10;
  ss[0] = 1;
  std::vector<std::vector<double> > biovar_test(1, {0.8, 0.9});
  std::vector<std::vector<double> > tlag_test(1, {2.4, 1.7});
  TORSTEN_CPT_PARAM_OVERLOAD_TEST(pmx_solve_onecpt, time, amt, rate, ii, evid, cmt, addl, ss,
                                  pMatrix, biovar, tlag, 1e-6, 1e-6);
}

TEST_F(TorstenOneCptTest, events_specific_data) {

  nt = 11;
  pMatrix.resize(nt);
  for (int i = 0; i < nt; i++) {
    pMatrix[i].resize(3);
    if (i < 6) pMatrix[i][0] = 10; // CL
    else pMatrix[i][0] = 50; // CL is piece-wise contant
    pMatrix[i][1] = 80; // Vc
    pMatrix[i][2] = 1.2; // ka
  }

  time.resize(nt);
  time[0] = 0.0;
  for(int i = 1; i < nt; i++) time[i] = time[i - 1] + 2.5;

  amt.resize(nt);
  amt[0] = 1000;
  std::fill(amt.begin() + 1, amt.end(), 0.0);

  rate.resize(nt);
  std::fill(rate.begin(), rate.end(), 0.0);

  cmt.resize(nt);
  cmt[0] = 1;
  std::fill(cmt.begin() + 1, cmt.end(), 2);

  evid.resize(nt);
  evid[0] = 1;
  std::fill(evid.begin() + 1, evid.end(), 0);

  ii.resize(nt);
  ii[0] = 12;
  std::fill(ii.begin() + 1, ii.end(), 0);

  addl.resize(nt);
  addl[0] = 1;
  std::fill(addl.begin() + 1, addl.end(), 0);

  ss.resize(nt);
  std::fill(ss.begin(), ss.end(), 0);

  MatrixXd x = pmx_solve_onecpt(time, amt, rate, ii, evid, cmt, addl, ss,
                                                     pMatrix, biovar, tlag);
  MatrixXd amounts(nt, 2);
  amounts << 1000.0, 0.0,
    4.978707e+01, 761.1109513,
    2.478752e+00, 594.7341503,
    1.234098e-01, 437.0034049,
    6.144212e-03, 319.8124495,
    5.488119e+02, 670.0046601,
    2.732374e+01, 323.4948561,
    1.360369e+00, 76.9219400,
    6.772877e-02, 16.5774607,
    3.372017e-03, 3.4974152,
    1.678828e-04, 0.7342228;
  MatrixXd xt = amounts.transpose();
  torsten::test::test_val(xt, x);

  TORSTEN_CPT_GRAD_THETA_TEST(pmx_solve_onecpt, time, amt, rate, ii, evid, cmt, addl, ss, pMatrix, biovar, tlag, 2e-5, 1e-6, 1e-4, 1e-5);
}

TEST_F(TorstenOneCptTest, single_iv_var) {
  using std::vector;
  amt[0] = 1200;
  rate[0] = 1200;

  MatrixXd x = pmx_solve_onecpt(time, amt, rate, ii, evid, cmt, addl, ss,
                                                     pMatrix, biovar, tlag);

  MatrixXd amounts(10, 2);
  amounts << 0.00000,   0.00000,
             259.18178,  40.38605,
             451.18836, 145.61440,
             593.43034, 296.56207,
             698.80579, 479.13371,
             517.68806, 642.57025,
             383.51275, 754.79790,
             284.11323, 829.36134,
             210.47626, 876.28631,
             19.09398, 844.11769;
  MatrixXd xt = amounts.transpose();
  torsten::test::test_val(xt, x);

  std::fill(rate.begin(), rate.end(), 340.0);
  TORSTEN_CPT_GRAD_RATE_TEST(pmx_solve_onecpt, time, amt, rate, ii, evid, cmt, addl, ss, pMatrix, biovar, tlag, 2e-5, 1e-6, 1e-6, 1e-6);
}

TEST_F(TorstenOneCptTest, single_iv_var_overload) {
  resize(2);
  amt[0] = 1200;
  rate[0] = 340;
  std::vector<stan::math::var> rate_v(stan::math::to_var(rate));
  TORSTEN_CPT_PARAM_OVERLOAD_TEST(pmx_solve_onecpt, time, amt, rate_v, ii, evid, cmt, addl, ss,
                                  pMatrix, biovar, tlag, 1e-6, 1e-6);
}

TEST_F(TorstenOneCptTest, single_iv_central_cmt_var) {
  using stan::math::var;
  cmt[0] = 2;  // IV infusion, not absorption from the gut
  rate[0] = 600;

  MatrixXd x = pmx_solve_onecpt(time, amt, rate, ii, evid, cmt, addl, ss,
                                                     pMatrix, biovar, tlag);

  MatrixXd amounts(10, 2);
  amounts << 
    0,           0,
    0, 147.6804745,
    0, 290.8172984,
    0, 429.5502653,
    0, 564.0148675,
    0, 694.3424289,
    0, 820.6602327,
    0, 893.3511610,
    0,  865.865635,
    0, 674.3368348;
  MatrixXd xt = amounts.transpose();
  torsten::test::test_val(xt, x);

  std::fill(rate.begin(), rate.end(), 660.);
  TORSTEN_CPT_GRAD_RATE_TEST(pmx_solve_onecpt, time, amt, rate, ii, evid, cmt, addl, ss, pMatrix, biovar, tlag, 2e-5, 1e-6, 1e-6, 1e-6);
}

TEST_F(TorstenOneCptTest, single_iv_central_cmt_var_overload) {
  resize(3);
  cmt[0] = 2;  // IV infusion, not absorption from the gut
  rate[0] = 660;
  std::vector<stan::math::var> rate_v(stan::math::to_var(rate));
  TORSTEN_CPT_PARAM_OVERLOAD_TEST(pmx_solve_onecpt, time, amt, rate_v, ii, evid, cmt, addl, ss,
                                  pMatrix, biovar, tlag, 1e-6, 1e-6);
}

/*
TEST(Torsten, pmx_solve_onecptModel_SS_rate_2) {
  // Test the special case where the infusion rate is longer than
  // the interdose interval.
  // THIS TEST FAILS.
  using std::vector;

  vector<vector<double> > pMatrix(1);
  pMatrix[0].resize(3);
  pMatrix[0][0] = 10;  // CL
  pMatrix[0][1] = 80;  // Vc
  pMatrix[0][2] = 1.2;  // ka

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
  time[0] = 0;
  for(int i = 1; i < 9; i++) time[i] = time[i - 1] + 0.25;
  time[9] = 4.0;

  vector<double> amt(10, 0);
  amt[0] = 1200;

  vector<double> rate(10, 0);
  rate[0] = 75;

  vector<int> cmt(10, 2);
  cmt[0] = 1;

  vector<int> evid(10, 0);
  evid[0] = 1;

  vector<double> ii(10, 0);
  ii[0] = 12;

  vector<int> addl(10, 0);
  addl[0] = 14;

  vector<int> ss(10, 0);

  MatrixXd x;
  x = pmx_solve_onecpt(time, amt, rate, ii, evid, cmt, addl, ss,
                    pMatrix, biovar, tlag);

  std::cout << x << std::endl;

  MatrixXd amounts(10, 2);
  amounts << 62.50420, 724.7889,
             78.70197, 723.4747,
             90.70158, 726.3310,
             99.59110, 732.1591,
             106.17663, 740.0744,
             111.05530, 749.4253,
             114.66951, 759.7325,
             117.34699, 770.6441,
             119.33051, 781.9027,
             124.48568, 870.0308;

  // expect_matrix_eq(amounts, x);

  // Test AutoDiff against FiniteDiff
  // test_pmx_solve_onecpt(time, amt, rate, ii, evid, cmt, addl, ss,
  //                   pMatrix, biovar, tlag, 1e-8, 5e-4);
} */
