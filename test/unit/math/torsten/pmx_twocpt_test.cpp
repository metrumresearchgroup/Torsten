#include <test/unit/math/torsten/test_util.hpp>
#include <test/unit/math/torsten/expect_near_matrix_eq.hpp>
#include <test/unit/math/torsten/pmx_twocpt_test_fixture.hpp>
#include <test/unit/math/torsten/pmx_twocpt_mpi_test_fixture.hpp>
#include <stan/math/torsten/pmx_solve_twocpt.hpp>
#include <stan/math/torsten/pmx_solve_bdf.hpp>
#include <stan/math/torsten/pmx_solve_rk45.hpp>
#include <stan/math/torsten/pmx_twocpt_model.hpp>
#include <stan/math/torsten/to_var.hpp>
#include <gtest/gtest.h>
#include <vector>

using std::vector;
using Eigen::Matrix;
using Eigen::MatrixXd;
using Eigen::Dynamic;
using stan::math::var;
using torsten::pmx_solve_twocpt;

TEST_F(TorstenTwoCptTest, single_bolus_tlag) {
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

  TORSTEN_CPT_GRAD_TLAG_TEST(pmx_solve_twocpt, time, amt, rate, ii, evid, cmt, addl, ss, pMatrix, biovar, tlag, 2e-5, 1e-6, 1e-9, 1e-9);
}

TEST_F(TorstenTwoCptTest, multiple_bolus) {
  Matrix<double, Dynamic, Dynamic> x = pmx_solve_twocpt(time, amt, rate, ii, evid, cmt, addl, ss, pMatrix, biovar, tlag);
  Matrix<double, Dynamic, Dynamic> amounts(10, 3);
  amounts << 1000.0, 0.0, 0.0,
    740.818221, 238.3713, 12.75775,
    548.811636, 379.8439, 43.55827,
    406.569660, 455.3096, 83.95657,
    301.194212, 486.6965, 128.32332,
    223.130160, 489.4507, 173.01118,
    165.298888, 474.3491, 215.75441,
    122.456428, 448.8192, 255.23842,
    90.717953, 417.9001, 290.79297,
    8.229747, 200.8720, 441.38985;
  MatrixXd xt = amounts.transpose();
  torsten::test::test_val(xt, x);

  std::vector<std::vector<double>> biovar_test(1, {0.8, 0.9, 0.9});
  std::vector<std::vector<double>> tlag_test(1, {0.4, 0.8, 0.8});
  TORSTEN_CPT_GRAD_THETA_TEST(pmx_solve_twocpt, time, amt, rate, ii, evid, cmt, addl, ss, pMatrix, biovar, tlag, 2e-5, 1e-6, 1e-5, 1e-6);
  TORSTEN_CPT_GRAD_BIOVAR_TEST(pmx_solve_twocpt, time, amt, rate, ii, evid, cmt, addl, ss, pMatrix, biovar_test, tlag, 2e-5, 1e-6, 1e-6, 1e-6);
  TORSTEN_CPT_GRAD_TLAG_TEST(pmx_solve_twocpt, time, amt, rate, ii, evid, cmt, addl, ss, pMatrix, biovar, tlag_test, 2e-5, 1e-6, 1e-6, 1e-6);
}

TEST_F(TorstenTwoCptTest, multiple_bolus_overload) {
  resize(3);
  ii[0] = 4.0;
  addl[0] = 1;
  std::vector<std::vector<double>> biovar_test(1, {0.8, 0.9, 0.9});
  std::vector<std::vector<double>> tlag_test(1, {0.4, 0.8, 0.8});
  TORSTEN_CPT_PARAM_OVERLOAD_TEST(pmx_solve_twocpt, time, amt, rate, ii, evid, cmt, addl, ss, pMatrix, biovar_test, tlag_test, 1e-6, 1e-6);
}

TEST_F(TorstenTwoCptTest, multiple_bolus_central_cmt) {
  cmt[0] = 2;

  std::vector<std::vector<double>> biovar_test(1, {0.8, 0.9, 0.9});
  std::vector<std::vector<double>> tlag_test(1, {0.4, 2.8, 2.8});
  TORSTEN_CPT_GRAD_THETA_TEST(pmx_solve_twocpt, time, amt, rate, ii, evid, cmt, addl, ss, pMatrix, biovar, tlag, 2e-5, 1e-6, 1e-4, 1e-5);
  TORSTEN_CPT_GRAD_BIOVAR_TEST(pmx_solve_twocpt, time, amt, rate, ii, evid, cmt, addl, ss, pMatrix, biovar_test, tlag, 2e-5, 1e-6, 1e-4, 1e-5);
  TORSTEN_CPT_GRAD_TLAG_TEST(pmx_solve_twocpt, time, amt, rate, ii, evid, cmt, addl, ss, pMatrix, biovar, tlag_test, 2e-5, 1e-6, 1e-4, 1e-5);

  using model_t = refactor::PMXTwoCptModel<double, double, double, double>;
  TORSTEN_CPT_ODE_GRAD_TEST(pmx_solve_twocpt, torsten::pmx_solve_bdf, model_t::f_, model_t::Ncmt, 
                            time, amt, rate, ii, evid, cmt, addl, ss, pMatrix, biovar, tlag, 1.E-2, 1.E-2);
}

TEST_F(TorstenTwoCptTest, multiple_addl_iv) {
  cmt[0] = 2;
  rate[0] = 350;
  addl[0] = 2;

  std::vector<std::vector<double> > biovar_test(1, {0.8, 0.9, 0.9});
  std::vector<std::vector<double> > tlag_test(1, {2.3, 2.8, 2.7});
  TORSTEN_CPT_GRAD_THETA_TEST(pmx_solve_twocpt, time, amt, rate, ii, evid, cmt, addl, ss, pMatrix, biovar, tlag, 2e-5, 1e-6, 1e-4, 1e-5);
  TORSTEN_CPT_GRAD_BIOVAR_TEST(pmx_solve_twocpt, time, amt, rate, ii, evid, cmt, addl, ss, pMatrix, biovar_test, tlag, 2e-5, 1e-6, 1e-4, 1e-5);
  TORSTEN_CPT_GRAD_TLAG_TEST(pmx_solve_twocpt, time, amt, rate, ii, evid, cmt, addl, ss, pMatrix, biovar, tlag_test, 2e-5, 1e-6, 1e-4, 1e-5);

  using model_t = refactor::PMXTwoCptModel<double, double, double, double>;
  TORSTEN_CPT_ODE_GRAD_TEST(pmx_solve_twocpt, torsten::pmx_solve_bdf, model_t::f_, model_t::Ncmt, 
                            time, amt, rate, ii, evid, cmt, addl, ss, pMatrix, biovar_test, tlag_test, 1.E-3, 1.E-3);
}

TEST_F(TorstenTwoCptTest, multiple_iv) {
  using model_t = refactor::PMXTwoCptModel<double, double, double, double>;

  rate[0] = 300.0;
  rate[2] = 300.0;
  rate[5] = 400.0;
  addl[0] = 0;

  std::vector<std::vector<double> > biovar_test(1, {0.8, 0.9, 0.9});
  std::vector<std::vector<double> > tlag_test(1, {2.3, 2.8, 2.7});
  TORSTEN_CPT_GRAD_THETA_TEST(pmx_solve_twocpt, time, amt, rate, ii, evid, cmt, addl, ss, pMatrix, biovar, tlag, 2e-5, 1e-6, 1e-4, 1e-5);
  TORSTEN_CPT_GRAD_BIOVAR_TEST(pmx_solve_twocpt, time, amt, rate, ii, evid, cmt, addl, ss, pMatrix, biovar_test, tlag, 2e-5, 1e-6, 1e-4, 1e-5);
  TORSTEN_CPT_GRAD_TLAG_TEST(pmx_solve_twocpt, time, amt, rate, ii, evid, cmt, addl, ss, pMatrix, biovar, tlag_test, 2e-5, 1e-6, 1e-4, 1e-5);

  using model_t = refactor::PMXTwoCptModel<double, double, double, double>;
  TORSTEN_CPT_ODE_GRAD_TEST(pmx_solve_twocpt, torsten::pmx_solve_bdf, model_t::f_, model_t::Ncmt, 
                            time, amt, rate, ii, evid, cmt, addl, ss, pMatrix, biovar_test, tlag_test, 1.E-3, 4.E-4);
}

TEST_F(TorstenPopulationPMXTwoCptTest, multiple_bolus_tlag) {
  tlag[0][0] = 1.7;  // tlag1
  tlag[0][1] = 0;  // tlag2
  tlag[0][2] = 0;  // tlag3

  TORSTEN_CPT_GRAD_THETA_TEST(pmx_solve_twocpt, time, amt, rate, ii, evid, cmt, addl, ss, pMatrix, biovar, tlag, 2e-5, 1e-6, 1e-4, 1e-5);
  TORSTEN_CPT_GRAD_BIOVAR_TEST(pmx_solve_twocpt, time, amt, rate, ii, evid, cmt, addl, ss, pMatrix, biovar, tlag, 2e-5, 1e-6, 1e-4, 1e-5);
  TORSTEN_CPT_GRAD_TLAG_TEST(pmx_solve_twocpt, time, amt, rate, ii, evid, cmt, addl, ss, pMatrix, biovar, tlag, 2e-5, 1e-6, 1e-4, 1e-5);
}

TEST_F(TorstenTwoCptTest, steady_state) {
  time[0] = 0.0;
  time[1] = 0.0;
  for(int i = 2; i < 10; i++) time[i] = time[i - 1] + 5;

  amt[0] = 1200;
  addl[0] = 10;
  ss[0] = 1;

  Matrix<double, Dynamic, Dynamic> x;
  x = torsten::pmx_solve_twocpt(time, amt, rate, ii, evid, cmt, addl, ss,
                             pMatrix, biovar, tlag);

  Matrix<double, Dynamic, Dynamic> amounts(10, 3);
  amounts << 1.200001e+03, 224.5332, 1196.900,
    1.200001e+03, 224.5332, 1196.900, 
    2.974504e+00, 360.2587, 1533.000,
    7.373059e-03, 244.3668, 1294.140,
    3.278849e+01, 548.9136, 1534.479,
    8.127453e-02, 270.3431, 1396.353,
    3.614333e+02, 799.6771, 1304.769,
    8.959035e-01, 316.6314, 1494.540,
    2.220723e-03, 234.0179, 1244.718,
    9.875702e+00, 432.4698, 1552.220;
  MatrixXd xt = amounts.transpose();

  for(int i = 0; i < amounts.rows(); i++) {
    EXPECT_NEAR(xt(0, i), x(0, i), std::max(xt(0, i), x(0, i)) * 1e-6);
    EXPECT_NEAR(xt(1, i), x(1, i), std::max(xt(1, i), x(1, i)) * 1e-6);
  }


  tlag[0][0] = 1.7;  // tlag1
  tlag[0][1] = 0;  // tlag2
  tlag[0][2] = 0;  // tlag3
  TORSTEN_CPT_GRAD_THETA_TEST(pmx_solve_twocpt, time, amt, rate, ii, evid, cmt, addl, ss, pMatrix, biovar, tlag, 2e-5, 1e-6, 1e-4, 1e-5);
  TORSTEN_CPT_GRAD_BIOVAR_TEST(pmx_solve_twocpt, time, amt, rate, ii, evid, cmt, addl, ss, pMatrix, biovar, tlag, 2e-5, 1e-6, 1e-4, 1e-5);
  TORSTEN_CPT_GRAD_TLAG_TEST(pmx_solve_twocpt, time, amt, rate, ii, evid, cmt, addl, ss, pMatrix, biovar, tlag, 2e-5, 1e-6, 1e-4, 1e-5);
}

TEST_F(TorstenTwoCptTest, steady_state_overload) {
  resize(3);
  time[0] = 0.0;
  time[1] = 0.0;
  for(int i = 2; i < nt; i++) time[i] = time[i - 1] + 5;

  amt[0] = 1200;
  addl[0] = 10;
  ss[0] = 1;

  tlag[0][0] = 1.7;  // tlag1
  tlag[0][1] = 0;  // tlag2
  tlag[0][2] = 0;  // tlag3

  TORSTEN_CPT_PARAM_OVERLOAD_TEST(pmx_solve_twocpt, time, amt, rate, ii, evid, cmt, addl, ss, pMatrix, biovar, tlag, 1e-6, 1e-6);  
}

TEST_F(TorstenTwoCptTest, multiple_steady_state_iv) {
  time[0] = 0.0;
  time[1] = 0.0;
  for(int i = 2; i < 10; i++) time[i] = time[i - 1] + 5;
  amt[0] = 1200;
  rate[0] = 150;
  addl[0] = 10;
  ss[0] = 1;

  Matrix<double, Dynamic, Dynamic> x;
  x = torsten::pmx_solve_twocpt(time, amt, rate, ii, evid, cmt, addl, ss,
                             pMatrix, biovar, tlag);

  Matrix<double, Dynamic, Dynamic> amounts(10, 3);
  amounts << 1.028649, 286.5656, 1391.610,
    1.028649, 286.5656, 1391.610,
    124.692706, 452.4021, 1377.667,
    11.338982, 367.1773, 1461.416,
    121.612641, 410.2024, 1340.203,
    124.991604, 477.3286, 1452.499,
    87.660547, 315.1768, 1352.746,
    124.907445, 463.2236, 1402.095,
    3.415236, 318.7214, 1432.451,
    123.979747, 436.1057, 1355.890;
  MatrixXd xt = amounts.transpose();

  for(int i = 0; i < amounts.rows(); i++) {
    EXPECT_NEAR(xt(0, i), x(0, i), std::max(xt(0, i), x(0, i)) * 1e-6);
    EXPECT_NEAR(xt(1, i), x(1, i), std::max(xt(1, i), x(1, i)) * 1e-6);
  }

  std::vector<std::vector<double> > biovar_test(1, {0.8, 0.9, 0.7});
  std::vector<std::vector<double> > tlag_test(1, {2.4, 1.7, 1.7});
  TORSTEN_CPT_GRAD_THETA_TEST(pmx_solve_twocpt, time, amt, rate, ii, evid, cmt, addl, ss, pMatrix, biovar, tlag, 2e-5, 1e-6, 1e-4, 1e-5);
  TORSTEN_CPT_GRAD_BIOVAR_TEST(pmx_solve_twocpt, time, amt, rate, ii, evid, cmt, addl, ss, pMatrix, biovar_test, tlag, 2e-5, 1e-6, 1e-4, 1e-5);
  TORSTEN_CPT_GRAD_TLAG_TEST(pmx_solve_twocpt, time, amt, rate, ii, evid, cmt, addl, ss, pMatrix, biovar, tlag_test, 2e-5, 1e-6, 1e-4, 1e-5);
}

TEST_F(TorstenTwoCptTest, multiple_steady_state_iv_overload) {
  resize(3);
  time[0] = 0.0;
  time[1] = 0.0;
  for(int i = 2; i < nt; i++) time[i] = time[i - 1] + 5;
  addl[0] = 1;
  amt[0] = 1200;
  rate[0] = 150;
  addl[0] = 10;
  ss[0] = 1;

  std::vector<std::vector<double> > biovar_test(1, {0.8, 0.9, 0.7});
  std::vector<std::vector<double> > tlag_test(1, {2.4, 1.7, 1.7});
  TORSTEN_CPT_PARAM_OVERLOAD_TEST(pmx_solve_twocpt, time, amt, rate, ii, evid, cmt, addl, ss,
                                  pMatrix, biovar, tlag, 1e-6, 1e-6);  
}

TEST_F(TorstenTwoCptTest, events_specific_data) {
  nt = 11;
  pMatrix.resize(nt);
  for (int i = 0; i < nt; i++) {
    pMatrix[i].resize(5);
    if (i < 6) pMatrix[i][0] = 5; // CL
    else pMatrix[i][0] = 50;  // CL is piece-wise constant
    pMatrix[i][1] = 8;  // Q
    pMatrix[i][2] = 20;  // Vc
    pMatrix[i][3] = 70;  // Vp
    pMatrix[i][4] = 1.2;  // ka
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

  Matrix<double, Dynamic, Dynamic> x = torsten::pmx_solve_twocpt(time, amt, rate, ii, evid, cmt, addl, ss,
                                                              pMatrix, biovar, tlag);

  Matrix<double, Dynamic, Dynamic> amounts(nt, 3);
  amounts << 1.000000e+03,   0.000000,   0.0000,
    4.978707e+01, 352.089056, 349.4148,
    2.478752e+00, 146.871246, 458.3010,
    1.234098e-01,  93.537648, 442.6420,
    6.144212e-03,  77.732083, 405.7800,
    5.488119e+02, 449.105589, 412.0337,
    2.732374e+01,  36.675537, 430.0023,
    1.360369e+00,  14.886990, 341.6754,
    6.772877e-02,  10.966107, 267.7033,
    3.372017e-03,   8.549649, 209.5604,
    1.678828e-04,   6.690631, 164.0364;
  MatrixXd xt = amounts.transpose();
  torsten::test::test_val(xt, x);

  TORSTEN_CPT_GRAD_THETA_TEST(pmx_solve_twocpt, time, amt, rate, ii, evid, cmt, addl, ss, pMatrix, biovar, tlag, 2e-5, 1e-6, 1e-4, 1e-5);
}

TEST_F(TorstenTwoCptTest, multiple_iv_var) {
  using std::vector;

  pMatrix[0][0] = 5;  // CL
  pMatrix[0][1] = 8;  // Q
  pMatrix[0][2] = 35;  // Vc
  pMatrix[0][3] = 105;  // Vp
  pMatrix[0][4] = 1.2;  // ka
  
  amt[0] = 1200;
  rate[0] = 1200;

  Matrix<double, Dynamic, Dynamic> x = torsten::pmx_solve_twocpt(time, amt, rate, ii, evid, cmt, addl, ss,
                                                              pMatrix, biovar, tlag);

  Matrix<double, Dynamic, Dynamic> amounts(10, 3);
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
  torsten::test::test_val(xt, x);

  std::fill(rate.begin(), rate.end(), 340);
  TORSTEN_CPT_GRAD_RATE_TEST(pmx_solve_twocpt, time, amt, rate, ii, evid, cmt, addl, ss, pMatrix, biovar, tlag, 2e-5, 1e-6, 1e-6, 1e-6);
}

TEST_F(TorstenTwoCptTest, multiple_iv_var_overload) {
  resize(3);
  time[0] = 0.0;
  time[1] = 0.0;
  for(int i = 2; i < nt; i++) time[i] = time[i - 1] + 5;

  pMatrix[0][0] = 5;  // CL
  pMatrix[0][1] = 8;  // Q
  pMatrix[0][2] = 35;  // Vc
  pMatrix[0][3] = 105;  // Vp
  pMatrix[0][4] = 1.2;  // ka
  
  amt[0] = 1200;
  rate[0] = 340;
  std::vector<stan::math::var> rate_v(stan::math::to_var(rate));
  TORSTEN_CPT_PARAM_OVERLOAD_TEST(pmx_solve_twocpt, time, amt, rate_v, ii, evid, cmt, addl, ss,
                                  pMatrix, biovar, tlag, 1e-6, 1e-6);
}

TEST_F(TorstenTwoCptTest, single_iv_central_cmt_var) {
  using std::vector;
  using stan::math::var;

  pMatrix[0][2] = 35;  // Vc
  pMatrix[0][3] = 105;  // Vp

  amt[0] = 1200;
  rate[0] = 700;
  cmt[0] = 2;

  Matrix<double, Dynamic, Dynamic> x = torsten::pmx_solve_twocpt(time, amt, rate, ii, evid, cmt, addl, ss,
                                                              pMatrix, biovar, tlag);
  Matrix<double, Dynamic, Dynamic> amounts(10, 3);
  amounts << 0.00000,   0.00000,   0.0000000,
    0, 167.150925, 4.81832410,
    0, 319.651326,  18.583668,
    0, 458.954120, 40.3399483,
    0, 586.366896, 69.2271570,
    0, 703.066462,  104.47175,
    0, 810.111927, 145.377981,
    0, 883.621482, 191.218629,
    0, 809.159915, 235.475039,
    0, 425.085221, 450.714833;
  MatrixXd xt = amounts.transpose();
  torsten::test::test_val(xt, x);

  std::fill(rate.begin(), rate.end(), 780);
  TORSTEN_CPT_GRAD_RATE_TEST(pmx_solve_twocpt, time, amt, rate, ii, evid, cmt, addl, ss, pMatrix, biovar, tlag, 2e-5, 1e-6, 1e-6, 1e-6);
}

TEST_F(TorstenTwoCptTest, single_iv_central_cmt_var_overload) {
  resize(3);
  time[0] = 0.0;
  time[1] = 0.0;
  for(int i = 2; i < nt; i++) time[i] = time[i - 1] + 5;

  pMatrix[0][2] = 35;  // Vc
  pMatrix[0][3] = 105;  // Vp

  amt[0] = 1200;
  cmt[0] = 2;
  rate[0] = 780;
  std::vector<stan::math::var> rate_v(stan::math::to_var(rate));
  TORSTEN_CPT_PARAM_OVERLOAD_TEST(pmx_solve_twocpt, time, amt, rate_v, ii, evid, cmt, addl, ss,
                                  pMatrix, biovar, tlag, 1e-6, 1e-6);
}
