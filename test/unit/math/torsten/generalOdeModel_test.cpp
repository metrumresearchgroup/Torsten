#include <stan/math/rev/mat.hpp>  // FIX ME - includes should be more specific
#include <gtest/gtest.h>
#include <test/unit/math/prim/mat/fun/expect_near_matrix_eq.hpp>
#include <test/unit/math/prim/mat/fun/expect_matrix_eq.hpp>
#include <test/unit/math/torsten/util_generalOdeModel.hpp>

  // Developer's note: for the autodiff test, the rk45 agrees
  // more closely with finite diff than bdf by an order of
  // magnitude.
  template <typename T0, typename T1, typename T2, typename T3>
  inline
  std::vector<typename boost::math::tools::promote_args<T0, T1, T2, T3>::type>
  oneCptModelODE(const T0& t,
                 const std::vector<T1>& x,
                 const std::vector<T2>& parms,
                 const std::vector<T3>& rate,
                 const std::vector<int>& dummy, std::ostream* pstream__) {
    typedef typename boost::math::tools::promote_args<T0, T1, T2, T3>::type scalar;

    scalar CL = parms[0], V1 = parms[1], ka = parms[2], k10 = CL / V1;
    std::vector<scalar> y(2, 0);

    y[0] = -ka * x[0];
    y[1] = ka * x[0] - k10 * x[1];

    return y;
  }

  struct oneCptModelODE_functor {
    template <typename T0, typename T1, typename T2, typename T3>
    inline
    std::vector<typename boost::math::tools::promote_args<T0, T1, T2, T3>::type>
    operator()(const T0& t,
               const std::vector<T1>& x,
               const std::vector<T2>& parms,
               const std::vector<T3>& rate,
               const std::vector<int>& dummy, std::ostream* pstream__) const {
      return oneCptModelODE(t, x, parms, rate, dummy, pstream__);
    }
  };

class TorstenPKOneCptODETest : public testing::Test {
  void SetUp() {
    // make sure memory's clean before starting each test
    stan::math::recover_memory();
  }
  
public:
  TorstenPKOneCptODETest() :
  pMatrix{ {10, 80, 1.2 } },
  nCmt(2),
  biovar{ { 1, 1 } },
  tlag{ { 0, 0 } },
  time(10, 0.0),
  amt(10, 0),
  rate(10, 0),
  cmt(10, 2),
  evid(10, 0),
  ii(10, 0),
  addl(10, 0),
  ss(10, 0)
  {
    time[0] = 0.0;
    time[1] = 0.0;
    for(int i = 2; i < 10; i++) time[i] = time[i - 1] + 5;
    cmt[0] = 1;
    evid[0] = 1;
  }
  std::vector<std::vector<double> > pMatrix;  // CL, VC, Ka
  const int nCmt;  // F1, F2
  std::vector<std::vector<double> > biovar;
  std::vector<std::vector<double> > tlag;
  std::vector<double> time;
  std::vector<double> amt;
  std::vector<double> rate;
  std::vector<int> cmt;
  std::vector<int> evid;
  std::vector<double> ii;
  std::vector<int> addl;
  std::vector<int> ss;
  double rel_tol = 1e-8, abs_tol = 1e-8;
  long int max_num_steps = 1e8;

  struct oneCptModelODE_functor {
    template <typename T0, typename T1, typename T2, typename T3>
    inline
    std::vector<typename boost::math::tools::promote_args<T0, T1, T2, T3>::type>
    operator()(const T0& t,
               const std::vector<T1>& x,
               const std::vector<T2>& parms,
               const std::vector<T3>& rate,
               const std::vector<int>& dummy, std::ostream* pstream__) const {
      return oneCptModelODE(t, x, parms, rate, dummy, pstream__);
    }
  };
};

TEST_F(TorstenPKOneCptODETest, steady_state_repeated) {
//  Steady state induced by multiple bolus doses (SS = 1, rate = 0)
  using Eigen::Matrix;
  using Eigen::Dynamic;
  using torsten::generalOdeModel_rk45;
  using torsten::generalOdeModel_bdf; 

  amt[0] = 1200;    
  ss[0] = 1;
  ii[0] = 12;
  addl[0] = 10;

  Matrix<double, Eigen::Dynamic, Eigen::Dynamic> x_rk45{
    generalOdeModel_rk45(oneCptModelODE_functor(), nCmt,
                         time, amt, rate, ii, evid, cmt, addl, ss,
                         pMatrix, biovar, tlag,
                         0,
                         rel_tol, abs_tol, max_num_steps) };
  
  Matrix<double, Eigen::Dynamic, Eigen::Dynamic> x_bdf{
    generalOdeModel_bdf(oneCptModelODE_functor(), nCmt,
                        time, amt, rate, ii, evid, cmt, addl, ss,
                        pMatrix, biovar, tlag,
                        0,
                        rel_tol, abs_tol, max_num_steps) };

  Matrix<double, Dynamic, Dynamic> amounts(10, 2);
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

  for (int i = 0; i < amounts.rows(); i++) {
    EXPECT_NEAR(amounts(i, 0), x_rk45(i, 0), std::max(amounts(i, 0), x_rk45(i, 0)) * 1e-6);
    EXPECT_NEAR(amounts(i, 1), x_rk45(i, 1), std::max(amounts(i, 1), x_rk45(i, 1)) * 1e-6);
    EXPECT_NEAR(amounts(i, 0), x_bdf(i, 0), std::max(amounts(i, 0), x_bdf(i, 0)) * 1e-5);
    EXPECT_NEAR(amounts(i, 1), x_bdf(i, 1), std::max(amounts(i, 1), x_bdf(i, 1)) * 1e-5);
  }

  // Test AutoDiff against FiniteDiff
  double diff = 1e-8, diff2 = 5e-3;
  test_generalOdeModel(oneCptModelODE_functor(), nCmt,
                       time, amt, rate, ii, evid, cmt, addl, ss,
                       pMatrix, biovar, tlag,
                       rel_tol, abs_tol, max_num_steps, diff, diff2, "rk45");
  test_generalOdeModel(oneCptModelODE_functor(), nCmt,
                       time, amt, rate, ii, evid, cmt, addl, ss,
                       pMatrix, biovar, tlag,
                       rel_tol, abs_tol, max_num_steps, diff, diff2, "bdf");
}

TEST_F(TorstenPKOneCptODETest, steady_state_single_infusion) {
  // Steady state with constant rate infusion (SS = 1, rate != 0, ii = 0)
  using std::vector;
  using Eigen::Matrix;
  using Eigen::Dynamic;

  time[0] = 0.0;
  for(int i = 1; i < 10; i++) time[i] = time[i - 1] + 2.5;
  rate[0] = 150;
  ss[0] = 1;

  Matrix<double, Dynamic, Dynamic> x_rk45;
  x_rk45 = torsten::generalOdeModel_rk45(oneCptModelODE_functor(), nCmt,
                                time, amt, rate, ii, evid, cmt, addl, ss,
                                pMatrix, biovar, tlag,
                                0,
                                rel_tol, abs_tol, max_num_steps);

  Matrix<double, Dynamic, Dynamic> x_bdf;
  x_bdf = torsten::generalOdeModel_bdf(oneCptModelODE_functor(), nCmt,
                              time, amt, rate, ii, evid, cmt, addl, ss,
                              pMatrix, biovar, tlag,
                              0,
                              rel_tol, abs_tol, max_num_steps);

  // Couldn't get solution from mrgsolve, so comparing analytical and
  // numerical solutions as a provisional unit test.
  Matrix<double, Dynamic, Dynamic> x_an;
  x_an = torsten::PKModelOneCpt(time, amt, rate, ii, evid, cmt, addl, ss,
                       pMatrix, biovar, tlag);

  double err;
  for (int i = 0; i < x_an.rows(); i++) {
    // relative error is determined empirically
    // (seems to be relatively high -- could be the method is not that precised)
    err = std::max(std::max(x_an(i, 0), x_rk45(i, 0)) * 3.5e-3, 1e-11);
    EXPECT_NEAR(x_an(i, 0), x_rk45(i, 0), err);

    err = std::max(std::max(x_an(i, 1), x_rk45(i, 1)) * 3.5e-3, 1e-11);
    EXPECT_NEAR(x_an(i, 1), x_rk45(i, 1), err);

    err = std::max(std::max(x_an(i, 0), x_bdf(i, 0)) * 1e-2, 1e-11);
    EXPECT_NEAR(x_an(i, 0), x_bdf(i, 0), err);

    err = std::max(std::max(x_an(i, 1), x_bdf(i, 1)) * 1e-2, 1e-11);
    EXPECT_NEAR(x_an(i, 1), x_bdf(i, 1), err);
  }

  // Test AutoDiff against FiniteDiff
  double diff = 1e-8, diff2 = 5e-3;
  test_generalOdeModel(oneCptModelODE_functor(), nCmt,
                       time, amt, rate, ii, evid, cmt, addl, ss,
                       pMatrix, biovar, tlag,
                       rel_tol, abs_tol, max_num_steps, diff, diff2, "rk45");
  test_generalOdeModel(oneCptModelODE_functor(), nCmt,
                       time, amt, rate, ii, evid, cmt, addl, ss,
                       pMatrix, biovar, tlag,
                       rel_tol, abs_tol, max_num_steps, diff, diff2, "bdf");
}

TEST_F(TorstenPKOneCptODETest, steady_state_multi_infusions) {
// TEST(Torsten, genCpt_One_SS_3) {
  // Steady state with multiple truncated infusions
  // (SS = 1, rate != 0, ii > 0)
  using std::vector;
  using Eigen::Matrix;
  using Eigen::Dynamic;

  time[0] = 0.0;
  for(int i = 1; i < 10; i++) time[i] = time[i - 1] + 2.5;
  amt[0] = 1200;
  rate[0] = 150;
  ii[0] = 16;
  ss[0] = 1;

  Matrix<double, Dynamic, Dynamic> x_rk45;
  x_rk45 = torsten::generalOdeModel_rk45(oneCptModelODE_functor(), nCmt,
                                time, amt, rate, ii, evid, cmt, addl, ss,
                                pMatrix, biovar, tlag,
                                0,
                                rel_tol, abs_tol, max_num_steps);

  Matrix<double, Dynamic, Dynamic> x_bdf;
  x_bdf = torsten::generalOdeModel_bdf(oneCptModelODE_functor(), nCmt,
                              time, amt, rate, ii, evid, cmt, addl, ss,
                              pMatrix, biovar, tlag,
                              0,
                              rel_tol, abs_tol, max_num_steps);

  Matrix<double, Dynamic, Dynamic> amounts(10, 2);
  amounts << 8.465519e-03, 360.2470,
             1.187770e+02, 490.4911,
             1.246902e+02, 676.1759,
             1.249846e+02, 816.5263,
             1.133898e+01, 750.0054,
             5.645344e-01, 557.3459,
             2.810651e-02, 408.1926,
             1.399341e-03, 298.6615,
             6.966908e-05, 218.5065,
             3.468619e-06, 159.8628;

  for (int i = 0; i < amounts.rows(); i++) {
    EXPECT_NEAR(amounts(i, 0), x_rk45(i, 0), std::max(amounts(i, 0), x_rk45(i, 0)) * 1e-3);
    EXPECT_NEAR(amounts(i, 1), x_rk45(i, 1), std::max(amounts(i, 1), x_rk45(i, 1)) * 1e-3);
    EXPECT_NEAR(amounts(i, 0), x_bdf(i, 0), std::max(amounts(i, 0), x_bdf(i, 0)) * 1e-3);
    EXPECT_NEAR(amounts(i, 1), x_bdf(i, 1), std::max(amounts(i, 1), x_bdf(i, 1)) * 1e-3);
  }

  // Test AutoDiff against FiniteDiff
  // Currently, torsten does not handle multiple truncated infusions case when
  // amt * F is a parameter (this scenario returns an exception and is not
  // tested here).
  double diff = 1e-8, diff2 = 5e-3;
  test_generalOdeModel(oneCptModelODE_functor(), nCmt,
                       time, amt, rate, ii, evid, cmt, addl, ss,
                       pMatrix, biovar, tlag,
                       rel_tol, abs_tol, max_num_steps, diff, diff2, "rk45",
                       2);

  // diff_bdf2 determined empirically
  double diff_bdf2 = 1e-2;
  test_generalOdeModel(oneCptModelODE_functor(), nCmt,
                       time, amt, rate, ii, evid, cmt, addl, ss,
                       pMatrix, biovar, tlag,
                       rel_tol, abs_tol, max_num_steps, diff, diff_bdf2,
                       "bdf", 2);
}

TEST_F(TorstenPKOneCptODETest, MultipleDose) {
// TEST(Torsten, genCpt_One_MultipleDose) {
  using std::vector;
  using Eigen::Matrix;
  using Eigen::Dynamic;

  double rel_err = 1e-6;

  time[0] = 0.0;
  for(int i = 1; i < 9; i++) time[i] = time[i - 1] + 0.25;
  time[9] = 4.0;

  amt[0] = 1000;
  ii[0] = 12;
  addl[0] = 14;

  Matrix<double, Eigen::Dynamic, Eigen::Dynamic> x_rk45;
  x_rk45 = torsten::generalOdeModel_rk45(oneCptModelODE_functor(), nCmt,
                                time, amt, rate, ii, evid, cmt, addl, ss,
                                pMatrix, biovar, tlag,
                                0,
                                rel_tol, abs_tol, max_num_steps);

  Matrix<double, Eigen::Dynamic, Eigen::Dynamic> x_bdf;
  x_bdf = torsten::generalOdeModel_bdf(oneCptModelODE_functor(), nCmt,
                              time, amt, rate, ii, evid, cmt, addl, ss,
                              pMatrix, biovar, tlag,
                              0,
                              rel_tol, abs_tol, max_num_steps);

  Matrix<double, Dynamic, Dynamic> amounts(10, nCmt);
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

  expect_near_matrix_eq(amounts, x_rk45, rel_err);
  expect_near_matrix_eq(amounts, x_bdf, rel_err);

  // Test AutoDiff against FiniteDiff
  double diff = 1e-8, diff2 = 5e-3;
  test_generalOdeModel(oneCptModelODE_functor(), nCmt,
                       time, amt, rate, ii, evid, cmt, addl, ss,
                       pMatrix, biovar, tlag,
                       rel_tol, abs_tol, max_num_steps, diff, diff2, "rk45");
  test_generalOdeModel(oneCptModelODE_functor(), nCmt,
                       time, amt, rate, ii, evid, cmt, addl, ss,
                       pMatrix, biovar, tlag,
                       rel_tol, abs_tol, max_num_steps, diff, diff2, "bdf");
}

TEST_F(TorstenPKOneCptODETest, MultipleDose_overload) {
  using std::vector;
  using Eigen::Matrix;
  using Eigen::Dynamic;
  using torsten::generalOdeModel_rk45;
  using torsten::generalOdeModel_bdf;
  using F = oneCptModelODE_functor;

  time[0] = 0.0;
  for(int i = 1; i < 9; i++) time[i] = time[i - 1] + 0.25;
  time[9] = 4.0;

  amt[0] = 1000;
  ii[0] = 12;
  addl[0] = 14;

  const double rel_err = 1e-4;
  Matrix<double, Dynamic, Dynamic> amounts(10, 2);
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

  Matrix<double, Eigen::Dynamic, Eigen::Dynamic> x;
  x = generalOdeModel_rk45(F(), nCmt,
                           time, amt, rate, ii, evid, cmt, addl, ss,
                           pMatrix[0], biovar, tlag,
                           0,
                           rel_tol, abs_tol, max_num_steps);
  expect_near_matrix_eq(amounts, x, rel_err);


  x = generalOdeModel_rk45(F(), nCmt,
                                    time, amt, rate, ii, evid, cmt, addl, ss,
                                    pMatrix[0], biovar[0], tlag,
                                    0,
                                    rel_tol, abs_tol, max_num_steps);
  expect_near_matrix_eq(amounts, x, rel_err);

  x = generalOdeModel_rk45(F(), nCmt,
                                    time, amt, rate, ii, evid, cmt, addl, ss,
                                    pMatrix[0], biovar[0], tlag[0],
                                    0,
                                    rel_tol, abs_tol, max_num_steps);
  expect_near_matrix_eq(amounts, x, rel_err);

  x = generalOdeModel_rk45(F(), nCmt,
                                    time, amt, rate, ii, evid, cmt, addl, ss,
                                    pMatrix[0], biovar, tlag[0],
                                    0,
                                    rel_tol, abs_tol, max_num_steps);
  expect_near_matrix_eq(amounts, x, rel_err);

  x = generalOdeModel_rk45(F(), nCmt,
                                    time, amt, rate, ii, evid, cmt, addl, ss,
                                    pMatrix, biovar[0], tlag,
                                    0,
                                    rel_tol, abs_tol, max_num_steps);
  expect_near_matrix_eq(amounts, x, rel_err);

  x = generalOdeModel_rk45(F(), nCmt,
                                    time, amt, rate, ii, evid, cmt, addl, ss,
                                    pMatrix, biovar[0], tlag[0],
                                    0,
                                    rel_tol, abs_tol, max_num_steps);
  expect_near_matrix_eq(amounts, x, rel_err);

  x = generalOdeModel_rk45(F(), nCmt,
                                    time, amt, rate, ii, evid, cmt, addl, ss,
                                    pMatrix, biovar, tlag[0],
                                    0,
                                    rel_tol, abs_tol, max_num_steps);
  expect_near_matrix_eq(amounts, x, rel_err);

  x = generalOdeModel_bdf(F(), nCmt,
                                    time, amt, rate, ii, evid, cmt, addl, ss,
                                    pMatrix[0], biovar, tlag,
                                    0,
                                    rel_tol, abs_tol, max_num_steps);
  expect_near_matrix_eq(amounts, x, rel_err);

  x = generalOdeModel_bdf(F(), nCmt,
                                    time, amt, rate, ii, evid, cmt, addl, ss,
                                    pMatrix[0], biovar[0], tlag,
                                    0,
                                    rel_tol, abs_tol, max_num_steps);
  expect_near_matrix_eq(amounts, x, rel_err);

  x = generalOdeModel_bdf(F(), nCmt,
                                    time, amt, rate, ii, evid, cmt, addl, ss,
                                    pMatrix[0], biovar[0], tlag[0],
                                    0,
                                    rel_tol, abs_tol, max_num_steps);
  expect_near_matrix_eq(amounts, x, rel_err);

  x = generalOdeModel_bdf(F(), nCmt,
                                    time, amt, rate, ii, evid, cmt, addl, ss,
                                    pMatrix[0], biovar, tlag[0],
                                    0,
                                    rel_tol, abs_tol, max_num_steps);
  expect_near_matrix_eq(amounts, x, rel_err);

  x = generalOdeModel_bdf(F(), nCmt,
                                    time, amt, rate, ii, evid, cmt, addl, ss,
                                    pMatrix, biovar[0], tlag,
                                    0,
                                    rel_tol, abs_tol, max_num_steps);
  expect_near_matrix_eq(amounts, x, rel_err);

  x = generalOdeModel_bdf(F(), nCmt,
                                    time, amt, rate, ii, evid, cmt, addl, ss,
                                    pMatrix, biovar[0], tlag[0],
                                    0,
                                    rel_tol, abs_tol, max_num_steps);
  expect_near_matrix_eq(amounts, x, rel_err);

  x = generalOdeModel_bdf(F(), nCmt,
                                    time, amt, rate, ii, evid, cmt, addl, ss,
                                    pMatrix, biovar, tlag[0],
                                    0,
                                    rel_tol, abs_tol, max_num_steps);
  expect_near_matrix_eq(amounts, x, rel_err);
}

TEST_F(TorstenPKOneCptODETest, signature) {
  using stan::math::var;
  using std::vector;
  using Eigen::Matrix;
  using Eigen::Dynamic;

  double rel_err = 1e-4;

  vector<vector<var> > pMatrix_v(1);
  pMatrix_v[0].resize(3);
  pMatrix_v[0][0] = 10;  // CL
  pMatrix_v[0][1] = 80;  // Vc
  pMatrix_v[0][2] = 1.2;  // ka

  vector<vector<var> > biovar_v(1);
  biovar_v[0].resize(nCmt);
  biovar_v[0][0] = 1;  // F1
  biovar_v[0][1] = 1;  // F2

  vector<vector<var> > tlag_v(1);
  tlag_v[0].resize(nCmt);
  tlag_v[0][0] = 0;  // tlag1
  tlag_v[0][1] = 0;  // tlag2

  time[0] = 0.0;
  for(int i = 1; i < 9; i++) time[i] = time[i - 1] + 0.25;
  time[9] = 4.0;

  amt[0] = 1000;
  ii[0] = 12;
  addl[0] = 14;

  Matrix<double, Dynamic, Dynamic> amounts(10, 2);
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

  // RK45
  Matrix<var, Dynamic, Dynamic> res;

#ifndef PK_ONECPT_SIGNATURE_TEST
#define PK_ONECPT_SIGNATURE_TEST(X, PMAT, BIOVAR, TLAG)                 \
  X = torsten::generalOdeModel_rk45(oneCptModelODE_functor(), nCmt,     \
                                    time, amt, rate, ii, evid, cmt, addl, ss, \
                                    PMAT, BIOVAR, TLAG,                 \
                                    0,                                  \
                                    rel_tol, abs_tol, max_num_steps);   \
  for (int j = 0; j < X.size(); j++)                                    \
    EXPECT_NEAR(amounts(j), X(j).val(),                                 \
                std::max(amounts(j), X(j).val()) * rel_err);\

  PK_ONECPT_SIGNATURE_TEST(res, pMatrix_v[0] , biovar      , tlag);
  PK_ONECPT_SIGNATURE_TEST(res, pMatrix_v[0] , biovar_v    , tlag);
  PK_ONECPT_SIGNATURE_TEST(res, pMatrix_v[0] , biovar      , tlag_v);
  PK_ONECPT_SIGNATURE_TEST(res, pMatrix_v[0] , biovar_v    , tlag_v);
  PK_ONECPT_SIGNATURE_TEST(res, pMatrix[0]   , biovar_v    , tlag);
  PK_ONECPT_SIGNATURE_TEST(res, pMatrix[0]   , biovar_v    , tlag_v);
  PK_ONECPT_SIGNATURE_TEST(res, pMatrix[0]   , biovar      , tlag_v);
  PK_ONECPT_SIGNATURE_TEST(res, pMatrix_v[0] , biovar[0]   , tlag);
  PK_ONECPT_SIGNATURE_TEST(res, pMatrix_v[0] , biovar_v[0] , tlag);
  PK_ONECPT_SIGNATURE_TEST(res, pMatrix_v[0] , biovar[0]   , tlag_v);
  PK_ONECPT_SIGNATURE_TEST(res, pMatrix_v[0] , biovar_v[0] , tlag_v);
  PK_ONECPT_SIGNATURE_TEST(res, pMatrix[0]   , biovar_v[0] , tlag);
  PK_ONECPT_SIGNATURE_TEST(res, pMatrix[0]   , biovar_v[0] , tlag_v);
  PK_ONECPT_SIGNATURE_TEST(res, pMatrix[0]   , biovar[0]   , tlag_v);
  PK_ONECPT_SIGNATURE_TEST(res, pMatrix_v[0] , biovar      , tlag[0]);
  PK_ONECPT_SIGNATURE_TEST(res, pMatrix_v[0] , biovar_v    , tlag[0]);
  PK_ONECPT_SIGNATURE_TEST(res, pMatrix_v[0] , biovar      , tlag_v[0]);
  PK_ONECPT_SIGNATURE_TEST(res, pMatrix_v[0] , biovar_v    , tlag_v[0]);
  PK_ONECPT_SIGNATURE_TEST(res, pMatrix[0]   , biovar_v    , tlag[0]);
  PK_ONECPT_SIGNATURE_TEST(res, pMatrix[0]   , biovar_v    , tlag_v[0]);
  PK_ONECPT_SIGNATURE_TEST(res, pMatrix[0]   , biovar      , tlag_v[0]);
  PK_ONECPT_SIGNATURE_TEST(res, pMatrix_v[0] , biovar[0]   , tlag[0]);
  PK_ONECPT_SIGNATURE_TEST(res, pMatrix_v[0] , biovar_v[0] , tlag[0]);
  PK_ONECPT_SIGNATURE_TEST(res, pMatrix_v[0] , biovar[0]   , tlag_v[0]);
  PK_ONECPT_SIGNATURE_TEST(res, pMatrix_v[0] , biovar_v[0] , tlag_v[0]);
  PK_ONECPT_SIGNATURE_TEST(res, pMatrix[0]   , biovar_v[0] , tlag[0]);
  PK_ONECPT_SIGNATURE_TEST(res, pMatrix[0]   , biovar_v[0] , tlag_v[0]);
  PK_ONECPT_SIGNATURE_TEST(res, pMatrix[0]   , biovar[0]   , tlag_v[0]);
  PK_ONECPT_SIGNATURE_TEST(res, pMatrix_v    , biovar[0]   , tlag[0]);
  PK_ONECPT_SIGNATURE_TEST(res, pMatrix_v    , biovar_v[0] , tlag[0]);
  PK_ONECPT_SIGNATURE_TEST(res, pMatrix_v    , biovar[0]   , tlag_v[0]);
  PK_ONECPT_SIGNATURE_TEST(res, pMatrix_v    , biovar_v[0] , tlag_v[0]);
  PK_ONECPT_SIGNATURE_TEST(res, pMatrix      , biovar_v[0] , tlag[0]);
  PK_ONECPT_SIGNATURE_TEST(res, pMatrix      , biovar_v[0] , tlag_v[0]);
  PK_ONECPT_SIGNATURE_TEST(res, pMatrix      , biovar[0]   , tlag_v[0]);
  PK_ONECPT_SIGNATURE_TEST(res, pMatrix_v    , biovar      , tlag[0]);
  PK_ONECPT_SIGNATURE_TEST(res, pMatrix_v    , biovar_v    , tlag[0]);
  PK_ONECPT_SIGNATURE_TEST(res, pMatrix_v    , biovar      , tlag_v[0]);
  PK_ONECPT_SIGNATURE_TEST(res, pMatrix_v    , biovar_v    , tlag_v[0]);
  PK_ONECPT_SIGNATURE_TEST(res, pMatrix      , biovar_v    , tlag[0]);
  PK_ONECPT_SIGNATURE_TEST(res, pMatrix      , biovar_v    , tlag_v[0]);
  PK_ONECPT_SIGNATURE_TEST(res, pMatrix      , biovar      , tlag_v[0]);

#undef PK_ONECPT_SIGNATURE_TEST
#endif

  // BDF
#ifndef PK_ONECPT_SIGNATURE_TEST
#define PK_ONECPT_SIGNATURE_TEST(X, PMAT, BIOVAR, TLAG)                 \
  X = torsten::generalOdeModel_bdf(oneCptModelODE_functor(), nCmt,      \
                                    time, amt, rate, ii, evid, cmt, addl, ss, \
                                    PMAT, BIOVAR, TLAG,                 \
                                    0,                                  \
                                    rel_tol, abs_tol, max_num_steps);   \
  for (int j = 0; j < X.size(); j++)                                    \
    EXPECT_NEAR(amounts(j), X(j).val(),                                 \
                std::max(amounts(j), X(j).val()) * rel_err);

  PK_ONECPT_SIGNATURE_TEST(res, pMatrix_v[0] , biovar      , tlag);
  PK_ONECPT_SIGNATURE_TEST(res, pMatrix_v[0] , biovar_v    , tlag);
  PK_ONECPT_SIGNATURE_TEST(res, pMatrix_v[0] , biovar      , tlag_v);
  PK_ONECPT_SIGNATURE_TEST(res, pMatrix_v[0] , biovar_v    , tlag_v);
  PK_ONECPT_SIGNATURE_TEST(res, pMatrix[0]   , biovar_v    , tlag);
  PK_ONECPT_SIGNATURE_TEST(res, pMatrix[0]   , biovar_v    , tlag_v);
  PK_ONECPT_SIGNATURE_TEST(res, pMatrix[0]   , biovar      , tlag_v);
  PK_ONECPT_SIGNATURE_TEST(res, pMatrix_v[0] , biovar[0]   , tlag);
  PK_ONECPT_SIGNATURE_TEST(res, pMatrix_v[0] , biovar_v[0] , tlag);
  PK_ONECPT_SIGNATURE_TEST(res, pMatrix_v[0] , biovar[0]   , tlag_v);
  PK_ONECPT_SIGNATURE_TEST(res, pMatrix_v[0] , biovar_v[0] , tlag_v);
  PK_ONECPT_SIGNATURE_TEST(res, pMatrix[0]   , biovar_v[0] , tlag);
  PK_ONECPT_SIGNATURE_TEST(res, pMatrix[0]   , biovar_v[0] , tlag_v);
  PK_ONECPT_SIGNATURE_TEST(res, pMatrix[0]   , biovar[0]   , tlag_v);
  PK_ONECPT_SIGNATURE_TEST(res, pMatrix_v[0] , biovar      , tlag[0]);
  PK_ONECPT_SIGNATURE_TEST(res, pMatrix_v[0] , biovar_v    , tlag[0]);
  PK_ONECPT_SIGNATURE_TEST(res, pMatrix_v[0] , biovar      , tlag_v[0]);
  PK_ONECPT_SIGNATURE_TEST(res, pMatrix_v[0] , biovar_v    , tlag_v[0]);
  PK_ONECPT_SIGNATURE_TEST(res, pMatrix[0]   , biovar_v    , tlag[0]);
  PK_ONECPT_SIGNATURE_TEST(res, pMatrix[0]   , biovar_v    , tlag_v[0]);
  PK_ONECPT_SIGNATURE_TEST(res, pMatrix[0]   , biovar      , tlag_v[0]);
  PK_ONECPT_SIGNATURE_TEST(res, pMatrix_v[0] , biovar[0]   , tlag[0]);
  PK_ONECPT_SIGNATURE_TEST(res, pMatrix_v[0] , biovar_v[0] , tlag[0]);
  PK_ONECPT_SIGNATURE_TEST(res, pMatrix_v[0] , biovar[0]   , tlag_v[0]);
  PK_ONECPT_SIGNATURE_TEST(res, pMatrix_v[0] , biovar_v[0] , tlag_v[0]);
  PK_ONECPT_SIGNATURE_TEST(res, pMatrix[0]   , biovar_v[0] , tlag[0]);
  PK_ONECPT_SIGNATURE_TEST(res, pMatrix[0]   , biovar_v[0] , tlag_v[0]);
  PK_ONECPT_SIGNATURE_TEST(res, pMatrix[0]   , biovar[0]   , tlag_v[0]);
  PK_ONECPT_SIGNATURE_TEST(res, pMatrix_v    , biovar[0]   , tlag[0]);
  PK_ONECPT_SIGNATURE_TEST(res, pMatrix_v    , biovar_v[0] , tlag[0]);
  PK_ONECPT_SIGNATURE_TEST(res, pMatrix_v    , biovar[0]   , tlag_v[0]);
  PK_ONECPT_SIGNATURE_TEST(res, pMatrix_v    , biovar_v[0] , tlag_v[0]);
  PK_ONECPT_SIGNATURE_TEST(res, pMatrix      , biovar_v[0] , tlag[0]);
  PK_ONECPT_SIGNATURE_TEST(res, pMatrix      , biovar_v[0] , tlag_v[0]);
  PK_ONECPT_SIGNATURE_TEST(res, pMatrix      , biovar[0]   , tlag_v[0]);
  PK_ONECPT_SIGNATURE_TEST(res, pMatrix_v    , biovar      , tlag[0]);
  PK_ONECPT_SIGNATURE_TEST(res, pMatrix_v    , biovar_v    , tlag[0]);
  PK_ONECPT_SIGNATURE_TEST(res, pMatrix_v    , biovar      , tlag_v[0]);
  PK_ONECPT_SIGNATURE_TEST(res, pMatrix_v    , biovar_v    , tlag_v[0]);
  PK_ONECPT_SIGNATURE_TEST(res, pMatrix      , biovar_v    , tlag[0]);
  PK_ONECPT_SIGNATURE_TEST(res, pMatrix      , biovar_v    , tlag_v[0]);
  PK_ONECPT_SIGNATURE_TEST(res, pMatrix      , biovar      , tlag_v[0]);

#undef PK_ONECPT_SIGNATURE_TEST
#endif

  // CHECK - do I need an AD test for every function signature ?
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

TEST_F(TorstenPKOneCptODETest, abstime_SingleDose) {
  using std::vector;
  using Eigen::Matrix;
  using Eigen::Dynamic;

  double rel_err = 1e-6;

  pMatrix[0].push_back(2);
  pMatrix[0].push_back(1);
  time[0] = 0.0;
  for(int i = 1; i < 9; i++) time[i] = time[i - 1] + 0.25;
  time[9] = 4.0;

  amt[0] = 1000;
  ii[0] = 12;
  addl[0] = 14;

  Matrix<double, Eigen::Dynamic, Eigen::Dynamic> x_rk45;
  x_rk45 = torsten::generalOdeModel_rk45(oneCptModelODE_abstime_functor(), 2,
                                time, amt, rate, ii, evid, cmt, addl, ss,
                                pMatrix, biovar, tlag,
                                0,
                                rel_tol, abs_tol, max_num_steps);

  Matrix<double, Eigen::Dynamic, Eigen::Dynamic> x_bdf;
  x_bdf = torsten::generalOdeModel_bdf(oneCptModelODE_abstime_functor(), 2,
                              time, amt, rate, ii, evid, cmt, addl, ss,
                              pMatrix, biovar, tlag,
                              0,
                              rel_tol, abs_tol, max_num_steps);

  Matrix<double, Dynamic, Dynamic> amounts(10, 2);
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

  expect_near_matrix_eq(amounts, x_rk45, rel_err);
  expect_near_matrix_eq(amounts, x_bdf, rel_err);

  // Test AutoDiff against FiniteDiff
   double diff = 1e-8, diff2 = .25; // CHECK - diff2 seems pretty high!!
   test_generalOdeModel(oneCptModelODE_abstime_functor(), nCmt,
                        time, amt, rate, ii, evid, cmt, addl, ss,
                        pMatrix, biovar, tlag,
                        rel_tol, abs_tol, max_num_steps, diff, diff2, "rk45");
   test_generalOdeModel(oneCptModelODE_abstime_functor(), nCmt,
                        time, amt, rate, ii, evid, cmt, addl, ss,
                        pMatrix, biovar, tlag,
                        rel_tol, abs_tol, max_num_steps, diff, diff2, "bdf");
}

TEST(Torsten, genCptOne_MultipleDoses_timePara) {
  using std::vector;
  using Eigen::Matrix;
  using Eigen::Dynamic;

  double rel_err_rk45 = 1e-6;
  double rel_err_bdf = 1e-4;

  int nEvent = 11;
  vector<vector<double> > pMatrix(nEvent);
  for (int i = 0; i < nEvent; i++) {
    pMatrix[i].resize(5);
    if (i < 6) pMatrix[i][0] = 10; // CL
    else pMatrix[i][0] = 50;
    pMatrix[i][1] = 80; // Vc
    pMatrix[i][2] = 1.2; // ka
  }

  int nCmt = 2;
  vector<vector<double> > biovar(1);
  biovar[0].resize(nCmt);
  biovar[0][0] = 1;  // F1
  biovar[0][1] = 1;  // F2

  vector<vector<double> > tlag(1);
  tlag[0].resize(nCmt);
  tlag[0][0] = 0;  // tlag1
  tlag[0][1] = 0;  // tlag2

	vector<double> time(nEvent);
	time[0] = 0.0;
	for(int i = 1; i < nEvent; i++) time[i] = time[i - 1] + 2.5;

	vector<double> amt(nEvent, 0);
	amt[0] = 1000;

	vector<double> rate(nEvent, 0);

	vector<int> cmt(nEvent, 2);
	cmt[0] = 1;

	vector<int> evid(nEvent, 0);
	evid[0] = 1;

	vector<double> ii(nEvent, 0);
	ii[0] = 12;

	vector<int> addl(nEvent, 0);
	addl[0] = 1;

	vector<int> ss(nEvent, 0);

  double rel_tol = 1e-8, abs_tol = 1e-8;
   long int max_num_steps = 1e8;
   Matrix<double, Eigen::Dynamic, Eigen::Dynamic> x_rk45;
   x_rk45 = torsten::generalOdeModel_rk45(oneCptModelODE_functor(), nCmt,
                                 time, amt, rate, ii, evid, cmt, addl, ss,
                                 pMatrix, biovar, tlag,
                                 0,
                                 rel_tol, abs_tol, max_num_steps);

  Matrix<double, Eigen::Dynamic, Eigen::Dynamic> x_bdf;
  x_bdf = torsten::generalOdeModel_bdf(oneCptModelODE_functor(), nCmt,
                              time, amt, rate, ii, evid, cmt, addl, ss,
                              pMatrix, biovar, tlag,
                              0,
                              rel_tol, abs_tol, max_num_steps);

	Matrix<double, Dynamic, Dynamic> amounts(nEvent, 2);
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

  expect_near_matrix_eq(amounts, x_rk45, rel_err_rk45);
  expect_near_matrix_eq(amounts, x_bdf, rel_err_bdf);

  // Test AutoDiff against FiniteDiff
  double diff = 1e-8, diff2 = 2e-2;
  test_generalOdeModel(oneCptModelODE_functor(), nCmt,
                       time, amt, rate, ii, evid, cmt, addl, ss,
                       pMatrix, biovar, tlag,
                       rel_tol, abs_tol, max_num_steps, diff, diff2, "rk45");
  test_generalOdeModel(oneCptModelODE_functor(), nCmt,
                       time, amt, rate, ii, evid, cmt, addl, ss,
                       pMatrix, biovar, tlag,
                       rel_tol, abs_tol, max_num_steps, diff, diff2, "bdf");
}

TEST_F(TorstenPKOneCptODETest, Rate) {
  using std::vector;
  using Eigen::Matrix;
  using Eigen::Dynamic;

  time[0] = 0;
  for(int i = 1; i < 9; i++) time[i] = time[i - 1] + 0.25;
  time[9] = 4.0;

  amt[0] = 1200;
  rate[0] = 1200;
  ii[0] = 12;
  addl[0] = 14;

  rel_tol = 1e-6;
  abs_tol = 1e-6;
  max_num_steps = 1e6;

  Matrix<double, Dynamic, Dynamic>
    x_rk45 = torsten::generalOdeModel_rk45(oneCptModelODE_functor(), nCmt,
                                  time, amt, rate, ii, evid, cmt, addl, ss,
                                  pMatrix, biovar, tlag,
                                  0,
                                  rel_tol, abs_tol, max_num_steps);

  Matrix<double, Dynamic, Dynamic>
    x_bdf = torsten::generalOdeModel_bdf(oneCptModelODE_functor(), nCmt,
                                time, amt, rate, ii, evid, cmt, addl, ss,
                                pMatrix, biovar, tlag,
                                0,
                                rel_tol, abs_tol, max_num_steps);

  Matrix<double, Dynamic, Dynamic> amounts(10, 2);
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

  // rel err determined empirically
  double rel_err_rk45 = 1e-6, rel_err_bdf = 1e-5;
  expect_near_matrix_eq(amounts, x_rk45, rel_err_rk45);
  expect_near_matrix_eq(amounts, x_bdf, rel_err_bdf);

  // Test Autodiff
  double diff = 1e-8, diff2 = 2e-2;
  test_generalOdeModel(oneCptModelODE_functor(), nCmt,
                       time, amt, rate, ii, evid, cmt, addl, ss,
                       pMatrix, biovar, tlag,
                       rel_tol, abs_tol, max_num_steps, diff, diff2, "rk45");
  test_generalOdeModel(oneCptModelODE_functor(), nCmt,
                       time, amt, rate, ii, evid, cmt, addl, ss,
                       pMatrix, biovar, tlag,
                       rel_tol, abs_tol, max_num_steps, diff, diff2, "bdf");
}


struct twoCptModelODE_functor {
  template <typename T0, typename T1, typename T2, typename T3>
  inline
  std::vector<typename boost::math::tools::promote_args<T0, T1, T2, T3>::type>
  operator()(const T0& t,
             const std::vector<T1>& x,
             const std::vector<T2>& parms,
             const std::vector<T3>& rate,
             const std::vector<int>& dummy, std::ostream* pstream__) const {
    typedef typename boost::math::tools::promote_args<T0, T1, T2, T3>::type scalar;

    scalar
      CL = parms[0],
      Q = parms[1],
      V1 = parms[2],
      V2 = parms[3],
      ka = parms[4],
      k10 = CL / V1,
      k12 = Q / V1,
      k21 = Q / V2;

    std::vector<scalar> y(3, 0);
    y[0] = -ka * x[0];
    y[1] = ka * x[0] - (k10 + k12) * x[1] + k21 * x[2];
    y[2] = k12 * x[1] - k21 * x[2];

    return y;
  }
};


TEST(Torsten, generalTwoCptModel_Rate) {
  using std::vector;
  using Eigen::Matrix;
  using Eigen::Dynamic;

  vector<vector<double> > pMatrix(1);
  pMatrix[0].resize(5);
  pMatrix[0][0] = 5;  // CL
  pMatrix[0][1] = 8;  // Q
  pMatrix[0][2] = 35;  // Vc
  pMatrix[0][3] = 105;  // Vp
  pMatrix[0][4] = 1.2;  // ka

  vector<vector<double> > biovar(1);
  biovar[0].resize(3);
  biovar[0][0] = 1;  // F1
  biovar[0][1] = 1;  // F2
  biovar[0][2] = 1;  // F3

  vector<vector<double> > tlag(1);
  tlag[0].resize(3);
  tlag[0][0] = 0;  // tlag1
  tlag[0][1] = 0;  // tlag2
  tlag[0][2] = 0;  // tlag3

  vector<double> time(10);
  time[0] = 0;
  for(int i = 1; i < 9; i++) time[i] = time[i - 1] + 0.25;
  time[9] = 4.0;

  vector<double> amt(10, 0);
  amt[0] = 1200;

  vector<double> rate(10, 0);
  rate[0] = 1200;  // non-zero rate causes error with jacobians

  vector<int> cmt(10, 2);
  cmt[0] = 1;

  vector<int> evid(10, 0);
  evid[0] = 1;

  vector<double> ii(10, 0);
  ii[0] = 12;

  vector<int> addl(10, 0);
  addl[0] = 14;

  vector<int> ss(10, 0);

  double rel_tol = 1e-6, abs_tol = 1e-6;
  long int max_num_steps = 1e6;

  int nCmt = 3;

  Matrix<double, Dynamic, Dynamic> x_rk45, x_bdf;
  x_rk45 = torsten::generalOdeModel_rk45(twoCptModelODE_functor(), nCmt,
                                time, amt, rate, ii, evid, cmt, addl, ss,
                                pMatrix, biovar, tlag,
                                0,
                                rel_tol, abs_tol, max_num_steps);
  x_bdf = torsten::generalOdeModel_bdf(twoCptModelODE_functor(), nCmt,
                              time, amt, rate, ii, evid, cmt, addl, ss,
                              pMatrix, biovar, tlag,
                              0,
                              rel_tol, abs_tol, max_num_steps);

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

  // relative error determined empirically
  double rel_err_rk45 = 1e-6, rel_err_bdf = 1e-4;
  expect_near_matrix_eq(amounts, x_rk45, rel_err_rk45);
  expect_near_matrix_eq(amounts, x_bdf, rel_err_bdf);

  // Test Autodiff
  double diff = 1e-8, diff2 = 2e-2;
  test_generalOdeModel(twoCptModelODE_functor(), nCmt,
                       time, amt, rate, ii, evid, cmt, addl, ss,
                       pMatrix, biovar, tlag,
                       rel_tol, abs_tol, max_num_steps, diff, diff2, "rk45");
  test_generalOdeModel(twoCptModelODE_functor(), nCmt,
                       time, amt, rate, ii, evid, cmt, addl, ss,
                       pMatrix, biovar, tlag,
                       rel_tol, abs_tol, max_num_steps, diff, diff2, "bdf");
  }

struct FK_functor {
   // parms contains both the PK and the PD parameters.
   // x contains both the PK and the PD states.
  template <typename T0, typename T1, typename T2, typename T3>
  inline
  std::vector<typename boost::math::tools::promote_args<T0, T1, T2, T3>::type>
  operator()(const T0& t,
             const std::vector<T1>& x,
             const std::vector<T2>& parms,
             const std::vector<T3>& x_r,
             const std::vector<int>& x_i,
             std::ostream* pstream__) const {
    using Eigen::Matrix;
    using Eigen::Dynamic;
    typedef typename boost::math::tools::promote_args<T0, T1, T2, T3>::type scalar;

    // PK variables
    scalar
      CL = parms[0],
      Q = parms[1],
      VC = parms[2],
      VP = parms[3],
      ka = parms[4],
      k10 = CL / VC,
      k12 = Q / VC,
      k21 = Q / VP;

    // PD variables
    scalar
      MTT = parms[5],
      circ0 = parms[6],
      alpha = parms[7],
      gamma = parms[8],
      ktr = 4 / MTT,
      prol = x[3] + circ0,
      transit1 = x[4] + circ0,
      transit2 = x[5] + circ0,
      transit3 = x[6] + circ0,
      circ = x[7] + circ0;

    std::vector<scalar> dxdt(8);
    dxdt[0] = -ka * x[0];
    dxdt[1] = ka * x[0] - (k10 + k12) * x[1] + k21 * x[2];
    dxdt[2] = k12 * x[1] - k21 * x[2];

    scalar conc = x[1] / VC;
    scalar Edrug = alpha * conc;

    dxdt[3] = ktr * prol * ((1 - Edrug) * pow((circ0 / circ), gamma) - 1);
    dxdt[4] = ktr * (prol - transit1);
    dxdt[5] = ktr * (transit1 - transit2);
    dxdt[6] = ktr * (transit2 - transit3);
    dxdt[7] = ktr * (transit3 - circ);

    return dxdt;
  }
};

TEST(Torsten, genCpt_FK_SS) {
  using std::vector;
  using Eigen::Matrix;
  using Eigen::Dynamic;

  vector<double> theta(9);
  theta[0] = 10;  // CL
  theta[1] = 15;  // Q
  theta[2] = 35;  // Vc
  theta[3] = 105;  // Vp
  theta[4] = 2.0;  // ka
  theta[5] = 125;  // MTT
  theta[6] = 5;  // Circ0
  theta[7] = 3e-4;  // alpha
  theta[8] = 0.17;  // gamma

  vector<vector<double> > theta_v(1, theta);

  int nCmt = 8;
  vector<double> biovar(nCmt);
  for (int i = 0; i < nCmt; i++)
    biovar[i] = 1;

  vector<vector<double> > biovar_v(1, biovar);

  vector<double> tlag(nCmt);
  for (int i = 0; i < nCmt; i++)
    tlag[i] = 0;

  vector<vector<double> > tlag_v(1, tlag);

  vector<double> time(10);
  time[0] = 0;
  for(int i = 1; i < 10; i++) time[i] = time[i - 1] + 1.25;

  vector<double> amt(10, 0);
  amt[0] = 80 * 1000;

  vector<double> rate(10, 0);

  vector<int> cmt(10, 2);
  cmt[0] = 1;

  vector<int> evid(10, 0);
  evid[0] = 1;

  vector<double> ii(10, 0);
  ii[0] = 12;

  vector<int> addl(10, 0);

  vector<int> ss(10, 0);
  ss[0] = 1;

  double rel_tol = 1e-6, abs_tol = 1e-6;
  long int max_num_steps = 1e6;

  Matrix<double, Dynamic, Dynamic> x_rk45, x_bdf;
  x_rk45 = torsten::generalOdeModel_rk45(FK_functor(), nCmt,
                                time, amt, rate, ii, evid, cmt, addl, ss,
                                theta, biovar, tlag,
                                0,
                                rel_tol, abs_tol, max_num_steps);

  x_bdf = torsten::generalOdeModel_bdf(FK_functor(), nCmt,
                              time, amt, rate, ii, evid, cmt, addl, ss,
                              theta, biovar, tlag,
                              0,
                              rel_tol, abs_tol, max_num_steps);

  Matrix<double, Dynamic, Dynamic> amounts(10, 8);
  amounts << 8.000000e+04, 11996.63, 55694.35, -3.636308, -3.653620,  -3.653933, -3.653748, -3.653622,
    6.566800e+03, 53123.67, 70649.28, -3.650990, -3.653172, -3.653910, -3.653755, -3.653627,
    5.390358e+02, 34202.00, 80161.15, -3.662446, -3.653349, -3.653883, -3.653761, -3.653632,
    4.424675e+01, 23849.69, 80884.40, -3.665321, -3.653782, -3.653870, -3.653765, -3.653637,
    3.631995e+00, 19166.83, 78031.24, -3.664114, -3.654219, -3.653876, -3.653769, -3.653642,
    2.981323e-01, 16799.55, 74020.00, -3.660988, -3.654550, -3.653896, -3.653774, -3.653647,
    2.447219e-02, 15333.26, 69764.65, -3.656791, -3.654722, -3.653926, -3.653779, -3.653653,
    2.008801e-03, 14233.96, 65591.05, -3.651854, -3.654708, -3.653957, -3.653786, -3.653658,
    1.648918e-04, 13303.26, 61607.92, -3.646317, -3.654488, -3.653983, -3.653793, -3.653663,
    1.353552e-05, 12466.56, 57845.10, -3.640244, -3.654050, -3.653995, -3.653801, -3.653668;

  // relative error determined empirically (12%)
  double rel_err_rk45 = 1.2e-2, rel_err_bdf = 1.2e-2;
  expect_near_matrix_eq(amounts, x_rk45, rel_err_rk45);
  expect_near_matrix_eq(amounts, x_bdf, rel_err_bdf);

  // Test Autodiff
  double diff = 1e-8, diff2 = 2e-2;
  test_generalOdeModel(FK_functor(), nCmt,
                       time, amt, rate, ii, evid, cmt, addl, ss,
                       theta_v, biovar_v, tlag_v,
                       rel_tol, abs_tol, max_num_steps, diff, diff2, "rk45");

  std::cout << "WARNING: GRADIENT TESTS FOR GENERAL_ODE_BDF FAILS."
            << " SEE ISSUE 45."
            << std::endl;
  
  // gradients do not get properly evaluated in the bdf case!!!  
  // test_generalOdeModel(FK_functor(), nCmt,
  //                      time, amt, rate, ii, evid, cmt, addl, ss,
  //                      theta_v, biovar_v, tlag_v,
  //                      rel_tol, abs_tol, max_num_steps, diff, diff2, "bdf");
}
