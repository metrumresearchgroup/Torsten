#include <gtest/gtest.h>
#include <stan/math/rev/mat.hpp>  // FIX ME - includes should be more specific
#include <gtest/gtest.h>
#include <test/unit/math/torsten/expect_near_matrix_eq.hpp>
#include <test/unit/math/torsten/expect_matrix_eq.hpp>
#include <stan/math/torsten/linOdeModel.hpp>
#include <stan/math/torsten/pk_linode_model.hpp>
#include <test/unit/math/torsten/util_linOdeModel.hpp>
#include <vector>

using std::vector;
using Eigen::Matrix;
using Eigen::Dynamic;

class TorstenLinODEOneCptTest : public testing::Test {
  void SetUp() {
    // make sure memory's clean before starting each test
    stan::math::recover_memory();
  }
public:
  TorstenLinODEOneCptTest() :
    CL(10),
    Vc(80),
    ka(1.2),
    k10(CL / Vc),
    system(2, 2),
    system_array(1),
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
    system << -ka, 0, ka, -k10;
    system_array[0] = system;

    time[0] = 0.0;
    time[1] = 0.0;
    for(int i = 2; i < 10; i++) time[i] = time[i - 1] + 5;

    amt[0] = 1200;
    cmt[0] = 1;
    evid[0] = 1;
    ii[0] = 12;
    addl[0] = 10;
    SetUp();
  }

  double CL, Vc, ka, k10;
  Matrix<double, Dynamic, Dynamic> system;
  vector<Matrix<double, Dynamic, Dynamic> > system_array;
  vector<vector<double> > biovar;
  vector<vector<double> > tlag;
  vector<double> time;
  vector<double> amt;
  vector<double> rate;
  vector<int> cmt;
  vector<int> evid;
  vector<double> ii;
  vector<int> addl;
  vector<int> ss;
};

TEST_F(TorstenLinODEOneCptTest, steady_state) {
  ss[0] = 1;

  Matrix<double, Dynamic, Dynamic> x;
  x = torsten::linOdeModel(time, amt, rate, ii, evid, cmt, addl, ss,
                           system_array, biovar, tlag);

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

  for(int i = 0; i < amounts.rows(); i++) {
    EXPECT_NEAR(amounts(i, 0), x(i, 0), std::max(amounts(i, 0), x(i, 0)) * 1e-6);
    EXPECT_NEAR(amounts(i, 1), x(i, 1), std::max(amounts(i, 1), x(i, 1)) * 1e-6);
  }

  // Test auto-diff
  double diff = 1e-8, diff2 = 1e-4;
  test_linOdeModel(time, amt, rate, ii, evid, cmt, addl, ss,
                   system_array, biovar, tlag, diff, diff2);
}

TEST_F(TorstenLinODEOneCptTest, steady_state_overloads) {
  ss[0] = 1;

  Matrix<double, Dynamic, Dynamic> x_122, x_112, x_111, x_121, x_212,
    x_211, x_221;

  x_122 = torsten::linOdeModel(time, amt, rate, ii, evid, cmt, addl, ss,
                               system_array[0], biovar, tlag);
  x_112 = torsten::linOdeModel(time, amt, rate, ii, evid, cmt, addl, ss,
                               system_array[0], biovar[0], tlag);
  x_111 = torsten::linOdeModel(time, amt, rate, ii, evid, cmt, addl, ss,
                               system_array[0], biovar[0], tlag[0]);
  x_121 = torsten::linOdeModel(time, amt, rate, ii, evid, cmt, addl, ss,
                               system_array[0], biovar, tlag[0]);
  x_212 = torsten::linOdeModel(time, amt, rate, ii, evid, cmt, addl, ss,
                               system_array, biovar[0], tlag);
  x_211 = torsten::linOdeModel(time, amt, rate, ii, evid, cmt, addl, ss,
                               system_array, biovar[0], tlag[0]);
  x_221 = torsten::linOdeModel(time, amt, rate, ii, evid, cmt, addl, ss,
                               system_array, biovar, tlag[0]);

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

  for(int i = 0; i < amounts.rows(); i++) {
    EXPECT_NEAR(amounts(i, 0), x_122(i, 0),
                std::max(amounts(i, 0), x_122(i, 0)) * 1e-6);
    EXPECT_NEAR(amounts(i, 1), x_122(i, 1),
                std::max(amounts(i, 1), x_122(i, 1)) * 1e-6);
		
    EXPECT_NEAR(amounts(i, 0), x_112(i, 0),
                std::max(amounts(i, 0), x_112(i, 0)) * 1e-6);
    EXPECT_NEAR(amounts(i, 1), x_112(i, 1),
                std::max(amounts(i, 1), x_112(i, 1)) * 1e-6);
		
    EXPECT_NEAR(amounts(i, 0), x_111(i, 0),
                std::max(amounts(i, 0), x_111(i, 0)) * 1e-6);
    EXPECT_NEAR(amounts(i, 1), x_111(i, 1),
                std::max(amounts(i, 1), x_111(i, 1)) * 1e-6);
		
    EXPECT_NEAR(amounts(i, 0), x_121(i, 0),
                std::max(amounts(i, 0), x_121(i, 0)) * 1e-6);
    EXPECT_NEAR(amounts(i, 1), x_121(i, 1),
                std::max(amounts(i, 1), x_121(i, 1)) * 1e-6);
		
    EXPECT_NEAR(amounts(i, 0), x_212(i, 0),
                std::max(amounts(i, 0), x_212(i, 0)) * 1e-6);
    EXPECT_NEAR(amounts(i, 1), x_212(i, 1),
                std::max(amounts(i, 1), x_212(i, 1)) * 1e-6);
		
    EXPECT_NEAR(amounts(i, 0), x_211(i, 0),
                std::max(amounts(i, 0), x_211(i, 0)) * 1e-6);
    EXPECT_NEAR(amounts(i, 1), x_211(i, 1),
                std::max(amounts(i, 1), x_211(i, 1)) * 1e-6);
		
    EXPECT_NEAR(amounts(i, 0), x_221(i, 0),
                std::max(amounts(i, 0), x_221(i, 0)) * 1e-6);
    EXPECT_NEAR(amounts(i, 1), x_221(i, 1),
                std::max(amounts(i, 1), x_221(i, 1)) * 1e-6);
  }
}

TEST_F(TorstenLinODEOneCptTest, signature) {
  using stan::math::var;

  var CL_v = 10, Vc_v = 80, ka_v = 1.2, k10_v = CL_v / Vc_v;
  Matrix<var, Dynamic, Dynamic> system_v(2, 2);
  system_v << -ka_v, 0, ka_v, -k10_v;
  vector<Matrix<var, Dynamic, Dynamic> > system_array_v(1, system_v);
  
  vector<vector<var> > biovar_v(1);
  biovar_v[0].resize(2);
  biovar_v[0][0] = 1;  // F1
  biovar_v[0][1] = 1;  // F2

  vector<vector<var> > tlag_v(1);
  tlag_v[0].resize(2);
  tlag_v[0][0] = 0;  // tlag 1
  tlag_v[0][1] = 0;  // tlag 2

  ss[0] = 1;
  
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
  
  Matrix<var, Dynamic, Dynamic> res;
  const double rel_err = 1e-6;
#ifndef PK_ONECPT_ODE_SIGNATURE_TEST
#define PK_ONECPT_ODE_SIGNATURE_TEST(X, PMAT, BIOVAR, TLAG)             \
  X = torsten::linOdeModel(time, amt, rate, ii, evid, cmt, addl, ss,    \
                           PMAT, BIOVAR, TLAG);                         \
  for (int j = 0; j < X.size(); j++)                                    \
    EXPECT_NEAR(amounts(j), X(j).val(),                                 \
                std::max(amounts(j), X(j).val()) * rel_err);

  PK_ONECPT_ODE_SIGNATURE_TEST(res, system_array_v[0] , biovar      , tlag);
  PK_ONECPT_ODE_SIGNATURE_TEST(res, system_array_v[0] , biovar      , tlag);
  PK_ONECPT_ODE_SIGNATURE_TEST(res, system_array_v[0] , biovar_v    , tlag);
  PK_ONECPT_ODE_SIGNATURE_TEST(res, system_array_v[0] , biovar      , tlag_v);
  PK_ONECPT_ODE_SIGNATURE_TEST(res, system_array_v[0] , biovar_v    , tlag_v);
  PK_ONECPT_ODE_SIGNATURE_TEST(res, system_array[0]   , biovar_v    , tlag);
  PK_ONECPT_ODE_SIGNATURE_TEST(res, system_array[0]   , biovar_v    , tlag_v);
  PK_ONECPT_ODE_SIGNATURE_TEST(res, system_array[0]   , biovar      , tlag_v);
  PK_ONECPT_ODE_SIGNATURE_TEST(res, system_array_v[0] , biovar[0]   , tlag);
  PK_ONECPT_ODE_SIGNATURE_TEST(res, system_array_v[0] , biovar_v[0] , tlag);
  PK_ONECPT_ODE_SIGNATURE_TEST(res, system_array_v[0] , biovar[0]   , tlag_v);
  PK_ONECPT_ODE_SIGNATURE_TEST(res, system_array_v[0] , biovar_v[0] , tlag_v);
  PK_ONECPT_ODE_SIGNATURE_TEST(res, system_array[0]   , biovar_v[0] , tlag);
  PK_ONECPT_ODE_SIGNATURE_TEST(res, system_array[0]   , biovar_v[0] , tlag_v);
  PK_ONECPT_ODE_SIGNATURE_TEST(res, system_array[0]   , biovar[0]   , tlag_v);
  PK_ONECPT_ODE_SIGNATURE_TEST(res, system_array_v[0] , biovar      , tlag[0]);
  PK_ONECPT_ODE_SIGNATURE_TEST(res, system_array_v[0] , biovar_v    , tlag[0]);
  PK_ONECPT_ODE_SIGNATURE_TEST(res, system_array_v[0] , biovar      , tlag_v[0]);
  PK_ONECPT_ODE_SIGNATURE_TEST(res, system_array_v[0] , biovar_v    , tlag_v[0]);
  PK_ONECPT_ODE_SIGNATURE_TEST(res, system_array[0]   , biovar_v    , tlag[0]);
  PK_ONECPT_ODE_SIGNATURE_TEST(res, system_array[0]   , biovar_v    , tlag_v[0]);
  PK_ONECPT_ODE_SIGNATURE_TEST(res, system_array[0]   , biovar      , tlag_v[0]);
  PK_ONECPT_ODE_SIGNATURE_TEST(res, system_array_v[0] , biovar[0]   , tlag[0]);
  PK_ONECPT_ODE_SIGNATURE_TEST(res, system_array_v[0] , biovar_v[0] , tlag[0]);
  PK_ONECPT_ODE_SIGNATURE_TEST(res, system_array_v[0] , biovar[0]   , tlag_v[0]);
  PK_ONECPT_ODE_SIGNATURE_TEST(res, system_array_v[0] , biovar_v[0] , tlag_v[0]);
  PK_ONECPT_ODE_SIGNATURE_TEST(res, system_array[0]   , biovar_v[0] , tlag[0]);
  PK_ONECPT_ODE_SIGNATURE_TEST(res, system_array[0]   , biovar_v[0] , tlag_v[0]);
  PK_ONECPT_ODE_SIGNATURE_TEST(res, system_array[0]   , biovar[0]   , tlag_v[0]);
  PK_ONECPT_ODE_SIGNATURE_TEST(res, system_array_v    , biovar[0]   , tlag);
  PK_ONECPT_ODE_SIGNATURE_TEST(res, system_array_v    , biovar_v[0] , tlag);
  PK_ONECPT_ODE_SIGNATURE_TEST(res, system_array_v    , biovar[0]   , tlag_v);
  PK_ONECPT_ODE_SIGNATURE_TEST(res, system_array_v    , biovar_v[0] , tlag_v);
  PK_ONECPT_ODE_SIGNATURE_TEST(res, system_array      , biovar_v[0] , tlag);
  PK_ONECPT_ODE_SIGNATURE_TEST(res, system_array      , biovar_v[0] , tlag_v);
  PK_ONECPT_ODE_SIGNATURE_TEST(res, system_array      , biovar[0]   , tlag_v);
  PK_ONECPT_ODE_SIGNATURE_TEST(res, system_array_v    , biovar[0]   , tlag[0]);
  PK_ONECPT_ODE_SIGNATURE_TEST(res, system_array_v    , biovar_v[0] , tlag[0]);
  PK_ONECPT_ODE_SIGNATURE_TEST(res, system_array_v    , biovar[0]   , tlag_v[0]);
  PK_ONECPT_ODE_SIGNATURE_TEST(res, system_array_v    , biovar_v[0] , tlag_v[0]);
  PK_ONECPT_ODE_SIGNATURE_TEST(res, system_array      , biovar_v[0] , tlag[0]);
  PK_ONECPT_ODE_SIGNATURE_TEST(res, system_array      , biovar_v[0] , tlag_v[0]);
  PK_ONECPT_ODE_SIGNATURE_TEST(res, system_array      , biovar[0]   , tlag_v[0]);
  PK_ONECPT_ODE_SIGNATURE_TEST(res, system_array_v    , biovar      , tlag[0]);
  PK_ONECPT_ODE_SIGNATURE_TEST(res, system_array_v    , biovar_v    , tlag[0]);
  PK_ONECPT_ODE_SIGNATURE_TEST(res, system_array_v    , biovar      , tlag_v[0]);
  PK_ONECPT_ODE_SIGNATURE_TEST(res, system_array_v    , biovar_v    , tlag_v[0]);
  PK_ONECPT_ODE_SIGNATURE_TEST(res, system_array      , biovar_v    , tlag[0]);
  PK_ONECPT_ODE_SIGNATURE_TEST(res, system_array      , biovar_v    , tlag_v[0]);
  PK_ONECPT_ODE_SIGNATURE_TEST(res, system_array      , biovar      , tlag_v[0]);
  
  // CHECK - do I need an AD test for every function signature ?
#endif

}

TEST_F(TorstenLinODEOneCptTest, steady_state_rate) {
  rate[0] = 150;
  ss[0] = 1;

  Matrix<double, Dynamic, Dynamic> x;
  x = torsten::linOdeModel(time, amt, rate, ii, evid, cmt, addl, ss,
                           system_array, biovar, tlag);

  Matrix<double, Dynamic, Dynamic> amounts(10, 2);
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

  for(int i = 0; i < amounts.rows(); i++) {
    EXPECT_NEAR(amounts(i, 0), x(i, 0), std::max(amounts(i, 0), x(i, 0)) * 1e-6);
    EXPECT_NEAR(amounts(i, 1), x(i, 1), std::max(amounts(i, 1), x(i, 1)) * 1e-6);
  }

  double diff = 1e-8, diff2 = 1e-4;
  test_linOdeModel(time, amt, rate, ii, evid, cmt, addl, ss,
                   system_array, biovar, tlag, diff, diff2);

}

TEST(Torsten, linOne_MultipleDoses_timePara) {

    int nEvent = 11;
    
    vector<vector<double> > biovar(1);
    biovar[0].resize(2);
    biovar[0][0] = 1;  // F1
    biovar[0][1] = 1;  // F2
    
    vector<vector<double> > tlag(1);
    tlag[0].resize(2);
    tlag[0][0] = 0;  // tlag1
    tlag[0][1] = 0;  // tlag2
    
  double CL_1 = 10, CL_2 = 50, Vc = 80, ka = 1.2,
      k10_1 = CL_1 / Vc, k10_2 = CL_2 / Vc;
	Matrix<double, Dynamic, Dynamic> system_1(2, 2), system_2(2, 2);
	system_1 << -ka, 0, ka, -k10_1;
	system_2 << -ka, 0, ka, -k10_2;
	
	vector<Matrix<double, Dynamic, Dynamic> > system_array(nEvent, system_1);
	for (int i = 6; i < nEvent; i++)
	  system_array[i] = system_2;

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

	Matrix<double, Dynamic, Dynamic> x;
	x = torsten::linOdeModel(time, amt, rate, ii, evid, cmt, addl, ss,
                 system_array, biovar, tlag);

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

	expect_matrix_eq(amounts, x);

	double diff = 1e-8, diff2 = 1e-4;
	test_linOdeModel(time, amt, rate, ii, evid, cmt, addl, ss,
		              system_array, biovar, tlag, diff, diff2);
}

TEST_F(TorstenLinODEOneCptTest, rate) {
  using std::vector;

  time[0] = 0;
  for(int i = 1; i < 9; i++) time[i] = time[i - 1] + 0.25;
  time[9] = 4.0;

  rate[0] = 1200;

  addl[0] = 10;

  Matrix<double, Dynamic, Dynamic>
    x = torsten::linOdeModel(time, amt, rate, ii, evid, cmt, addl, ss,
                    system_array, biovar, tlag);

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

  expect_matrix_eq(amounts, x);
}

TEST_F(TorstenLinODEOneCptTest, amt_var) {
  using std::vector;
  using stan::math::var;

  time[0] = 0;
  for(int i = 1; i < 9; i++) time[i] = time[i - 1] + 0.25;
  time[9] = 4.0;

  rate[0] = 0;
  amt[1] = 1000;
  std::vector<stan::math::var> amt_v {stan::math::to_var(amt)};

  std::vector<double> g;
  Matrix<var, Dynamic, Dynamic> x;
  Matrix<double, Dynamic, Dynamic> x1, x2;
  const double h = 0.01;
  auto test_it = [&](Matrix<var, Dynamic, Dynamic>& res) {
    for (int i = 0; i < res.rows(); ++i) {
      res(i, 1).grad(amt_v, g);
      EXPECT_NEAR(g[0], (x1(i, 1) - x2(i, 1))/(2 * h), 1e-6);
      stan::math::set_zero_all_adjoints();

      // for gut there is a non-smooth point, skip it
      res(i, 0).grad(amt_v, g);
      EXPECT_NEAR(g[0], (x1(i, 0) - x2(i, 0))/(2 * h), 1e-6);
      stan::math::set_zero_all_adjoints();
    }
  };

  vector<Matrix<var, Dynamic, Dynamic> > system_array_v(1);
  system_array_v[0] = stan::math::to_var(system_array[0]);
  vector<vector<var> > biovar_v{ { 1, 1 } };
  vector<vector<var> > tlag_v{ { 0, 0 } };


  amt[0] += h;
  x1 = torsten::linOdeModel(time, amt, rate, ii, evid, cmt, addl, ss,
                            system_array, biovar, tlag);
  amt[0] -= 2 * h;
  x2 = torsten::linOdeModel(time, amt, rate, ii, evid, cmt, addl, ss,
                            system_array, biovar, tlag);

  x = torsten::linOdeModel(time, amt_v, rate, ii, evid, cmt, addl, ss,
                           system_array, biovar, tlag);

#ifndef LINODE_ONECPT_AMT_VAR_TEST
#define LINODE_ONECPT_AMT_VAR_TEST(SYS, BIOVAR, TLAG)                   \
  x = torsten::linOdeModel(time, amt_v, rate, ii, evid, cmt, addl, ss,  \
                           SYS, BIOVAR, TLAG);                          \
  test_it(x);

  LINODE_ONECPT_AMT_VAR_TEST(system_array   , biovar   , tlag);
  LINODE_ONECPT_AMT_VAR_TEST(system_array_v , biovar   , tlag);
  LINODE_ONECPT_AMT_VAR_TEST(system_array_v , biovar_v , tlag);
  LINODE_ONECPT_AMT_VAR_TEST(system_array_v , biovar_v , tlag_v);
  LINODE_ONECPT_AMT_VAR_TEST(system_array   , biovar_v , tlag);
  LINODE_ONECPT_AMT_VAR_TEST(system_array   , biovar_v , tlag_v);
  LINODE_ONECPT_AMT_VAR_TEST(system_array   , biovar   , tlag_v);
  LINODE_ONECPT_AMT_VAR_TEST(system_array_v , biovar_v , tlag_v);

#undef LINODE_ONECPT_AMT_VAR_TEST
#endif
}

TEST_F(TorstenLinODEOneCptTest, rate_var) {
  using std::vector;
  using stan::math::var;

  time[0] = 0;
  for(int i = 1; i < 9; i++) time[i] = time[i - 1] + 0.25;
  time[9] = 4.0;

  rate[0] = 1200;
  const double t_cutoff = amt[0]/rate[0];
  std::vector<stan::math::var> rate_v {stan::math::to_var(rate)};

  const double h = 0.001;
  Matrix<var, Dynamic, Dynamic> x;
  Matrix<double, Dynamic, Dynamic> x1, x2;

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

  std::vector<double> g;
  const double tol = 1e-6;
  auto test_it = [&](Matrix<var, Dynamic, Dynamic>& res) {
    // test val
    expect_matrix_eq(amounts, stan::math::value_of(res));

    // test adj
    for (int i = 0; i < res.rows(); ++i) {
      res(i, 1).grad(rate_v, g);
      stan::math::set_zero_all_adjoints();
      EXPECT_NEAR(g[0], (x1(i, 1) - x2(i, 1))/(2 * h), tol);

      // skip the non-smooth location
      res(i, 0).grad(rate_v, g);
      stan::math::set_zero_all_adjoints();
      if (time[i] != t_cutoff) {
        EXPECT_NEAR(g[0], (x1(i, 0) - x2(i, 0))/(2 * h), tol);
      }
    }
  };

  vector<Matrix<var, Dynamic, Dynamic> > system_array_v(1);
  system_array_v[0] = stan::math::to_var(system_array[0]);
  vector<vector<var> > biovar_v{ { 1, 1 } };
  vector<vector<var> > tlag_v{ { 0, 0 } };

  rate[0] += h;
  x1 = torsten::linOdeModel(time, amt, rate, ii, evid, cmt, addl, ss,
                            system_array, biovar, tlag);
  rate[0] -= 2 * h;
  x2 = torsten::linOdeModel(time, amt, rate, ii, evid, cmt, addl, ss,
                            system_array, biovar, tlag);
  x = torsten::linOdeModel(time, amt, rate_v, ii, evid, cmt, addl, ss,
                           system_array, biovar, tlag);

#ifndef LINODE_ONECPT_RATE_VAR_TEST
#define LINODE_ONECPT_RATE_VAR_TEST(SYS, BIOVAR, TLAG)                  \
  x = torsten::linOdeModel(time, amt, rate_v, ii, evid, cmt, addl, ss,  \
                           SYS, BIOVAR, TLAG);                          \
  test_it(x);

  LINODE_ONECPT_RATE_VAR_TEST(system_array   , biovar   , tlag);
  LINODE_ONECPT_RATE_VAR_TEST(system_array_v , biovar   , tlag);
  LINODE_ONECPT_RATE_VAR_TEST(system_array_v , biovar_v , tlag);
  LINODE_ONECPT_RATE_VAR_TEST(system_array_v , biovar_v , tlag_v);
  LINODE_ONECPT_RATE_VAR_TEST(system_array   , biovar_v , tlag);
  LINODE_ONECPT_RATE_VAR_TEST(system_array   , biovar_v , tlag_v);
  LINODE_ONECPT_RATE_VAR_TEST(system_array   , biovar   , tlag_v);
  LINODE_ONECPT_RATE_VAR_TEST(system_array_v , biovar_v , tlag_v);

#undef LINODE_ONECPT_RATE_VAR_TEST
#endif
}

class TorstenLinODETwoCptTest : public testing::Test {
  void SetUp() {
    // make sure memory's clean before starting each test
    stan::math::recover_memory();
  }
public:
  TorstenLinODETwoCptTest() :
    CL(5),
    Q(8),
    V2(20),
    V3(70),
    ka(1.2),
    k10(CL / V2),
    k12(Q / V2),
    k21(Q / V3),
    system(3, 3),
    system_array(1),
    biovar{ { 1, 1, 1 } },
    tlag{ { 0, 0, 0 } },
    time(10, 0.0),
    amt(10, 0),
    rate(10, 0),
    cmt(10, 3),
    evid(10, 0),
    ii(10, 0),
    addl(10, 0),
    ss(10, 0)
  {
    system <<
      -ka , 0           , 0,
      ka  , -(k10 + k12), k21,
      0   , k12         , -k21;
    system_array[0] = system;

    time[0] = 0;
    for(int i = 1; i < 9; i++) time[i] = time[i - 1] + 0.25;
    time[9] = 4.0;

    amt[0] = 1000;
    cmt[0] = 1;
    evid[0] = 1;
    ii[0] = 12;
    addl[0] = 10;
    SetUp();
  }

  double CL, Q, V2, V3, ka, k10, k12, k21;
  Matrix<double, Dynamic, Dynamic> system;
  vector<Matrix<double, Dynamic, Dynamic> > system_array;
  vector<vector<double> > biovar;
  vector<vector<double> > tlag;
  vector<double> time;
  vector<double> amt;
  vector<double> rate;
  vector<int> cmt;
  vector<int> evid;
  vector<double> ii;
  vector<int> addl;
  vector<int> ss;
};

TEST_F(TorstenLinODETwoCptTest, MultipleDoses) {
  addl[0] = 14;
  Matrix<double, Dynamic, Dynamic> x;
  x = torsten::linOdeModel(time, amt, rate, ii, evid, cmt, addl, ss,
                           system_array, biovar, tlag);

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

  expect_matrix_eq(amounts, x);

}

TEST_F(TorstenLinODETwoCptTest, steady_state) {
  time[0] = 0.0;
  time[1] = 0.0;
  for(int i = 2; i < 10; i++) time[i] = time[i - 1] + 5;

  amt[0] = 1200;
  ss[0] = 1;

  Matrix<double, Dynamic, Dynamic> x;
  x = torsten::linOdeModel(time, amt, rate, ii, evid, cmt, addl, ss,
                           system_array, biovar, tlag);

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

  for(int i = 0; i < amounts.rows(); i++) {
    EXPECT_NEAR(amounts(i, 0), x(i, 0), std::max(amounts(i, 0), x(i, 0)) * 1e-6);
    EXPECT_NEAR(amounts(i, 1), x(i, 1), std::max(amounts(i, 1), x(i, 1)) * 1e-6);
  }
}

TEST_F(TorstenLinODETwoCptTest, rate) {
  using std::vector;

  V2 = 35;
  V3 = 105;
  k10 = CL / V2;
  k12 = Q / V2;
  k21 = Q / V3;

  amt[0] = 1200;
  rate[0] = 1200;
  addl[0] = 14;

  system <<
    -ka , 0           , 0,
    ka  , -(k10 + k12), k21,
    0   , k12         , -k21;
  system_array[0] = system;

  Matrix<double, Dynamic, Dynamic> x;
  x = torsten::linOdeModel(time, amt, rate, ii, evid, cmt, addl, ss,
                           system_array, biovar, tlag);

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

  expect_matrix_eq(amounts, x);
}

TEST_F(TorstenLinODETwoCptTest, rate_var) {
  using std::vector;
  using stan::math::var;

  V2 = 35;
  V3 = 105;
  k10 = CL / V2;
  k12 = Q / V2;
  k21 = Q / V3;

  amt[0] = 1200;
  rate[0] = 1200;
  const double t_cutoff = amt[0]/rate[0];
  std::vector<stan::math::var> rate_v {stan::math::to_var(rate)};

  const double h = 0.001;
  Matrix<var, Dynamic, Dynamic> x;
  Matrix<double, Dynamic, Dynamic> x1, x2;

  system <<
    -ka , 0           , 0,
    ka  , -(k10 + k12), k21,
    0   , k12         , -k21;
  system_array[0] = system;

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

  std::vector<double> g;
  const double tol = 1e-6;
  auto test_it = [&](Matrix<var, Dynamic, Dynamic>& res) {
    // test val
    expect_matrix_eq(amounts, stan::math::value_of(res));

    // test adj
    for (int i = 0; i < res.rows(); ++i) {
      res(i, 1).grad(rate_v, g);
      stan::math::set_zero_all_adjoints();
      EXPECT_NEAR(g[0], (x1(i, 1) - x2(i, 1))/(2 * h), tol);

      // skip the non-smooth location
      res(i, 0).grad(rate_v, g);
      stan::math::set_zero_all_adjoints();
      if (time[i] != t_cutoff) {
        EXPECT_NEAR(g[0], (x1(i, 0) - x2(i, 0))/(2 * h), tol);
      }
    }
  };

  vector<Matrix<var, Dynamic, Dynamic> > system_array_v(1);
  system_array_v[0] = stan::math::to_var(system_array[0]);
  vector<vector<var> > biovar_v{ { 1, 1, 1 } };
  vector<vector<var> > tlag_v{ { 0, 0, 0 } };

  rate[0] += h;
  x1 = torsten::linOdeModel(time, amt, rate, ii, evid, cmt, addl, ss,
                            system_array, biovar, tlag);
  rate[0] -= 2 * h;
  x2 = torsten::linOdeModel(time, amt, rate, ii, evid, cmt, addl, ss,
                            system_array, biovar, tlag);
  x = torsten::linOdeModel(time, amt, rate_v, ii, evid, cmt, addl, ss,
                           system_array, biovar, tlag);


#ifndef LINODE_TWOCPT_RATE_VAR_TEST
#define LINODE_TWOCPT_RATE_VAR_TEST(SYS, BIOVAR, TLAG)                  \
  x = torsten::linOdeModel(time, amt, rate_v, ii, evid, cmt, addl, ss,  \
                           SYS, BIOVAR, TLAG);                          \
  test_it(x);

  LINODE_TWOCPT_RATE_VAR_TEST(system_array   , biovar   , tlag);
  LINODE_TWOCPT_RATE_VAR_TEST(system_array_v , biovar   , tlag);
  LINODE_TWOCPT_RATE_VAR_TEST(system_array_v , biovar_v , tlag);
  LINODE_TWOCPT_RATE_VAR_TEST(system_array_v , biovar_v , tlag_v);
  LINODE_TWOCPT_RATE_VAR_TEST(system_array   , biovar_v , tlag);
  LINODE_TWOCPT_RATE_VAR_TEST(system_array   , biovar_v , tlag_v);
  LINODE_TWOCPT_RATE_VAR_TEST(system_array   , biovar   , tlag_v);
  LINODE_TWOCPT_RATE_VAR_TEST(system_array_v , biovar_v , tlag_v);

#undef LINODE_TWOCPT_RATE_VAR_TEST
#endif

}
