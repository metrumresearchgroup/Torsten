#include <stan/math/rev/mat.hpp>  // FIX ME - includes should be more specific
#include <gtest/gtest.h>
#include <test/unit/math/torsten/expect_near_matrix_eq.hpp>
#include <test/unit/math/torsten/expect_matrix_eq.hpp>
#include <test/unit/math/torsten/util_generalOdeModel.hpp>
#include <stan/math/torsten/PKModelOneCpt.hpp>
#include <stan/math/torsten/pk_onecpt_model.hpp>
#include <stan/math/torsten/pk_twocpt_model.hpp>
#include <gtest/gtest.h>
#include <stan/math/rev/mat.hpp>  // FIX ME - include should be more specific
#include <test/unit/math/torsten/expect_matrix_eq.hpp>
#include <test/unit/math/torsten/util_PKModelOneCpt.hpp>
#include <vector>

using std::vector;
using Eigen::Matrix;
using Eigen::Dynamic;

class TorstenPKOneCptTest : public testing::Test {
  void SetUp() {
    // make sure memory's clean before starting each test
    stan::math::recover_memory();
  }
public:
  TorstenPKOneCptTest() :
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
    time[0] = 0;
    for(int i = 1; i < 9; i++) time[i] = time[i - 1] + 0.25;
    time[9] = 4.0;
    amt[0] = 1000;
    cmt[0] = 1;    
    evid[0] = 1;
    ii[0] = 12;
    addl[0] = 14;
    SetUp();
  }

  vector<vector<double> > pMatrix;  // CL, VC, Ka
  const int nCmt;  // F1, F2
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

TEST_F(TorstenPKOneCptTest, MultipleDoses) {
  Matrix<double, Dynamic, Dynamic> x;
  x = torsten::PKModelOneCpt(time, amt, rate, ii, evid, cmt, addl, ss, pMatrix, biovar, tlag);

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

  expect_matrix_eq(amounts, x);

  // Test AutoDiff against FiniteDiff
  test_PKModelOneCpt(time, amt, rate, ii, evid, cmt, addl, ss,
                     pMatrix, biovar, tlag, 1e-8, 1e-4);
}

TEST_F(TorstenPKOneCptTest, MultipleDoses_IV) {
  cmt[0] = 2;  // IV infusion, not absorption from the gut
  Matrix<double, Dynamic, Dynamic> x;
  x = torsten::PKModelOneCpt(time, amt, rate, ii, evid, cmt, addl, ss, pMatrix, biovar, tlag);

  Matrix<double, Dynamic, Dynamic> amounts(10, 2);
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

  expect_matrix_eq(amounts, x);

  // Test AutoDiff against FiniteDiff
  test_PKModelOneCpt(time, amt, rate, ii, evid, cmt, addl, ss,
                     pMatrix, biovar, tlag, 1e-8, 1e-4);
}

TEST_F(TorstenPKOneCptTest, MultipleDoses_lagTime) {
  tlag[0].resize(nCmt);
  tlag[0][0] = 5;  // tlag1
  tlag[0][1] = 0;  // tlag2

  time[0] = 0;
  for(int i = 1; i < 10; i++) time[i] = time[i - 1] + 1;

  amt[0] = 1200;
  
  Matrix<double, Dynamic, Dynamic> x;
  x = torsten::PKModelOneCpt(time, amt, rate, ii, evid, cmt, addl, ss, pMatrix, biovar, tlag);
  
  Matrix<double, Dynamic, Dynamic> amounts(10, 2);
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

  expect_matrix_eq(amounts, x);

  // Test AutoDiff against FiniteDiff
  test_PKModelOneCpt(time, amt, rate, ii, evid, cmt, addl, ss,
                     pMatrix, biovar, tlag, 1e-8, 1e-4);
}

TEST_F(TorstenPKOneCptTest, MultipleDoses_function_overload) {
	Matrix<double, Dynamic, Dynamic> x_122, x_112, x_111, x_121, x_212,
	                                 x_211, x_221;
	x_122 = torsten::PKModelOneCpt(time, amt, rate, ii, evid, cmt, addl, ss,
                       pMatrix[0], biovar, tlag);
	x_112 = torsten::PKModelOneCpt(time, amt, rate, ii, evid, cmt, addl, ss,
                       pMatrix[0], biovar[0], tlag);
	x_111 = torsten::PKModelOneCpt(time, amt, rate, ii, evid, cmt, addl, ss,
                       pMatrix[0], biovar[0], tlag[0]);
	x_121 = torsten::PKModelOneCpt(time, amt, rate, ii, evid, cmt, addl, ss,
                       pMatrix[0], biovar, tlag[0]);
	x_212 = torsten::PKModelOneCpt(time, amt, rate, ii, evid, cmt, addl, ss,
                       pMatrix, biovar[0], tlag);
	x_211 = torsten::PKModelOneCpt(time, amt, rate, ii, evid, cmt, addl, ss,
                       pMatrix, biovar[0], tlag[0]);
	x_221 = torsten::PKModelOneCpt(time, amt, rate, ii, evid, cmt, addl, ss,
                       pMatrix, biovar, tlag[0]);

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

	expect_matrix_eq(amounts, x_122);
	expect_matrix_eq(amounts, x_112);
	expect_matrix_eq(amounts, x_111);
	expect_matrix_eq(amounts, x_121);
	expect_matrix_eq(amounts, x_212);
	expect_matrix_eq(amounts, x_211);
	expect_matrix_eq(amounts, x_221);

// CHECK - do I need an AD test for every function signature ?
}

TEST_F(TorstenPKOneCptTest, signature) {
  using stan::math::var;
  vector<vector<var> > pMatrix_v{ {10, 80, 1.2 } };
  vector<vector<var> > biovar_v { { 1, 1 } };
  vector<vector<var> > tlag_v{ { 0, 0 } };

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

  vector<Matrix<var, Dynamic, Dynamic> > x_122(7);
  x_122[0] = torsten::PKModelOneCpt(time, amt, rate, ii, evid, cmt, addl, ss,
                           pMatrix_v[0], biovar, tlag);
  x_122[1] = torsten::PKModelOneCpt(time, amt, rate, ii, evid, cmt, addl, ss,
                           pMatrix_v[0], biovar_v, tlag);
  x_122[2] = torsten::PKModelOneCpt(time, amt, rate, ii, evid, cmt, addl, ss,
                           pMatrix_v[0], biovar, tlag_v);
  x_122[3] = torsten::PKModelOneCpt(time, amt, rate, ii, evid, cmt, addl, ss,
                           pMatrix_v[0], biovar_v, tlag_v);
  x_122[4] = torsten::PKModelOneCpt(time, amt, rate, ii, evid, cmt, addl, ss,
                           pMatrix[0], biovar_v, tlag);
  x_122[5] = torsten::PKModelOneCpt(time, amt, rate, ii, evid, cmt, addl, ss,
                           pMatrix[0], biovar_v, tlag_v);
  x_122[6] = torsten::PKModelOneCpt(time, amt, rate, ii, evid, cmt, addl, ss,
                           pMatrix[0], biovar, tlag_v);

  for (size_t i = 0; i < x_122.size(); i++)
    for (int j = 0; j < x_122[i].rows(); j++)
      for (int k = 0; k < x_122[i].cols(); k++)
        EXPECT_FLOAT_EQ(amounts(j, k), x_122[i](j, k).val());


  vector<Matrix<var, Dynamic, Dynamic> > x_112(7);
  x_112[0] = torsten::PKModelOneCpt(time, amt, rate, ii, evid, cmt, addl, ss,
                           pMatrix_v[0], biovar[0], tlag);
  x_112[1] = torsten::PKModelOneCpt(time, amt, rate, ii, evid, cmt, addl, ss,
                           pMatrix_v[0], biovar_v[0], tlag);
  x_112[2] = torsten::PKModelOneCpt(time, amt, rate, ii, evid, cmt, addl, ss,
                           pMatrix_v[0], biovar[0], tlag_v);
  x_112[3] = torsten::PKModelOneCpt(time, amt, rate, ii, evid, cmt, addl, ss,
                           pMatrix_v[0], biovar_v[0], tlag_v);
  x_112[4] = torsten::PKModelOneCpt(time, amt, rate, ii, evid, cmt, addl, ss,
                           pMatrix[0], biovar_v[0], tlag);
  x_112[5] = torsten::PKModelOneCpt(time, amt, rate, ii, evid, cmt, addl, ss,
                           pMatrix[0], biovar_v[0], tlag_v);
  x_112[6] = torsten::PKModelOneCpt(time, amt, rate, ii, evid, cmt, addl, ss,
                           pMatrix[0], biovar[0], tlag_v);

  for (size_t i = 0; i < x_112.size(); i++)
    for (int j = 0; j < x_112[i].rows(); j++)
      for (int k = 0; k < x_112[i].cols(); k++)
        EXPECT_FLOAT_EQ(amounts(j, k), x_112[i](j, k).val());


  vector<Matrix<var, Dynamic, Dynamic> > x_121(7);
  x_121[0] = torsten::PKModelOneCpt(time, amt, rate, ii, evid, cmt, addl, ss,
                           pMatrix_v[0], biovar, tlag[0]);
  x_121[1] = torsten::PKModelOneCpt(time, amt, rate, ii, evid, cmt, addl, ss,
                           pMatrix_v[0], biovar_v, tlag[0]);
  x_121[2] = torsten::PKModelOneCpt(time, amt, rate, ii, evid, cmt, addl, ss,
                           pMatrix_v[0], biovar, tlag_v[0]);
  x_121[3] = torsten::PKModelOneCpt(time, amt, rate, ii, evid, cmt, addl, ss,
                           pMatrix_v[0], biovar_v, tlag_v[0]);
  x_121[4] = torsten::PKModelOneCpt(time, amt, rate, ii, evid, cmt, addl, ss,
                           pMatrix[0], biovar_v, tlag[0]);
  x_121[5] = torsten::PKModelOneCpt(time, amt, rate, ii, evid, cmt, addl, ss,
                           pMatrix[0], biovar_v, tlag_v[0]);
  x_121[6] = torsten::PKModelOneCpt(time, amt, rate, ii, evid, cmt, addl, ss,
                           pMatrix[0], biovar, tlag_v[0]);

  for (size_t i = 0; i < x_121.size(); i++)
    for (int j = 0; j < x_121[i].rows(); j++)
      for (int k = 0; k < x_121[i].cols(); k++)
        EXPECT_FLOAT_EQ(amounts(j, k), x_121[i](j, k).val());


  vector<Matrix<var, Dynamic, Dynamic> > x_111(7);
  x_111[0] = torsten::PKModelOneCpt(time, amt, rate, ii, evid, cmt, addl, ss,
                           pMatrix_v[0], biovar[0], tlag[0]);
  x_111[1] = torsten::PKModelOneCpt(time, amt, rate, ii, evid, cmt, addl, ss,
                           pMatrix_v[0], biovar_v[0], tlag[0]);
  x_111[2] = torsten::PKModelOneCpt(time, amt, rate, ii, evid, cmt, addl, ss,
                           pMatrix_v[0], biovar[0], tlag_v[0]);
  x_111[3] = torsten::PKModelOneCpt(time, amt, rate, ii, evid, cmt, addl, ss,
                           pMatrix_v[0], biovar_v[0], tlag_v[0]);
  x_111[4] = torsten::PKModelOneCpt(time, amt, rate, ii, evid, cmt, addl, ss,
                           pMatrix[0], biovar_v[0], tlag[0]);
  x_111[5] = torsten::PKModelOneCpt(time, amt, rate, ii, evid, cmt, addl, ss,
                           pMatrix[0], biovar_v[0], tlag_v[0]);
  x_111[6] = torsten::PKModelOneCpt(time, amt, rate, ii, evid, cmt, addl, ss,
                           pMatrix[0], biovar[0], tlag_v[0]);

  for (size_t i = 0; i < x_111.size(); i++)
    for (int j = 0; j < x_111[i].rows(); j++)
      for (int k = 0; k < x_111[i].cols(); k++)
        EXPECT_FLOAT_EQ(amounts(j, k), x_111[i](j, k).val());


  vector<Matrix<var, Dynamic, Dynamic> > x_212(7);
  x_212[0] = torsten::PKModelOneCpt(time, amt, rate, ii, evid, cmt, addl, ss,
                           pMatrix_v, biovar[0], tlag);
  x_212[1] = torsten::PKModelOneCpt(time, amt, rate, ii, evid, cmt, addl, ss,
                           pMatrix_v, biovar_v[0], tlag);
  x_212[2] = torsten::PKModelOneCpt(time, amt, rate, ii, evid, cmt, addl, ss,
                           pMatrix_v, biovar[0], tlag_v);
  x_212[3] = torsten::PKModelOneCpt(time, amt, rate, ii, evid, cmt, addl, ss,
                           pMatrix_v, biovar_v[0], tlag_v);
  x_212[4] = torsten::PKModelOneCpt(time, amt, rate, ii, evid, cmt, addl, ss,
                           pMatrix, biovar_v[0], tlag);
  x_212[5] = torsten::PKModelOneCpt(time, amt, rate, ii, evid, cmt, addl, ss,
                           pMatrix, biovar_v[0], tlag_v);
  x_212[6] = torsten::PKModelOneCpt(time, amt, rate, ii, evid, cmt, addl, ss,
                           pMatrix, biovar[0], tlag_v);

  for (size_t i = 0; i < x_212.size(); i++)
    for (int j = 0; j < x_212[i].rows(); j++)
      for (int k = 0; k < x_212[i].cols(); k++)
        EXPECT_FLOAT_EQ(amounts(j, k), x_212[i](j, k).val());


  vector<Matrix<var, Dynamic, Dynamic> > x_211(7);
  x_211[0] = torsten::PKModelOneCpt(time, amt, rate, ii, evid, cmt, addl, ss,
                           pMatrix_v, biovar[0], tlag[0]);
  x_211[1] = torsten::PKModelOneCpt(time, amt, rate, ii, evid, cmt, addl, ss,
                           pMatrix_v, biovar_v[0], tlag[0]);
  x_211[2] = torsten::PKModelOneCpt(time, amt, rate, ii, evid, cmt, addl, ss,
                           pMatrix_v, biovar[0], tlag_v[0]);
  x_211[3] = torsten::PKModelOneCpt(time, amt, rate, ii, evid, cmt, addl, ss,
                           pMatrix_v, biovar_v[0], tlag_v[0]);
  x_211[4] = torsten::PKModelOneCpt(time, amt, rate, ii, evid, cmt, addl, ss,
                           pMatrix, biovar_v[0], tlag[0]);
  x_211[5] = torsten::PKModelOneCpt(time, amt, rate, ii, evid, cmt, addl, ss,
                           pMatrix, biovar_v[0], tlag_v[0]);
  x_211[6] = torsten::PKModelOneCpt(time, amt, rate, ii, evid, cmt, addl, ss,
                           pMatrix, biovar[0], tlag_v[0]);

  for (size_t i = 0; i < x_211.size(); i++)
    for (int j = 0; j < x_211[i].rows(); j++)
      for (int k = 0; k < x_211[i].cols(); k++)
        EXPECT_FLOAT_EQ(amounts(j, k), x_211[i](j, k).val());


  vector<Matrix<var, Dynamic, Dynamic> > x_221(7);
  x_221[0] = torsten::PKModelOneCpt(time, amt, rate, ii, evid, cmt, addl, ss,
                           pMatrix_v, biovar, tlag[0]);
  x_221[1] = torsten::PKModelOneCpt(time, amt, rate, ii, evid, cmt, addl, ss,
                           pMatrix_v, biovar_v, tlag[0]);
  x_221[2] = torsten::PKModelOneCpt(time, amt, rate, ii, evid, cmt, addl, ss,
                           pMatrix_v, biovar, tlag_v[0]);
  x_221[3] = torsten::PKModelOneCpt(time, amt, rate, ii, evid, cmt, addl, ss,
                           pMatrix_v, biovar_v, tlag_v[0]);
  x_221[4] = torsten::PKModelOneCpt(time, amt, rate, ii, evid, cmt, addl, ss,
                           pMatrix, biovar_v, tlag[0]);
  x_221[5] = torsten::PKModelOneCpt(time, amt, rate, ii, evid, cmt, addl, ss,
                           pMatrix, biovar_v, tlag_v[0]);
  x_221[6] = torsten::PKModelOneCpt(time, amt, rate, ii, evid, cmt, addl, ss,
                           pMatrix, biovar, tlag_v[0]);

  for (size_t i = 0; i < x_221.size(); i++)
    for (int j = 0; j < x_221[i].rows(); j++)
      for (int k = 0; k < x_221[i].cols(); k++)
        EXPECT_FLOAT_EQ(amounts(j, k), x_221[i](j, k).val());

  // CHECK - do I need an AD test for every function signature ?
}

TEST_F(TorstenPKOneCptTest, Steady_State) {
  time[0] = 0.0;
  time[1] = 0.0;
  for(int i = 2; i < 10; i++) time[i] = time[i - 1] + 5;
  amt[0] = 1200;
  addl[0] = 10;
  ss[0] = 1;

  Matrix<double, Dynamic, Dynamic> x;
  x = torsten::PKModelOneCpt(time, amt, rate, ii, evid, cmt, addl, ss,
                             pMatrix, biovar, tlag);

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
  test_PKModelOneCpt(time, amt, rate, ii, evid, cmt, addl, ss,
                     pMatrix, biovar, tlag, 1e-8, 1e-4);
}

TEST_F(TorstenPKOneCptTest, Steady_State_rate) {
  time[0] = 0.0;
  time[1] = 0.0;
  for(int i = 2; i < 10; i++) time[i] = time[i - 1] + 5;
  amt[0] = 1200;
  rate[0] = 150;
  addl[0] = 10;
  ss[0] = 1;
  Matrix<double, Dynamic, Dynamic> x;
  x = torsten::PKModelOneCpt(time, amt, rate, ii, evid, cmt, addl, ss,
                             pMatrix, biovar, tlag);

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

  // Test AutoDiff against FiniteDiff
  test_PKModelOneCpt(time, amt, rate, ii, evid, cmt, addl, ss,
                     pMatrix, biovar, tlag, 1e-8, 1e-4);
}

TEST(TorstenPKOneCpt, events_specific_data) {

  int nEvent = 11;
	vector<vector<double> > pMatrix(nEvent);
	for (int i = 0; i < nEvent; i++) {
	  pMatrix[i].resize(3);
	  if (i < 6) pMatrix[i][0] = 10; // CL
	  else pMatrix[i][0] = 50; // CL is piece-wise contant
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

	Matrix<double, Dynamic, Dynamic> x;
	x = torsten::PKModelOneCpt(time, amt, rate, ii, evid, cmt, addl, ss,
                   pMatrix, biovar, tlag);

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

	// Test AutoDiff against FiniteDiff
    test_PKModelOneCpt(time, amt, rate, ii, evid, cmt, addl, ss,
                       pMatrix, biovar, tlag, 1e-8, 1e-4);
}

TEST_F(TorstenPKOneCptTest, rate) {
  using std::vector;
  amt[0] = 1200;
  rate[0] = 1200;

  Matrix<double, Dynamic, Dynamic> x;
  x = torsten::PKModelOneCpt(time, amt, rate, ii, evid, cmt, addl, ss,
                    pMatrix, biovar, tlag);

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

  // Test AutoDiff against FiniteDiff
  test_PKModelOneCpt(time, amt, rate, ii, evid, cmt, addl, ss,
                     pMatrix, biovar, tlag, 1e-8, 5e-4);
}

TEST_F(TorstenPKOneCptTest, multiple_trunc_rate_var) {
  using std::vector;
  using stan::math::var;

  amt[0] = 1200;
  rate[0] = 1200;

  const double t_cutoff = amt[0]/rate[0];

  std::vector<stan::math::var> rate_v {stan::math::to_var(rate)};

  // grad: test against 2nd-order finite diff
  const double h = 0.001;
  Matrix<double, Dynamic, Dynamic> x1, x2;
  rate[0] += h;
  x1 = torsten::PKModelOneCpt(time, amt, rate, ii, evid, cmt, addl, ss,
                              pMatrix, biovar, tlag);
  rate[0] -= 2 * h;
  x2 = torsten::PKModelOneCpt(time, amt, rate, ii, evid, cmt, addl, ss,
                              pMatrix, biovar, tlag);

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
  auto test_it = [&](Matrix<var, Dynamic, Dynamic>& res) {
    // test val
    expect_matrix_eq(amounts, stan::math::value_of(res));

    // test adj
    for (int i = 0; i < res.rows(); ++i) {
      res(i, 1).grad(rate_v, g);
      EXPECT_NEAR(g[0], (x1(i, 1) - x2(i, 1))/(2 * h), 1e-6);
      stan::math::set_zero_all_adjoints();

      // for gut there is a non-smooth point, skip it
      res(i, 0).grad(rate_v, g);
      if(time[i] != t_cutoff)
        EXPECT_NEAR(g[0], (x1(i, 0) - x2(i, 0))/(2 * h), 1e-6);
      stan::math::set_zero_all_adjoints();
    }
  };

  Matrix<var, Dynamic, Dynamic> x;
  vector<vector<var> > pMatrix_v{ {10, 80, 1.2 } };
  vector<vector<var> > biovar_v{ { 1, 1 } };
  vector<vector<var> > tlag_v{ { 0, 0 } };

#ifndef ONECPT_RATE_VAR_TEST
#define ONECPT_RATE_VAR_TEST(PMAT, BIOVAR, TLAG)                        \
  x = torsten::PKModelOneCpt(time, amt, rate_v, ii, evid, cmt, addl, ss, \
                             PMAT, BIOVAR, TLAG);                       \
  test_it(x);

  ONECPT_RATE_VAR_TEST(pMatrix   , biovar   , tlag);
  ONECPT_RATE_VAR_TEST(pMatrix_v , biovar   , tlag);
  ONECPT_RATE_VAR_TEST(pMatrix_v , biovar_v , tlag);
  ONECPT_RATE_VAR_TEST(pMatrix_v , biovar_v , tlag_v);
  ONECPT_RATE_VAR_TEST(pMatrix   , biovar_v , tlag);
  ONECPT_RATE_VAR_TEST(pMatrix   , biovar_v , tlag_v);
  ONECPT_RATE_VAR_TEST(pMatrix   , biovar   , tlag_v);
  ONECPT_RATE_VAR_TEST(pMatrix_v , biovar_v , tlag_v);

#undef ONECPT_RATE_VAR_TEST
#endif

}

TEST_F(TorstenPKOneCptTest, MultipleDoses_IV_rate_var) {
  using stan::math::var;
  cmt[0] = 2;  // IV infusion, not absorption from the gut

  rate[0] = 600;
  std::vector<stan::math::var> rate_v {stan::math::to_var(rate)};

  Matrix<double, Dynamic, Dynamic> amounts(10, 2);
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

  // grad: test against 2nd-order finite diff
  const double h = 0.001;
  Matrix<double, Dynamic, Dynamic> x1, x2;
  rate[0] += h;
  x1 = torsten::PKModelOneCpt(time, amt, rate, ii, evid, cmt, addl, ss,
                              pMatrix, biovar, tlag);
  rate[0] -= 2 * h;
  x2 = torsten::PKModelOneCpt(time, amt, rate, ii, evid, cmt, addl, ss,
                              pMatrix, biovar, tlag);

  std::vector<double> g;
  auto test_it = [&](Matrix<var, Dynamic, Dynamic>& res) {
    // test val
    expect_matrix_eq(amounts, stan::math::value_of(res));

    // test adj
    for (int i = 0; i < res.rows(); ++i) {
      res(i, 1).grad(rate_v, g);
      EXPECT_NEAR(g[0], (x1(i, 1) - x2(i, 1))/(2 * h), 1e-6);
      stan::math::set_zero_all_adjoints();

      // for gut there is a non-smooth point, skip it
      res(i, 0).grad(rate_v, g);
      EXPECT_NEAR(g[0], (x1(i, 0) - x2(i, 0))/(2 * h), 1e-6);
      stan::math::set_zero_all_adjoints();
    }
  };

  Matrix<var, Dynamic, Dynamic> x;
  vector<vector<var> > pMatrix_v{ {10, 80, 1.2 } };
  vector<vector<var> > biovar_v{ { 1, 1 } };
  vector<vector<var> > tlag_v{ { 0, 0 } };

#ifndef ONECPT_RATE_VAR_TEST
#define ONECPT_RATE_VAR_TEST(PMAT, BIOVAR, TLAG)                        \
  x = torsten::PKModelOneCpt(time, amt, rate_v, ii, evid, cmt, addl, ss, \
                             PMAT, BIOVAR, TLAG);                       \
  test_it(x);

  ONECPT_RATE_VAR_TEST(pMatrix   , biovar   , tlag);
  ONECPT_RATE_VAR_TEST(pMatrix_v , biovar   , tlag);
  ONECPT_RATE_VAR_TEST(pMatrix_v , biovar_v , tlag);
  ONECPT_RATE_VAR_TEST(pMatrix_v , biovar_v , tlag_v);
  ONECPT_RATE_VAR_TEST(pMatrix   , biovar_v , tlag);
  ONECPT_RATE_VAR_TEST(pMatrix   , biovar_v , tlag_v);
  ONECPT_RATE_VAR_TEST(pMatrix   , biovar   , tlag_v);
  ONECPT_RATE_VAR_TEST(pMatrix_v , biovar_v , tlag_v);

#undef ONECPT_RATE_VAR_TEST
#endif
}

/*
TEST(Torsten, PKModelOneCptModel_SS_rate_2) {
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

  Matrix<double, Dynamic, Dynamic> x;
  x = PKModelOneCpt(time, amt, rate, ii, evid, cmt, addl, ss,
                    pMatrix, biovar, tlag);

  std::cout << x << std::endl;

  Matrix<double, Dynamic, Dynamic> amounts(10, 2);
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
  // test_PKModelOneCpt(time, amt, rate, ii, evid, cmt, addl, ss,
  //                   pMatrix, biovar, tlag, 1e-8, 5e-4);
} */
