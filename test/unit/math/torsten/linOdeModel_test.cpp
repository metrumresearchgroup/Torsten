#include <gtest/gtest.h>
#include <stan/math/rev/mat.hpp>  // FIX ME - include should be more specific
#include <test/unit/math/prim/mat/fun/expect_matrix_eq.hpp>
#include <test/unit/math/torsten/util_linOdeModel.hpp>
#include <vector>

using std::vector;
using Eigen::Matrix;
using Eigen::Dynamic;

TEST(Torsten, LinCpt_OneSS) {

  double CL = 10, Vc = 80, ka = 1.2, k10 = CL / Vc;
	Matrix<double, Dynamic, Dynamic> system(2, 2);
	system << -ka, 0, ka, -k10;
	vector<Matrix<double, Dynamic, Dynamic> > system_array(1, system); 

	vector<vector<double> > biovar(1);
	biovar[0].resize(2);
	biovar[0][0] = 1;  // F1
	biovar[0][1] = 1;  // F2
	
	vector<vector<double> > tlag(1);
	tlag[0].resize(2);
	tlag[0][0] = 0;  // tlag1
	tlag[0][1] = 0;  // tlag2

	vector<double> time(10);
	time[0] = 0.0;
	time[1] = 0.0;
	for(int i = 2; i < 10; i++) time[i] = time[i - 1] + 5;

	vector<double> amt(10, 0);
	amt[0] = 1200;

	vector<double> rate(10, 0);

	vector<int> cmt(10, 2);
	cmt[0] = 1;

	vector<int> evid(10, 0);
	evid[0] = 1;

	vector<double> ii(10, 0);
	ii[0] = 12;

	vector<int> addl(10, 0);
	addl[0] = 10;

	vector<int> ss(10, 0);
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

TEST(Torsten, LinCpt_OneSS_overloads) {

  double CL = 10, Vc = 80, ka = 1.2, k10 = CL / Vc;
  Matrix<double, Dynamic, Dynamic> system(2, 2);
  system << -ka, 0, ka, -k10;
  vector<Matrix<double, Dynamic, Dynamic> > system_array(1, system); 
  
  vector<vector<double> > biovar(1);
  biovar[0].resize(2);
  biovar[0][0] = 1;  // F1
  biovar[0][1] = 1;  // F2
  
  vector<vector<double> > tlag(1);
  tlag[0].resize(2);
  tlag[0][0] = 0;  // tlag1
  tlag[0][1] = 0;  // tlag2
  
  vector<double> time(10);
  time[0] = 0.0;
  time[1] = 0.0;
  for(int i = 2; i < 10; i++) time[i] = time[i - 1] + 5;

	vector<double> amt(10, 0);
	amt[0] = 1200;

	vector<double> rate(10, 0);

	vector<int> cmt(10, 2);
	cmt[0] = 1;

	vector<int> evid(10, 0);
	evid[0] = 1;

	vector<double> ii(10, 0);
	ii[0] = 12;

	vector<int> addl(10, 0);
	addl[0] = 10;

	vector<int> ss(10, 0);
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

TEST(Torsten, linOdeModel_signature_test) {
  using stan::math::var;
  
  double CL = 10, Vc = 80, ka = 1.2, k10 = CL / Vc;
  Matrix<double, Dynamic, Dynamic> system(2, 2);
  system << -ka, 0, ka, -k10;
  vector<Matrix<double, Dynamic, Dynamic> > system_array(1, system);
  
  vector<vector<double> > biovar(1);
  biovar[0].resize(2);
  biovar[0][0] = 1;  // F1
  biovar[0][1] = 1;  // F2

  vector<vector<double> > tlag(1);
  tlag[0].resize(2);
  tlag[0][0] = 0;  // tlag1
  tlag[0][1] = 0;  // tlag2

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

  vector<double> time(10);
  time[0] = 0.0;
  time[1] = 0.0;
  for(int i = 2; i < 10; i++) time[i] = time[i - 1] + 5;
  
  vector<double> amt(10, 0);
  amt[0] = 1200;
  
  vector<double> rate(10, 0);
  
  vector<int> cmt(10, 2);
  cmt[0] = 1;
  
  vector<int> evid(10, 0);
  evid[0] = 1;
  
  vector<double> ii(10, 0);
  ii[0] = 12;
  
  vector<int> addl(10, 0);
  addl[0] = 10;
  
  vector<int> ss(10, 0);
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
  
  vector<Matrix<var, Dynamic, Dynamic> > x_122(7);
  x_122[0] = torsten::linOdeModel(time, amt, rate, ii, evid, cmt, addl, ss,
                           system_array_v[0], biovar, tlag);
  x_122[1] = torsten::linOdeModel(time, amt, rate, ii, evid, cmt, addl, ss,
                           system_array_v[0], biovar_v, tlag);
  x_122[2] = torsten::linOdeModel(time, amt, rate, ii, evid, cmt, addl, ss,
                           system_array_v[0], biovar, tlag_v);
  x_122[3] = torsten::linOdeModel(time, amt, rate, ii, evid, cmt, addl, ss,
                           system_array_v[0], biovar_v, tlag_v);
  x_122[4] = torsten::linOdeModel(time, amt, rate, ii, evid, cmt, addl, ss,
                           system_array[0], biovar_v, tlag);
  x_122[5] = torsten::linOdeModel(time, amt, rate, ii, evid, cmt, addl, ss,
                           system_array[0], biovar_v, tlag_v);
  x_122[6] = torsten::linOdeModel(time, amt, rate, ii, evid, cmt, addl, ss,
                           system_array[0], biovar, tlag_v);
  
  for (size_t i = 0; i < x_122.size(); i++)
    for (int j = 0; j < x_122[i].rows(); j++)
      for (int k = 0; k < x_122[i].cols(); k++)
        EXPECT_NEAR(amounts(j, k), x_122[i](j, k).val(),
          std::max(amounts(j, k), x_122[i](j, k).val()) * 1e-6);
  
 
  vector<Matrix<var, Dynamic, Dynamic> > x_112(7);
  x_112[0] = torsten::linOdeModel(time, amt, rate, ii, evid, cmt, addl, ss,
                           system_array_v[0], biovar[0], tlag);
  x_112[1] = torsten::linOdeModel(time, amt, rate, ii, evid, cmt, addl, ss,
                           system_array_v[0], biovar_v[0], tlag);
  x_112[2] = torsten::linOdeModel(time, amt, rate, ii, evid, cmt, addl, ss,
                           system_array_v[0], biovar[0], tlag_v);
  x_112[3] = torsten::linOdeModel(time, amt, rate, ii, evid, cmt, addl, ss,
                           system_array_v[0], biovar_v[0], tlag_v);
  x_112[4] = torsten::linOdeModel(time, amt, rate, ii, evid, cmt, addl, ss,
                           system_array[0], biovar_v[0], tlag);
  x_112[5] = torsten::linOdeModel(time, amt, rate, ii, evid, cmt, addl, ss,
                           system_array[0], biovar_v[0], tlag_v);
  x_112[6] = torsten::linOdeModel(time, amt, rate, ii, evid, cmt, addl, ss,
                           system_array[0], biovar[0], tlag_v);
  
  for (size_t i = 0; i < x_112.size(); i++)
    for (int j = 0; j < x_112[i].rows(); j++)
      for (int k = 0; k < x_112[i].cols(); k++)
        EXPECT_NEAR(amounts(j, k), x_112[i](j, k).val(),
          std::max(amounts(j, k), x_112[i](j, k).val()) * 1e-6);
  
  
  vector<Matrix<var, Dynamic, Dynamic> > x_121(7);
  x_121[0] = torsten::linOdeModel(time, amt, rate, ii, evid, cmt, addl, ss,
                           system_array_v[0], biovar, tlag[0]);
  x_121[1] = torsten::linOdeModel(time, amt, rate, ii, evid, cmt, addl, ss,
                           system_array_v[0], biovar_v, tlag[0]);
  x_121[2] = torsten::linOdeModel(time, amt, rate, ii, evid, cmt, addl, ss,
                           system_array_v[0], biovar, tlag_v[0]);
  x_121[3] = torsten::linOdeModel(time, amt, rate, ii, evid, cmt, addl, ss,
                           system_array_v[0], biovar_v, tlag_v[0]);
  x_121[4] = torsten::linOdeModel(time, amt, rate, ii, evid, cmt, addl, ss,
                           system_array[0], biovar_v, tlag[0]);
  x_121[5] = torsten::linOdeModel(time, amt, rate, ii, evid, cmt, addl, ss,
                           system_array[0], biovar_v, tlag_v[0]);
  x_121[6] = torsten::linOdeModel(time, amt, rate, ii, evid, cmt, addl, ss,
                           system_array[0], biovar, tlag_v[0]);
  
  for (size_t i = 0; i < x_121.size(); i++)
    for (int j = 0; j < x_121[i].rows(); j++)
      for (int k = 0; k < x_121[i].cols(); k++)
        EXPECT_NEAR(amounts(j, k), x_121[i](j, k).val(),
          std::max(amounts(j, k), x_121[i](j, k).val()) * 1e-6);
  
  
  vector<Matrix<var, Dynamic, Dynamic> > x_111(7);
  x_111[0] = torsten::linOdeModel(time, amt, rate, ii, evid, cmt, addl, ss,
                           system_array_v[0], biovar[0], tlag[0]);
  x_111[1] = torsten::linOdeModel(time, amt, rate, ii, evid, cmt, addl, ss,
                           system_array_v[0], biovar_v[0], tlag[0]);
  x_111[2] = torsten::linOdeModel(time, amt, rate, ii, evid, cmt, addl, ss,
                           system_array_v[0], biovar[0], tlag_v[0]);
  x_111[3] = torsten::linOdeModel(time, amt, rate, ii, evid, cmt, addl, ss,
                           system_array_v[0], biovar_v[0], tlag_v[0]);
  x_111[4] = torsten::linOdeModel(time, amt, rate, ii, evid, cmt, addl, ss,
                           system_array[0], biovar_v[0], tlag[0]);
  x_111[5] = torsten::linOdeModel(time, amt, rate, ii, evid, cmt, addl, ss,
                           system_array[0], biovar_v[0], tlag_v[0]);
  x_111[6] = torsten::linOdeModel(time, amt, rate, ii, evid, cmt, addl, ss,
                           system_array[0], biovar[0], tlag_v[0]);
  
  for (size_t i = 0; i < x_111.size(); i++)
    for (int j = 0; j < x_111[i].rows(); j++)
      for (int k = 0; k < x_111[i].cols(); k++)
        EXPECT_NEAR(amounts(j, k), x_111[i](j, k).val(),
          std::max(amounts(j, k), x_111[i](j, k).val()) * 1e-6);
  
  
  vector<Matrix<var, Dynamic, Dynamic> > x_212(7);
  x_212[0] = torsten::linOdeModel(time, amt, rate, ii, evid, cmt, addl, ss,
                           system_array_v, biovar[0], tlag);
  x_212[1] = torsten::linOdeModel(time, amt, rate, ii, evid, cmt, addl, ss,
                           system_array_v, biovar_v[0], tlag);
  x_212[2] = torsten::linOdeModel(time, amt, rate, ii, evid, cmt, addl, ss,
                           system_array_v, biovar[0], tlag_v);
  x_212[3] = torsten::linOdeModel(time, amt, rate, ii, evid, cmt, addl, ss,
                           system_array_v, biovar_v[0], tlag_v);
  x_212[4] = torsten::linOdeModel(time, amt, rate, ii, evid, cmt, addl, ss,
                           system_array, biovar_v[0], tlag);
  x_212[5] = torsten::linOdeModel(time, amt, rate, ii, evid, cmt, addl, ss,
                           system_array, biovar_v[0], tlag_v);
  x_212[6] = torsten::linOdeModel(time, amt, rate, ii, evid, cmt, addl, ss,
                           system_array, biovar[0], tlag_v);
  
  for (size_t i = 0; i < x_212.size(); i++)
    for (int j = 0; j < x_212[i].rows(); j++)
      for (int k = 0; k < x_212[i].cols(); k++)
        EXPECT_NEAR(amounts(j, k), x_212[i](j, k).val(),
          std::max(amounts(j, k), x_212[i](j, k).val()) * 1e-6);
  
  
  vector<Matrix<var, Dynamic, Dynamic> > x_211(7);
  x_211[0] = torsten::linOdeModel(time, amt, rate, ii, evid, cmt, addl, ss,
                           system_array_v, biovar[0], tlag[0]);
  x_211[1] = torsten::linOdeModel(time, amt, rate, ii, evid, cmt, addl, ss,
                           system_array_v, biovar_v[0], tlag[0]);
  x_211[2] = torsten::linOdeModel(time, amt, rate, ii, evid, cmt, addl, ss,
                           system_array_v, biovar[0], tlag_v[0]);
  x_211[3] = torsten::linOdeModel(time, amt, rate, ii, evid, cmt, addl, ss,
                           system_array_v, biovar_v[0], tlag_v[0]);
  x_211[4] = torsten::linOdeModel(time, amt, rate, ii, evid, cmt, addl, ss,
                           system_array, biovar_v[0], tlag[0]);
  x_211[5] = torsten::linOdeModel(time, amt, rate, ii, evid, cmt, addl, ss,
                           system_array, biovar_v[0], tlag_v[0]);
  x_211[6] = torsten::linOdeModel(time, amt, rate, ii, evid, cmt, addl, ss,
                           system_array, biovar[0], tlag_v[0]);
  
  for (size_t i = 0; i < x_211.size(); i++)
    for (int j = 0; j < x_211[i].rows(); j++)
      for (int k = 0; k < x_211[i].cols(); k++)
        EXPECT_NEAR(amounts(j, k), x_211[i](j, k).val(),
          std::max(amounts(j, k), x_211[i](j, k).val()) * 1e-6);
  
  
  vector<Matrix<var, Dynamic, Dynamic> > x_221(7);
  x_221[0] = torsten::linOdeModel(time, amt, rate, ii, evid, cmt, addl, ss,
                           system_array_v, biovar, tlag[0]);
  x_221[1] = torsten::linOdeModel(time, amt, rate, ii, evid, cmt, addl, ss,
                           system_array_v, biovar_v, tlag[0]);
  x_221[2] = torsten::linOdeModel(time, amt, rate, ii, evid, cmt, addl, ss,
                           system_array_v, biovar, tlag_v[0]);
  x_221[3] = torsten::linOdeModel(time, amt, rate, ii, evid, cmt, addl, ss,
                           system_array_v, biovar_v, tlag_v[0]);
  x_221[4] = torsten::linOdeModel(time, amt, rate, ii, evid, cmt, addl, ss,
                           system_array, biovar_v, tlag[0]);
  x_221[5] = torsten::linOdeModel(time, amt, rate, ii, evid, cmt, addl, ss,
                           system_array, biovar_v, tlag_v[0]);
  x_221[6] = torsten::linOdeModel(time, amt, rate, ii, evid, cmt, addl, ss,
                           system_array, biovar, tlag_v[0]);

  for (size_t i = 0; i < x_221.size(); i++)
    for (int j = 0; j < x_221[i].rows(); j++)
      for (int k = 0; k < x_221[i].cols(); k++)
        EXPECT_NEAR(amounts(j, k), x_221[i](j, k).val(),
          std::max(amounts(j, k), x_221[i](j, k).val()) * 1e-6);
  
  // CHECK - do I need an AD test for every function signature ?
}

TEST(Torsten, linOdeModel_OneSS_rate) {

	double CL = 10, Vc = 80, ka = 1.2, k10 = CL / Vc;
	Matrix<double, Dynamic, Dynamic> system(2,2);
	system << -ka, 0, ka, -k10;
	vector<Matrix<double, Dynamic, Dynamic> > system_array(1, system);

	vector<vector<double> > biovar(1);
	biovar[0].resize(2);
	biovar[0][0] = 1;  // F1
	biovar[0][1] = 1;  // F2
	
	vector<vector<double> > tlag(1);
	tlag[0].resize(2);
	tlag[0][0] = 0;  // tlag1
	tlag[0][1] = 0;  // tlag2
	
	vector<double> time(10);
	time[0] = 0.0;
	time[1] = 0.0;
	for(int i = 2; i < 10; i++) time[i] = time[i - 1] + 5;

	vector<double> amt(10, 0);
	amt[0] = 1200;

	vector<double> rate(10, 0);
	rate[0] = 150;

	vector<int> cmt(10, 2);
	cmt[0] = 1;

	vector<int> evid(10, 0);
	evid[0] = 1;

	vector<double> ii(10, 0);
	ii[0] = 12;

	vector<int> addl(10, 0);
	addl[0] = 10;

	vector<int> ss(10, 0);
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

TEST(Torsten, linOde_rate) {
  using std::vector;
  
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

  double CL = 10, Vc = 80, ka = 1.2, k10 = CL / Vc;
	Matrix<double, Dynamic, Dynamic> system(2, 2);
	system << -ka, 0, ka, -k10;
	vector<Matrix<double, Dynamic, Dynamic> > system_array(1, system); 

  vector<double> time(10);
  time[0] = 0;
  for(int i = 1; i < 9; i++) time[i] = time[i - 1] + 0.25;
  time[9] = 4.0;

  vector<double> amt(10, 0);
  amt[0] = 1200;

  vector<double> rate(10, 0);
  rate[0] = 1200;

  vector<int> cmt(10, 2);
  cmt[0] = 1;

  vector<int> evid(10, 0);
  evid[0] = 1;

  vector<double> ii(10, 0);
  ii[0] = 12;

  vector<int> addl(10, 0);
  addl[0] = 14;

  vector<int> ss(10, 0);
  
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
