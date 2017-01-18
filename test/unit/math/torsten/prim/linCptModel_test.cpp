#include <stan/math/rev/mat.hpp>
#include <gtest/gtest.h>
#include <test/unit/math/prim/mat/fun/expect_matrix_eq.hpp>
#include <test/unit/math/torsten/prim/util_linOdeModel.hpp>

using std::vector;
using Eigen::Matrix;
using Eigen::Dynamic;

TEST(Torsten, LinCpt_OneSS) {

	double CL = 10, Vc = 80, ka = 1.2, k10 = CL / Vc;
	Matrix<double, Dynamic, Dynamic> system(2, 2);
	system << -ka, 0, ka, -k10;
	vector<Matrix<double, Dynamic, Dynamic> > system_array(1, system); 

	vector<vector<double> > pMatrix(1);
	pMatrix[0].resize(4);
	pMatrix[0][0] = 1; // F1
	pMatrix[0][1] = 1; // F2
	pMatrix[0][2] = 0; // tlag1
	pMatrix[0][3] = 0; // tlag2

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
	x = linCptModel(system_array, pMatrix, time, amt, rate, ii, evid, cmt, addl, ss);

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

	double diff = 1e-8, diff2 = 1e-4;
	test_linOdeModel(system_array, pMatrix, time, amt, rate, ii, evid, cmt, addl,
		                 ss, diff, diff2);
}

TEST(Torsten, LinCpt_OneSS_overloads) {

	double CL = 10, Vc = 80, ka = 1.2, k10 = CL / Vc;
	Matrix<double, Dynamic, Dynamic> system(2, 2);
	system << -ka, 0, ka, -k10;
	vector<Matrix<double, Dynamic, Dynamic> > system_array(1, system); 

	vector<vector<double> > pMatrix(1);
	pMatrix[0].resize(4);
	pMatrix[0][0] = 1; // F1
	pMatrix[0][1] = 1; // F2
	pMatrix[0][2] = 0; // tlag1
	pMatrix[0][3] = 0; // tlag2

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

	Matrix<double, Dynamic, Dynamic> x1, x2, x3;
	x1 = linCptModel(system_array, pMatrix[0],
	  time, amt, rate, ii, evid, cmt, addl, ss);
	x2 = linCptModel(system_array[0], pMatrix,
	  time, amt, rate, ii, evid, cmt, addl, ss);
	x3 = linCptModel(system_array[0], pMatrix[0],
	  time, amt, rate, ii, evid, cmt, addl, ss);

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
		EXPECT_NEAR(amounts(i, 0), x1(i, 0),
		  std::max(amounts(i, 0), x1(i, 0)) * 1e-6);
		EXPECT_NEAR(amounts(i, 1), x1(i, 1),
		  std::max(amounts(i, 1), x1(i, 1)) * 1e-6);
		  
		EXPECT_NEAR(amounts(i, 0), x2(i, 0),
		  std::max(amounts(i, 0), x2(i, 0)) * 1e-6);
		EXPECT_NEAR(amounts(i, 1), x2(i, 1),
		  std::max(amounts(i, 1), x2(i, 1)) * 1e-6);

		EXPECT_NEAR(amounts(i, 0), x3(i, 0),
		  std::max(amounts(i, 0), x3(i, 0)) * 1e-6);
		EXPECT_NEAR(amounts(i, 1), x3(i, 1),
		  std::max(amounts(i, 1), x3(i, 1)) * 1e-6);
	}

	double diff = 1e-8, diff2 = 1e-4;
	test_linOdeModel(system_array, pMatrix, time, amt, rate, ii, evid, cmt, addl,
		             ss, diff, diff2);
}

TEST(Torsten, linCptModel_OneSS_rate) {

	double CL = 10, Vc = 80, ka = 1.2, k10 = CL / Vc;
	Matrix<double, Dynamic, Dynamic> system(2,2);
	system << -ka, 0, ka, -k10;
	vector<Matrix<double, Dynamic, Dynamic> > system_array(1, system);

	vector<vector<double> > pMatrix(1);
	pMatrix[0].resize(4);
	pMatrix[0][0] = 1; // F1
	pMatrix[0][1] = 1; // F2
	pMatrix[0][2] = 0; // tlag1
	pMatrix[0][3] = 0; // tlag2

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
	x = linCptModel(system_array, pMatrix, time, amt, rate, ii, evid, cmt, addl, ss);

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
	test_linOdeModel(system_array, pMatrix, time, amt, rate, ii, evid, cmt, addl,
		             ss, diff, diff2);

}

TEST(Torsten, linOne_MultipleDoses_timePara) {

    int nEvent = 11;
    
	vector<vector<double> > pMatrix(1);
	pMatrix[0].resize(4);
	pMatrix[0][0] = 1; // F1
	pMatrix[0][1] = 1; // F2
	pMatrix[0][2] = 0; // tlag1
	pMatrix[0][3] = 0; // tlag2
    
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
	x = linCptModel(system_array, pMatrix, time, amt, rate, ii, evid, cmt, addl, ss);

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
	test_linOdeModel(system_array, pMatrix, time, amt, rate, ii, evid, cmt, addl,
		             ss, diff, diff2);
}

