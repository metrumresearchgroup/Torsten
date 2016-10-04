#include <stan/math/rev/mat.hpp>
#include <gtest/gtest.h>
#include <test/unit/math/prim/mat/fun/expect_matrix_eq.hpp>

using std::vector;
using Eigen::Matrix;
using Eigen::Dynamic;


TEST(Torsten, PKModelOneCpt_SingleDose) {

	double CL = 10, Vc = 80, ka = 1.2, k10 = CL / Vc;
	Matrix<double, Dynamic, Dynamic> system(2,2);
	system << -ka, 0,
	          ka, -k10;

	vector<Matrix<double, Dynamic, 1> > pMatrix(1);
	pMatrix[0].resize(7);
	pMatrix[0](0) = 10; // CL
	pMatrix[0](1) = 80; // Vc
	pMatrix[0](2) = 1.2; // ka
	pMatrix[0](3) = 1; // F1
	pMatrix[0](4) = 1; // F2
	pMatrix[0](5) = 0; // tlag1
	pMatrix[0](6) = 0; // tlag2

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

	Matrix<double, Dynamic, Dynamic> x;
	x = linCptModel(system, pMatrix, time, amt, rate, ii, evid, cmt, addl, ss);
	
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
}


TEST(Torsten, PKModelOneCpt_SS) {

	double CL = 10, Vc = 80, ka = 1.2, k10 = CL / Vc;
	Matrix<double, Dynamic, Dynamic> system(2,2);
	system << -ka, 0, ka, -k10;

	vector<Matrix<double, Dynamic, 1> > pMatrix(1);
	pMatrix[0].resize(4);
	pMatrix[0](0) = 1; // F1
	pMatrix[0](1) = 1; // F2
	pMatrix[0](2) = 0; // tlag1
	pMatrix[0](3) = 0; // tlag2

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
	x = linCptModel(system, pMatrix, time, amt, rate, ii, evid, cmt, addl, ss);

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
}


TEST(Torsten, PKModelOneCpt_SS_rate) {

	double CL = 10, Vc = 80, ka = 1.2, k10 = CL / Vc;
	Matrix<double, Dynamic, Dynamic> system(2,2);
	system << -ka, 0, ka, -k10;

	vector<Matrix<double, Dynamic, 1> > pMatrix(1);
	pMatrix[0].resize(4);
	pMatrix[0](0) = 1; // F1
	pMatrix[0](1) = 1; // F2
	pMatrix[0](2) = 0; // tlag1
	pMatrix[0](3) = 0; // tlag2

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
	x = linCptModel(system, pMatrix, time, amt, rate, ii, evid, cmt, addl, ss);

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
}
