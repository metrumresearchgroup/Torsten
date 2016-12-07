#include <gtest/gtest.h>
#include <stan/math/rev/mat.hpp>
#include <test/unit/math/rev/mat/fun/util.hpp>
#include <test/unit/math/prim/mat/fun/expect_matrix_eq.hpp>

using std::vector;
using Eigen::Matrix;
using Eigen::Dynamic;


TEST(Torsten, PKModelTwoCpt_MultipleDoses) {

	vector<vector<AVAR> > pMatrix(1);
	pMatrix[0].resize(11);
	pMatrix[0][0] = 5; // CL
	pMatrix[0][1] = 8; // Q
	pMatrix[0][2] = 20; // Vc
	pMatrix[0][3] = 70; // Vp
	pMatrix[0][4] = 1.2; // ka
	pMatrix[0][5] = 1; // F1
	pMatrix[0][6] = 1; // F2
	pMatrix[0][7] = 1; // F3
	pMatrix[0][8] = 0; // tlag1
	pMatrix[0][9] = 0; // tlag2
	pMatrix[0][10] = 0; // tlag3

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

	stan::math::matrix_v x;
	x = PKModelTwoCpt(pMatrix, time, amt, rate, ii, evid, cmt, addl, ss);
	
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
			   
	for (size_t i = 0; i < amounts.rows(); i++)
      for (size_t j = 0; j < amounts.cols(); j++)
        EXPECT_FLOAT_EQ(amounts(i, j), x(i, j).val());
}


TEST(Torsten, PKModelTwoCpt_MultipleDoses_overload) {

	vector<AVAR> pMatrix(11);
	pMatrix[0] = 5; // CL
	pMatrix[1] = 8; // Q
	pMatrix[2] = 20; // Vc
	pMatrix[3] = 70; // Vp
	pMatrix[4] = 1.2; // ka
	pMatrix[5] = 1; // F1
	pMatrix[6] = 1; // F2
	pMatrix[7] = 1; // F3
	pMatrix[8] = 0; // tlag1
	pMatrix[9] = 0; // tlag2
	pMatrix[10] = 0; // tlag3

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

	stan::math::matrix_v x;
	x = PKModelTwoCpt(pMatrix, time, amt, rate, ii, evid, cmt, addl, ss);
	
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
			   
	for (size_t i = 0; i < amounts.rows(); i++)
      for (size_t j = 0; j < amounts.cols(); j++)
        EXPECT_FLOAT_EQ(amounts(i, j), x(i, j).val());
}


TEST(Torsten, PKModelTwoCpt_SS) {

	vector<vector<AVAR> > pMatrix(1);
	pMatrix[0].resize(11);
	pMatrix[0][0] = 5; // CL
	pMatrix[0][1] = 8; // Q
	pMatrix[0][2] = 20; // Vc
	pMatrix[0][3] = 70; // Vp
	pMatrix[0][4] = 1.2; // ka
	pMatrix[0][5] = 1; // F1
	pMatrix[0][6] = 1; // F2
	pMatrix[0][7] = 1; // F3
	pMatrix[0][8] = 0; // tlag1
	pMatrix[0][9] = 0; // tlag2
	pMatrix[0][10] = 0; // tlag3
	
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

	stan::math::matrix_v x;
	x = PKModelTwoCpt(pMatrix, time, amt, rate, ii, evid, cmt, addl, ss);
	
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

	for (int i = 0; i < amounts.rows(); i++) {
		EXPECT_NEAR(amounts(i, 0), x(i, 0).val(),
		  std::max(amounts(i, 0), x(i, 0).val()) * 1e-6);
		EXPECT_NEAR(amounts(i, 1), x(i, 1).val(),
		  std::max(amounts(i, 1), x(i, 1).val()) * 1e-6);
	}
}

TEST(Torsten, PKModelTwoCpt_SS_rate) {
	
	vector<vector<AVAR> > pMatrix(1);
	pMatrix[0].resize(11);
	pMatrix[0][0] = 5; // CL
	pMatrix[0][1] = 8; // Q
	pMatrix[0][2] = 20; // Vc
	pMatrix[0][3] = 70; // Vp
	pMatrix[0][4] = 1.2; // ka
	pMatrix[0][5] = 1; // F1
	pMatrix[0][6] = 1; // F2
	pMatrix[0][7] = 1; // F3
	pMatrix[0][8] = 0; // tlag1
	pMatrix[0][9] = 0; // tlag2
	pMatrix[0][10] = 0; // tlag3

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

	stan::math::matrix_v x;
	x = PKModelTwoCpt(pMatrix, time, amt, rate, ii, evid, cmt, addl, ss);

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

	for (int i = 0; i < amounts.rows(); i++) {
		EXPECT_NEAR(amounts(i, 0), x(i, 0).val(),
		  std::max(amounts(i, 0), x(i, 0).val()) * 1e-6);
		EXPECT_NEAR(amounts(i, 1), x(i, 1).val(),
		  std::max(amounts(i, 1), x(i, 1).val()) * 1e-6);
	}
}

TEST(Torsten, PKModelTwoCpt_MultipleDoses_timePara) {

    int nEvent = 11;
    
	vector<vector<AVAR> > pMatrix(nEvent);
	
	for (int i = 0; i < nEvent; i++) {
	  pMatrix[i].resize(11);
	  if (i < 6) pMatrix[i][0] = 5; // CL
	  else pMatrix[i][0] = 50; // CL is piece-wise constant
	  pMatrix[i][1] = 8; // Q
	  pMatrix[i][2] = 20; // Vc
	  pMatrix[i][3] = 70; // Vp
	  pMatrix[i][4] = 1.2; // ka
	  pMatrix[i][5] = 1; // F1
	  pMatrix[i][6] = 1; // F2
	  pMatrix[i][7] = 1; // F3
	  pMatrix[i][8] = 0; // tlag1
	  pMatrix[i][9] = 0; // tlag2
	  pMatrix[i][10] = 0; // tlag3
	}

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

	stan::math::matrix_v x;
	x = PKModelTwoCpt(pMatrix, time, amt, rate, ii, evid, cmt, addl, ss);

	Matrix<double, Dynamic, Dynamic> amounts(nEvent, 3);
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
			   
    for (int i = 0; i < amounts.rows(); i++) {
		EXPECT_NEAR(amounts(i, 0), x(i, 0).val(),
		  std::max(amounts(i, 0), x(i, 0).val()) * 1e-6);
		EXPECT_NEAR(amounts(i, 1), x(i, 1).val(),
		  std::max(amounts(i, 1), x(i, 1).val()) * 1e-6);
	}
}
