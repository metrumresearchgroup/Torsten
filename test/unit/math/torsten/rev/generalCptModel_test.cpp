#include <stan/math/rev/mat.hpp>
#include <gtest/gtest.h>
#include <test/unit/math/rev/mat/fun/util.hpp>

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

TEST(Torsten, genCpt_One_SingleDose) {
  using std::vector;
  using Eigen::Matrix;
  using Eigen::Dynamic;

  double rel_err = 1e-6;
  
  vector<vector<AVAR> > pMatrix(1);
  pMatrix[0].resize(7);
  pMatrix[0][0] = 10; // CL
  pMatrix[0][1] = 80; // Vc
  pMatrix[0][2] = 1.2; // ka
  pMatrix[0][3] = 1; // F1
  pMatrix[0][4] = 1; // F2
  pMatrix[0][5] = 0; // tlag1
  pMatrix[0][6] = 0; // tlag2

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

  stan::math::matrix_v x_rk45;
  x_rk45 = generalCptModel_rk45(oneCptModelODE_functor(), 2,
                                pMatrix, time, amt, rate, ii, evid, cmt, addl, ss,
                                1e-8, 1e-8, 1e8);

  stan::math::matrix_v x_bdf;
  x_bdf = generalCptModel_bdf(oneCptModelODE_functor(), 2,
                              pMatrix, time, amt, rate, ii, evid, cmt, addl, ss,
                              1e-8, 1e-8, 1e8);
	
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

  for (int i = 0; i < amounts.rows(); i++)
    for (int j = 0; j < amounts.cols(); j++) {
      EXPECT_NEAR(amounts(i, j), x_rk45(i, j).val(),
        std::max(amounts(i, j), x_rk45(i, j).val()) * rel_err);
      EXPECT_NEAR(amounts(i, j), x_bdf(i, j).val(),
        std::max(amounts(i, j), x_bdf(i, j).val()) * rel_err);
  }
}

TEST(Torsten, genCpt_One_SingleDose_overload) {
  using std::vector;
  using Eigen::Matrix;
  using Eigen::Dynamic;

  double rel_err = 1e-6;
  
  vector<AVAR> pMatrix(7);
  pMatrix[0] = 10; // CL
  pMatrix[1] = 80; // Vc
  pMatrix[2] = 1.2; // ka
  pMatrix[3] = 1; // F1
  pMatrix[4] = 1; // F2
  pMatrix[5] = 0; // tlag1
  pMatrix[6] = 0; // tlag2

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

  stan::math::matrix_v x_rk45;
  x_rk45 = generalCptModel_rk45(oneCptModelODE_functor(), 2,
                                pMatrix, time, amt, rate, ii, evid, cmt, addl, ss,
                                1e-8, 1e-8, 1e8);

  stan::math::matrix_v x_bdf;
  x_bdf = generalCptModel_bdf(oneCptModelODE_functor(), 2,
                              pMatrix, time, amt, rate, ii, evid, cmt, addl, ss,
                              1e-8, 1e-8, 1e8);
	
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

  for (int i = 0; i < amounts.rows(); i++)
    for (int j = 0; j < amounts.cols(); j++) {
      EXPECT_NEAR(amounts(i, j), x_rk45(i, j).val(),
        std::max(amounts(i, j), x_rk45(i, j).val()) * rel_err);
      EXPECT_NEAR(amounts(i, j), x_bdf(i, j).val(),
        std::max(amounts(i, j), x_bdf(i, j).val()) * rel_err);
  }
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
  using Eigen::Dynamic;

  double rel_err = 1e-6;
  
  vector<vector<AVAR> > pMatrix(1);
  pMatrix[0].resize(9);
  pMatrix[0][0] = 10; // CL0
  pMatrix[0][1] = 80; // Vc
  pMatrix[0][2] = 1.2; // ka
  pMatrix[0][3] = 2; // CLSS
  pMatrix[0][4] = 1; // K 
  pMatrix[0][5] = 1; // F1
  pMatrix[0][6] = 1; // F2
  pMatrix[0][7] = 0; // tlag1
  pMatrix[0][8] = 0; // tlag2

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

  stan::math::matrix_v x_rk45;
  x_rk45 = generalCptModel_rk45(oneCptModelODE_abstime_functor(), 2,
                                pMatrix, time, amt, rate, ii, evid, cmt, addl, ss,
                                1e-8, 1e-8, 1e8);

  stan::math::matrix_v x_bdf;
  x_bdf = generalCptModel_bdf(oneCptModelODE_abstime_functor(), 2,
                              pMatrix, time, amt, rate, ii, evid, cmt, addl, ss,
                              1e-8, 1e-8, 1e8);

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
			  
  for (int i = 0; i < amounts.rows(); i++)
    for (int j = 0; j < amounts.cols(); j++) {
      EXPECT_NEAR(amounts(i, j), x_rk45(i, j).val(),
        std::max(amounts(i, j), x_rk45(i, j).val()) * rel_err);
      EXPECT_NEAR(amounts(i, j), x_bdf(i, j).val(),
        std::max(amounts(i, j), x_bdf(i, j).val()) * rel_err);
  }
}

TEST(Torsten, genCpOne_MultipleDoses_timePara) {
    double rel_err_rk45 = 1e-6;
    double rel_err_bdf = 1e-4;
    using std::vector;
    using Eigen::Matrix;
    using Eigen::Dynamic;

    int nEvent = 11;
	vector<vector<AVAR> > pMatrix(nEvent);
	
	for (int i = 0; i < nEvent; i++) {
	  pMatrix[i].resize(7);
	  if (i < 6) pMatrix[i][0] = 10; // CL
	  else pMatrix[i][0] = 50; // CL is piece-wise constant
	  pMatrix[i][1] = 80; // Vc
	  pMatrix[i][2] = 1.2; // ka
	  pMatrix[i][3] = 1; // F1
	  pMatrix[i][4] = 1; // F2
	  pMatrix[i][5] = 0; // tlag1
	  pMatrix[i][6] = 0; // tlag2
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

    stan::math::matrix_v x_rk45;
    x_rk45 = generalCptModel_rk45(oneCptModelODE_functor(), 2,
                                  pMatrix, time, amt, rate, ii, evid, cmt, addl, ss,
                                  1e-8, 1e-8, 1e8);
  
    stan::math::matrix_v x_bdf;
    x_bdf = generalCptModel_bdf(oneCptModelODE_functor(), 2,
                                pMatrix, time, amt, rate, ii, evid, cmt, addl, ss,
                                1e-8, 1e-8, 1e8);

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
			   
  for (int i = 0; i < amounts.rows(); i++)
    for (int j = 0; j < amounts.cols(); j++) {
      EXPECT_NEAR(amounts(i, j), x_rk45(i, j).val(),
        std::max(amounts(i, j), x_rk45(i, j).val()) * rel_err_rk45);
      EXPECT_NEAR(amounts(i, j), x_bdf(i, j).val(),
        std::max(amounts(i, j), x_bdf(i, j).val()) * rel_err_bdf);
  }
}
