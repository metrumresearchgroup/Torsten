#include <stan/math/rev/mat.hpp>  // FIX ME - includes should be more specific
#include <gtest/gtest.h>
#include <test/unit/math/prim/mat/fun/expect_near_matrix_eq.hpp>
#include <test/unit/math/prim/mat/fun/expect_matrix_eq.hpp>
#include <test/unit/math/torsten/prim/util_generalOdeModel.hpp>

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

TEST(Torsten, genCpt_One_SingleDose) {
  using std::vector;
  using Eigen::Matrix;
  using Eigen::Dynamic;

  double rel_err = 1e-6;

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

  double rel_tol = 1e-8, abs_tol = 1e-8;
  long int max_num_steps = 1e8;
  Matrix<double, Eigen::Dynamic, Eigen::Dynamic> x_rk45;
  x_rk45 = generalCptModel_rk45(oneCptModelODE_functor(), nCmt,
                                time, amt, rate, ii, evid, cmt, addl, ss,
                                pMatrix, biovar, tlag,
                                rel_tol, abs_tol, max_num_steps);

  Matrix<double, Eigen::Dynamic, Eigen::Dynamic> x_bdf;
  x_bdf = generalCptModel_bdf(oneCptModelODE_functor(), nCmt,
                              time, amt, rate, ii, evid, cmt, addl, ss,
                              pMatrix, biovar, tlag,
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
  double diff = 1e-8, diff2 = 3e-3;
  test_generalOdeModel(oneCptModelODE_functor(), nCmt,
                       time, amt, rate, ii, evid, cmt, addl, ss,
                       pMatrix, biovar, tlag,
                       rel_tol, abs_tol, max_num_steps, diff, diff2, "rk45");
  test_generalOdeModel(oneCptModelODE_functor(), nCmt,
                       time, amt, rate, ii, evid, cmt, addl, ss,
                       pMatrix, biovar, tlag,
                       rel_tol, abs_tol, max_num_steps, diff, diff2, "bdf");
  
}

TEST(Torsten, genCpt_One_SingleDose_overload) {
  using std::vector;
  using Eigen::Matrix;
  using Eigen::Dynamic;

  double rel_err = 1e-6;
  
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

  double rel_tol = 1e-8, abs_tol = 1e-8;
  long int max_num_steps = 1e8;
  
  Matrix<double, Eigen::Dynamic, Eigen::Dynamic> x_rk45_122, x_rk45_112, 
    x_rk45_111, x_rk45_121, x_rk45_212, x_rk45_211, x_rk45_221;
  x_rk45_122 = generalCptModel_rk45(oneCptModelODE_functor(), nCmt,
                                time, amt, rate, ii, evid, cmt, addl, ss,
                                pMatrix[0], biovar, tlag,
                                rel_tol, abs_tol, max_num_steps);

  x_rk45_112 = generalCptModel_rk45(oneCptModelODE_functor(), nCmt,
                                    time, amt, rate, ii, evid, cmt, addl, ss,
                                    pMatrix[0], biovar[0], tlag,
                                    rel_tol, abs_tol, max_num_steps);
  
  x_rk45_111 = generalCptModel_rk45(oneCptModelODE_functor(), nCmt,
                                    time, amt, rate, ii, evid, cmt, addl, ss,
                                    pMatrix[0], biovar[0], tlag[0],
                                    rel_tol, abs_tol, max_num_steps);
  
  x_rk45_121 = generalCptModel_rk45(oneCptModelODE_functor(), nCmt,
                                    time, amt, rate, ii, evid, cmt, addl, ss,
                                    pMatrix[0], biovar, tlag[0],
                                    rel_tol, abs_tol, max_num_steps);
  
  x_rk45_212 = generalCptModel_rk45(oneCptModelODE_functor(), nCmt,
                                    time, amt, rate, ii, evid, cmt, addl, ss,
                                    pMatrix, biovar[0], tlag,
                                    rel_tol, abs_tol, max_num_steps);
  
  x_rk45_211 = generalCptModel_rk45(oneCptModelODE_functor(), nCmt,
                                    time, amt, rate, ii, evid, cmt, addl, ss,
                                    pMatrix, biovar[0], tlag[0],
                                    rel_tol, abs_tol, max_num_steps);
  
  x_rk45_221 = generalCptModel_rk45(oneCptModelODE_functor(), nCmt,
                                    time, amt, rate, ii, evid, cmt, addl, ss,
                                    pMatrix, biovar, tlag[0],
                                    rel_tol, abs_tol, max_num_steps);


  Matrix<double, Eigen::Dynamic, Eigen::Dynamic> x_bdf_122, x_bdf_112, 
    x_bdf_111, x_bdf_121, x_bdf_212, x_bdf_211, x_bdf_221;
  x_bdf_122 = generalCptModel_bdf(oneCptModelODE_functor(), nCmt,
                                    time, amt, rate, ii, evid, cmt, addl, ss,
                                    pMatrix[0], biovar, tlag,
                                    rel_tol, abs_tol, max_num_steps);
  
  x_bdf_112 = generalCptModel_bdf(oneCptModelODE_functor(), nCmt,
                                    time, amt, rate, ii, evid, cmt, addl, ss,
                                    pMatrix[0], biovar[0], tlag,
                                    rel_tol, abs_tol, max_num_steps);
  
  x_bdf_111 = generalCptModel_bdf(oneCptModelODE_functor(), nCmt,
                                    time, amt, rate, ii, evid, cmt, addl, ss,
                                    pMatrix[0], biovar[0], tlag[0],
                                    rel_tol, abs_tol, max_num_steps);
  
  x_bdf_121 = generalCptModel_bdf(oneCptModelODE_functor(), nCmt,
                                    time, amt, rate, ii, evid, cmt, addl, ss,
                                    pMatrix[0], biovar, tlag[0],
                                    rel_tol, abs_tol, max_num_steps);
  
  x_bdf_212 = generalCptModel_bdf(oneCptModelODE_functor(), nCmt,
                                    time, amt, rate, ii, evid, cmt, addl, ss,
                                    pMatrix, biovar[0], tlag,
                                    rel_tol, abs_tol, max_num_steps);
  
  x_bdf_211 = generalCptModel_bdf(oneCptModelODE_functor(), nCmt,
                                    time, amt, rate, ii, evid, cmt, addl, ss,
                                    pMatrix, biovar[0], tlag[0],
                                    rel_tol, abs_tol, max_num_steps);
  
  x_bdf_221 = generalCptModel_bdf(oneCptModelODE_functor(), nCmt,
                                    time, amt, rate, ii, evid, cmt, addl, ss,
                                    pMatrix, biovar, tlag[0],
                                    rel_tol, abs_tol, max_num_steps);  
	
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

  expect_near_matrix_eq(amounts, x_rk45_122, rel_err);
  expect_near_matrix_eq(amounts, x_rk45_112, rel_err);
  expect_near_matrix_eq(amounts, x_rk45_111, rel_err);
  expect_near_matrix_eq(amounts, x_rk45_121, rel_err);
  expect_near_matrix_eq(amounts, x_rk45_212, rel_err);
  expect_near_matrix_eq(amounts, x_rk45_211, rel_err);
  expect_near_matrix_eq(amounts, x_rk45_221, rel_err);

  expect_near_matrix_eq(amounts, x_bdf_122, rel_err);
  expect_near_matrix_eq(amounts, x_bdf_112, rel_err);
  expect_near_matrix_eq(amounts, x_bdf_111, rel_err);
  expect_near_matrix_eq(amounts, x_bdf_121, rel_err);
  expect_near_matrix_eq(amounts, x_bdf_212, rel_err);
  expect_near_matrix_eq(amounts, x_bdf_211, rel_err);
  expect_near_matrix_eq(amounts, x_bdf_221, rel_err);
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

  vector<vector<double> > pMatrix(1);
  pMatrix[0].resize(5);
  pMatrix[0][0] = 10; // CL0
  pMatrix[0][1] = 80; // Vc
  pMatrix[0][2] = 1.2; // ka
  pMatrix[0][3] = 2; // CLSS
  pMatrix[0][4] = 1; // K
  
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

  double rel_tol = 1e-8, abs_tol = 1e-8;
  long int max_num_steps = 1e8;

  Matrix<double, Eigen::Dynamic, Eigen::Dynamic> x_rk45;
  x_rk45 = generalCptModel_rk45(oneCptModelODE_abstime_functor(), 2,
                                time, amt, rate, ii, evid, cmt, addl, ss,
                                pMatrix, biovar, tlag,
                                rel_tol, abs_tol, max_num_steps);
  
  Matrix<double, Eigen::Dynamic, Eigen::Dynamic> x_bdf;
  x_bdf = generalCptModel_bdf(oneCptModelODE_abstime_functor(), 2,
                              time, amt, rate, ii, evid, cmt, addl, ss,
                              pMatrix, biovar, tlag,
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
   x_rk45 = generalCptModel_rk45(oneCptModelODE_functor(), nCmt,
                                 time, amt, rate, ii, evid, cmt, addl, ss,
                                 pMatrix, biovar, tlag,
                                 rel_tol, abs_tol, max_num_steps);
  
  Matrix<double, Eigen::Dynamic, Eigen::Dynamic> x_bdf;
  x_bdf = generalCptModel_bdf(oneCptModelODE_functor(), nCmt,
                              time, amt, rate, ii, evid, cmt, addl, ss,
                              pMatrix, biovar, tlag,
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
  double diff = 1e-8, diff2 = 1e-2;
  test_generalOdeModel(oneCptModelODE_functor(), nCmt,
                       time, amt, rate, ii, evid, cmt, addl, ss,
                       pMatrix, biovar, tlag,
                       rel_tol, abs_tol, max_num_steps, diff, diff2, "rk45");
  test_generalOdeModel(oneCptModelODE_functor(), nCmt,
                       time, amt, rate, ii, evid, cmt, addl, ss,
                       pMatrix, biovar, tlag,
                       rel_tol, abs_tol, max_num_steps, diff, diff2, "bdf");
}
