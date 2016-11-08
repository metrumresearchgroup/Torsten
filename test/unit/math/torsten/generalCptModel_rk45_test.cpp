#include <stan/math/rev/mat.hpp>
#include <gtest/gtest.h>
#include <test/unit/math/prim/mat/fun/expect_matrix_eq.hpp>
#include <stan/math/torsten/GeneralCptModel_rk45.hpp>


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
  vector<scalar> y(2, 0);
  
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

TEST(Torsten, genCptrk45_One_SingleDose) {
  using std::vector;
  
  vector<Matrix<double, Eigen::Dynamic, 1> > pMatrix(1);
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

  Matrix<double, Eigen::Dynamic, Eigen::Dynamic> x;
  x = generalCptModel_rk45(oneCptModelODE_functor(), 2,
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
			  
  expect_matrix_eq(amounts, x);
}
