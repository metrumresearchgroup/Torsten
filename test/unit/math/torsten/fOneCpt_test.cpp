#include <stan/math/rev/mat.hpp>  // FIX ME - more specific
#include <gtest/gtest.h>
#include <test/unit/math/prim/mat/fun/expect_near_matrix_eq.hpp>
#include <test/unit/math/prim/mat/fun/expect_matrix_eq.hpp>
#include <test/unit/math/torsten/util_generalOdeModel.hpp>

TEST(Torsten, fOneCpt) {
  double dt = 0.25;

  int nParameters = 3;
  Eigen::VectorXd parameters(nParameters);
  parameters << 10, 80, 1.2;  // CL, V2, ka

  int nCmt = 2;
  Eigen::VectorXd init(nCmt);
  init << 1000, 0;  // initial dose in the gut

  std::vector<double> rate(2, 0);  // no rate

  Eigen::VectorXd pred;
  pred = fOneCpt(dt, parameters, init, rate);

  EXPECT_FLOAT_EQ(740.8182, pred(0));
  EXPECT_FLOAT_EQ(254.97490, pred(1));
}

/*
struct feedback_functor {
  template <typename T0, typename T1, typename T2, typename T3>
  inline
  std::vector<typename boost::math::tools::promote_args<T0, T1, T2, T3>::type>
  operator()(const T0& t,
             const std::vector<T1>& x,
             const std::vector<T2>& parms,
             const std::vector<T3>& rate,
             const std::vector<int>& dummy,
             std::ostream* pstream__) const {
  typedef
    typename boost::math::tools::promote_args<T0, T1, T2, T3>::type scalar;

  scalar VC = parms[2],
         MTT = parms[4],
         circ0 = parms[5],
         alpha = parms[6],
         gamma = parms[7],
         ktr = 4 / MTT,
         prol = x[3] + circ0,
         transit = x[4] + circ0,
         circ = x[5] + circ0;
  
  Eigen::Matrix<scalar, 1, Eigen::Dynamic> pred1, init(2);
  init(0) = parms[8];
  init(1) = parms[9];
  pred1 = fOneCpt(t, parameter, init, rate);
  
  scalar conc = pred(1) / VC;
  scalar Edrug = alpha * conc;

  std::vector<scalar> dxdt(3);
  dxdt[0] = ktr * prol * ((1 - Edrug) * ((circ0 / circ) ^ gamma) - 1);
  dxdt[1] = ktr * (prol - transit);
  dxdt[2] = ktr * (transit - circ);
  
  return dxdt;
}


TEST(Torsten, mixed_solver) {
  // write test (maybe required later)
} */

