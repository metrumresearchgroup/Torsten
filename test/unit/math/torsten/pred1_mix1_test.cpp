#include <stan/math/rev/mat.hpp>  // FIX ME - more specific
#include <gtest/gtest.h>
#include <test/unit/math/prim/mat/fun/expect_near_matrix_eq.hpp>
#include <test/unit/math/prim/mat/fun/expect_matrix_eq.hpp>
#include <test/unit/math/torsten/util_generalOdeModel.hpp>

struct ODE_functor {

  template <typename T0, typename T1, typename T2, typename T3>
  inline
  std::vector<typename boost::math::tools::promote_args<T0, T1, T2, T3>::type>
  operator()(const T0& t,
             const std::vector<T1>& x,
             const std::vector<T2>& x_pk;
             const std::vector<T3>& theta,
             const std::vector<double>& x_r,
             const std::vector<int>& x_i,
             std::ostream* pstream_) const {
    typedef typename boost::math::tools::promote_args<T0, T1, T2, T3>::type
      scalar;

    scalar VC = theta[1],
      Mtt = theta[3],
      Circ0 = theta[4],
      alpha = theta[5],
      gamma = theta[6],
      ktr = 4 / Mtt,
      prol = x[0] + circ0,
      transit = x[1] + circ0,
      circ = x[2] + circ0,
      conc = x_pk[1] / VC,
      Edrug = alpha * conc;

    std::vector<scalar> dxdt(3);
    dxdt[0] = ktr * prol * ((1 - Edrug) * pow(circ0 / circ, gamma) - 1);
    dxdt[1] = ktr * (prol - transit);
    dxdt[2] = ktr * (transit - circ);

    return dxdt;
  }

}


TEST(Torsten, pred1_mix) {
  double dt = 1;

  int nParameters = 3;
  Eigen::VectorXd parameters(nParameters);
  // CL, VC, ka, Mtt, Circ0, alpha, gamma
  parameters << 10, 35, 2.0, 125, 5, 3e-4, 0.17;

  // initialize Model Parameters object

  int nCmt = 2;
  Eigen::VectorXd init(nCmt);
  init << 80000, 0;  // initial dose in the gut

  std::vector<double> rate(2, 0);  // no rate

  Eigen::VectorXd pred;
  pred = Pred1_mix1(dt, parameters, init, rate);

  EXPECT_FLOAT_EQ(740.8182, pred(0));
  EXPECT_FLOAT_EQ(254.97490, pred(1));
}
