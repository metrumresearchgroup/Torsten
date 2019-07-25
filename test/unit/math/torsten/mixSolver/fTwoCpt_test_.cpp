#include <stan/math/rev/mat.hpp>  // FIX ME - more specific
#include <gtest/gtest.h>
#include <test/unit/math/torsten/expect_near_matrix_eq.hpp>
#include <test/unit/math/torsten/expect_matrix_eq.hpp>

TEST(Torsten, fTwoCpt) {
  double dt = 0.25;

  int nParameters = 5;
  Eigen::VectorXd parameters(nParameters);
  parameters << 5, 8, 20, 70, 1.2;  // CL, Q, VC, VP, ka

  int nCmt = 3;
  Eigen::VectorXd init(nCmt);
  init << 1000, 0, 0;  // initial dose in the gut

  std::vector<double> rate(3, 0);  // no rate

  Eigen::VectorXd pred;
  pred = torsten::fTwoCpt(dt, parameters, init, rate);

  EXPECT_FLOAT_EQ(740.8182, pred(0));
  EXPECT_FLOAT_EQ(238.3713, pred(1));
  EXPECT_FLOAT_EQ(12.75775, pred(2));
}


struct feedbackODE {
  template <typename T0, typename T1, typename T2, typename T3>
  inline
  std::vector<typename boost::math::tools::promote_args<T0, T1, T2, T3>::type>
  operator()(const T0& t,
             const std::vector<T1>& x,
             const std::vector<T2>& parms,
             const std::vector<T3>& rate,
             const std::vector<int>& dummy, std::ostream* pstream__) const {
    typedef typename boost::math::tools::promote_args<T0, T1, T2, T3>::type scalar;
    // PK variables
    scalar CL = parms[0], Q = parms[1], VC = parms[2], VP = parms[3],
      ka = parms[4], k10 = CL / VC, k12 = Q / VC,  k21 = Q / VP;

    // PD variables
    scalar MTT = parms[5], circ0 = parms[6], alpha = parms[7], gamma = parms[8],
       ktr = 4 / MTT, prol = x[3] + circ0, transit1 = x[4] + circ0,
       transit2 = x[5] + circ0, transit3 = x[6] + circ0,
       circ = x[7] + circ0;

    // return object
    std::vector<scalar> dxdt(8);

    dxdt[0] = -ka * x[0];
    dxdt[1] = ka * x[0] - (k10 + k12) * x[1] + k21 * x[2];
    dxdt[2] = k12 * x[1] - k21 * x[2];
    scalar conc = x[1] / VC;
    scalar Edrug = alpha * conc;

    dxdt[3] = ktr * prol * ((1 - Edrug) * (pow((circ0 / circ), gamma)) - 1);
    dxdt[4] = ktr * (prol - transit1);
    dxdt[5] = ktr * (transit1 - transit2);
    dxdt[6] = ktr * (transit2 - transit3);
    dxdt[7] = ktr * (transit3 - circ);

    return dxdt;
   }
};

struct feedbackODE_mixed {
  template <typename T0, typename T1, typename T2, typename T3>
  inline
  std::vector<typename boost::math::tools::promote_args<T0, T1, T2, T3>::type>
  operator()(const T0& t,
             const std::vector<T1>& x,
             const std::vector<T2>& parms,
             const std::vector<T3>& rate,
             const std::vector<int>& dummy, std::ostream* pstream__) const {
    using Eigen::Matrix;
    using Eigen::Dynamic;
    typedef typename boost::math::tools::promote_args<T0, T1, T2, T3>::type scalar;
    // PK variables
    scalar VC = parms[2];

    // PD variables
    scalar MTT = parms[5], circ0 = parms[6], alpha = parms[7], gamma = parms[8],
       ktr = 4 / MTT, prol = x[0] + circ0, transit1 = x[1] + circ0,
       transit2 = x[2] + circ0, transit3 = x[3] + circ0,
       circ = x[4] + circ0;

    Matrix<scalar, Dynamic, 1> initPK(3);
    initPK(0) = parms[9];
    initPK(1) = parms[10];
    initPK(2) = parms[11];

    Matrix<scalar, Dynamic, 1> parmsPK(5);
    for (int i = 0; i < parmsPK.size(); i++)
      parmsPK(i) = parms[i];


    Matrix<scalar, Dynamic, 1> predPK = torsten::fTwoCpt(t, parmsPK, initPK, rate);

    scalar conc = predPK(1) / VC;
    scalar Edrug = alpha * conc;

    // return object
    std::vector<scalar> dxdt(5);
    dxdt[0] = ktr * prol * ((1 - Edrug) * (pow((circ0 / circ), gamma)) - 1);
    dxdt[1] = ktr * (prol - transit1);
    dxdt[2] = ktr * (transit1 - transit2);
    dxdt[3] = ktr * (transit2 - transit3);
    dxdt[4] = ktr * (transit3 - circ);

    return dxdt;
   }
};

TEST(Torsten, ODEs) {
  using std::vector;
  using stan::math::var;
  using stan::math::integrate_ode_rk45;

  vector<var> init(8, 0);
  init[0] = 1000;
  double t0 = 0;
  vector<double> t(3);
  t[0] = 1;
  t[1] = 2;
  t[2] = 10;

  vector<var> parms(9);
  parms[0] = 10;
  parms[1] = 15;
  parms[2] = 35;
  parms[3] = 105;
  parms[4] = 2;
  parms[5] = 125;
  parms[6] = 5;
  parms[7] = 3e-4;
  parms[8] = 0.17;

  vector<double> rate(8, 0);

  vector<double> rdummy(1, 0);
  vector<int> idummy(1, 0);

  vector<vector<var> > temp = integrate_ode_rk45(feedbackODE(), init, t0, t, parms, rdummy, idummy);

  /* Mixed Solver */
  // Add init in the parameters
  parms.push_back(init[0]);
  parms.push_back(init[1]);
  parms.push_back(init[2]);

  vector<var> init_mixed(5, 0);
  vector<vector<var> > temp_mixed = integrate_ode_rk45(feedbackODE_mixed(), init_mixed,
                                                       t0, t, parms, rate, idummy);

  for (size_t i = 0; i < temp_mixed.size(); i++)
    for (size_t j = 0; j < temp_mixed[i].size(); j++)
      EXPECT_NEAR(temp[i][j + 3].val(), temp_mixed[i][j].val(), 1e-8);
}
