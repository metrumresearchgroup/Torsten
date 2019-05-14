#include <stan/math/rev/mat.hpp>  // FIX ME - more specific
#include <stan/math/torsten/PKModel/Pred/Pred1_mix2.hpp>
#include <stan/math/torsten/PKModel/functors/mix2_functor.hpp>
#include <gtest/gtest.h>

struct ODE_functor {
  template <typename T0, typename T1, typename T2, typename T3>
  inline
  std::vector<typename boost::math::tools::promote_args<T0, T1, T2, T3>::type>
  operator()(const T0& t,
             const std::vector<T1>& y,
             const std::vector<T2>& y_pk,
             const std::vector<T3>& theta,
             const std::vector<double>& x_r,
             const std::vector<int>& x_i,
             std::ostream* pstream_) const {
    typedef typename boost::math::tools::promote_args<T0, T1, T2, T3>::type
      scalar;

    scalar VC = theta[2],
      Mtt = theta[5],
      circ0 = theta[6],
      alpha = theta[7],
      gamma = theta[8],
      ktr = 4 / Mtt,
      prol = y[0] + circ0,
      transit = y[1] + circ0,
      circ = y[2] + circ0,
      conc = y_pk[1] / VC,
      Edrug = alpha * conc;

    std::vector<scalar> dxdt(3);
    dxdt[0] = ktr * prol * ((1 - Edrug) * pow(circ0 / circ, gamma) - 1);
    dxdt[1] = ktr * (prol - transit);
    dxdt[2] = ktr * (transit - circ);

    return dxdt;
  }
};


TEST(Torsten, pred1_mix) {
  using Eigen::Matrix;
  using Eigen::Dynamic;

  double dt = 1.0;

  int nParameters = 9;
  std::vector<double> parameters(nParameters);
  parameters[0] = 10;  // CL
  parameters[1] = 15;  // Q
  parameters[2] = 35;  // VC
  parameters[3] = 105;  // VP
  parameters[4] = 2.0;  // ka
  parameters[5] = 125;  // Mtt
  parameters[6] = 5;  // Circ0
  parameters[7] = 3e-4;  // alpha
  parameters[8] = 0.17;  // gamma

  std::vector<double> biovar_dummy(1, 0);
  std::vector<double> tlag_dummy(1, 0);
  Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> K(0, 0);

  // initialize Model Parameters object
  torsten::ModelParameters<double, double, double, double>
    parms(dt, parameters, biovar_dummy, tlag_dummy);

  int nOdes = 6;
  Matrix<double, 1, Dynamic> init(nOdes);
  init << 10000, 0, 0, 0, 0, 0;  // initial dose in the gut

  std::vector<double> rate(6, 0);  // no rate

  // Use default value of parameters
  double rel_tol = 1e-6, abs_tol = 1e-6;
  long int max_num_steps = 1e+6;

  typedef torsten::mix2_functor<ODE_functor> F0;
  torsten::Pred1_mix2<F0> Pred1(F0(ODE_functor()), rel_tol, abs_tol, max_num_steps, 0,
                       "rk45");
  Matrix<double, 1, Dynamic> pred = Pred1(dt, parms, init, rate);

  // Compare to results obtained with mrgsolve
  Eigen::VectorXd mrgResults(6);
  mrgResults << 1353.352829, 5597.489, 1787.0134,
                -0.0060546255, -7.847821e-05, -7.039447e-07;

  // PK compartments
  EXPECT_FLOAT_EQ(mrgResults(0), pred(0));
  EXPECT_FLOAT_EQ(mrgResults(1), pred(1));
  EXPECT_FLOAT_EQ(mrgResults(2), pred(2));

  // PD compartments
  // (are estimated less accurately)
  EXPECT_NEAR(mrgResults(3), pred(3), std::abs(mrgResults(3) * 1e-4));
  EXPECT_NEAR(mrgResults(4), pred(4), std::abs(mrgResults(4) * 1e-4));
  EXPECT_NEAR(mrgResults(5), pred(5), std::abs(mrgResults(5) * 1e-4));
}
