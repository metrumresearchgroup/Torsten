#include <stan/math/rev/mat.hpp>  // FIX ME - more specific
#include <stan/math/torsten/PKModel/Pred/Pred1_mix1.hpp>
#include <stan/math/torsten/PKModel/functors/mix1_functor.hpp>
#include <gtest/gtest.h>

struct ODE_functor {
  template <typename T0, typename T1, typename T2, typename T3>
  inline
  std::vector<typename boost::math::tools::promote_args<T0, T1, T2, T3>::type>
  operator()(const T0& t,
             const std::vector<T1>& x,
             const std::vector<T2>& x_pk,
             const std::vector<T3>& theta,
             const std::vector<double>& x_r,
             const std::vector<int>& x_i,
             std::ostream* pstream_) const {
    typedef typename boost::math::tools::promote_args<T0, T1, T2, T3>::type
      scalar;

    scalar VC = theta[1],
      Mtt = theta[3],
      circ0 = theta[4],
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
};


TEST(Torsten, pred1_mix) {
  using Eigen::Matrix;
  using Eigen::Dynamic;

  double dt = 1;

  int nParameters = 7;
  std::vector<double> parameters(nParameters);
  parameters[0] = 10;  // CL
  parameters[1] = 35;  // VC
  parameters[2] = 2.0;  // ka
  parameters[3] = 125;  // Mtt
  parameters[4] = 5;  // Circ0
  parameters[5] = 3e-4;  // alpha
  parameters[6] = 0.17;  // gamma

  std::vector<double> biovar_dummy(1, 0);
  std::vector<double> tlag_dummy(1, 0);
  Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> K(0, 0);

  // initialize Model Parameters object
  torsten::ModelParameters<double, double, double, double>
    parms(dt, parameters, biovar_dummy, tlag_dummy);

  int nOdes = 5;
  Matrix<double, 1, Dynamic> init(nOdes);
  init << 80000, 0, 0, 0, 0;  // initial dose in the gut

  std::vector<double> rate(5, 0);  // no rate

  // Construct pmetrics solver
  // Tune parameters to get better result (even if this makes the operation
  // more costly). FIX ME - is this the right approach?
  double rel_tol = 1e-8, abs_tol = 1e-8;
  long int max_num_steps = 1e+8;

  typedef torsten::mix1_functor<ODE_functor> F0;
  torsten::Pred1_mix1<F0> Pred1(F0(ODE_functor()), rel_tol, abs_tol, max_num_steps, 0,
                      "rk45");
  Matrix<double, 1, Dynamic> pred = Pred1(dt, parms, init, rate);

  // Compare to results obtained with mrgsolve
  Eigen::VectorXd mrgResults(5);
  mrgResults << 10826.82, 57506.59, -0.0556872,
                -0.000694898, -6.103923e-6;

  EXPECT_FLOAT_EQ(mrgResults(0), pred(0));
  EXPECT_FLOAT_EQ(mrgResults(1), pred(1));
  EXPECT_FLOAT_EQ(mrgResults(2), pred(2));
  EXPECT_FLOAT_EQ(mrgResults(3), pred(3));
  EXPECT_NEAR(mrgResults(4), pred(4), std::abs(mrgResults(4) * 1e-4));
  // NOTE - for the last value, get very low value for pred(4)
  // which may introduce double precision errors. Still, the two
  // results agree within a 1e-4 relative error.
}
