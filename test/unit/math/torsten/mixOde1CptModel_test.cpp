#include <stan/math/rev/mat.hpp>  // FIX ME - includes should be more specific
#include <gtest/gtest.h>
#include <test/unit/math/prim/mat/fun/expect_near_matrix_eq.hpp>
#include <test/unit/math/prim/mat/fun/expect_matrix_eq.hpp>

struct feebackODE {
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

    scalar
      VC = theta[1],
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

TEST(Torsten, mixOde1Cpt) {
  using std::vector;
  using Eigen::Matrix;
  using Eigen::Dynamic;

  int nEvents = 10;
  vector<double> time(nEvents);
  time[0] = 0.0;
  for (int i = 1; i < (nEvents - 1); i++) time[i] = time[i - 1] + 0.25;
  time[9] = 4.0;
  vector<double> amt(nEvents, 0);
  amt[0] = 10000;
  vector<double> rate(nEvents, 0);
  vector<int> cmt(nEvents, 2);
  cmt[0] = 1;
  vector<int> evid(nEvents, 0);
  evid[0] = 1;
  vector<double> ii(nEvents, 0);
  vector<int> addl(nEvents, 0);
  vector<int> ss(nEvents, 0);

  int nParameters = 7;
  vector<vector<double> > parameters(1);
  parameters[0].resize(nParameters);
  parameters[0][0] = 10;  // CL
  parameters[0][1] = 35;  // VC
  parameters[0][2] = 2.0;  // ka
  parameters[0][3] = 125;  // Mtt
  parameters[0][4] = 5;  // Circ0
  parameters[0][5] = 3e-4;  // alpha
  parameters[0][6] = 0.17;  // gamma

  int nOde = 5;
  vector<vector<double> > biovar(1);
  biovar[0].resize(nOde);
  for (int i = 0; i < nOde; i++) biovar[0][i] = 1;

  vector<vector<double> > tlag(1);
  tlag[0].resize(nOde);
  for (int i = 0; i < nOde; i++) tlag[0][i] = 0;

  int nPD = 3;
  double rel_tol = 1e8, abs_tol = 1e-8;
  long int max_num_steps = 1e8;
  Matrix<double, Dynamic, Dynamic>
    x = mixOde1CptModel_rk45(feebackODE(), nPD,
                             time, amt, rate, ii, evid, cmt, addl, ss,
                             parameters, biovar, tlag,
                             0,
                             rel_tol, abs_tol, max_num_steps);

  Matrix<double, Dynamic, Dynamic> amounts(nEvents, nOde);
  amounts << 10000, 0, 0, 0, 0,
    6065.306623, 3786.208, -0.0007126788, -1.985660e-06, -4.076891e-09,
    3678.794399, 5821.649, -0.0023972982, -1.390885e-05, -5.850529e-08,
    2231.301565, 6813.189, -0.0045843443, -4.140146e-05, -2.670477e-07,
    1353.352804, 7188.323, -0.0069950553, -8.713363e-05, -7.646749e-07,
    820.849975, 7205.188, -0.0094660436, -1.520317e-04, -1.698982e-06,
    497.870679, 7019.273, -0.0119035120, -2.360139e-04, -3.219375e-06,
    301.973833, 6723.888, -0.0142554906, -3.384324e-04, -5.470933e-06,
    183.156389, 6374.696, -0.0164950219, -4.583390e-04, -8.591089e-06,
    3.354626, 3716.663, -0.0300528403, -1.915294e-03, -7.813746e-05;

  double rel_err = 1e-5;
  expect_near_matrix_eq(amounts, x, rel_err);

  // Test AutoDiff against FiniteDiff
}
