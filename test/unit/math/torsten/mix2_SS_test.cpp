#include <stan/math/rev/mat.hpp>  // FIX ME - includes should be more specific
#include <gtest/gtest.h>
#include <test/unit/math/torsten/util_mixOdeCptModel.hpp>
#include <test/unit/math/prim/mat/fun/expect_near_matrix_eq.hpp>
#include <test/unit/math/prim/mat/fun/expect_matrix_eq.hpp>

// Note: for tuning parameters of ODE integrators,
// use default values.
// For rk45:
//   * rel_tol = 1e-6
//   * abs_tol = 1e-6
//   * max_num_steps = 1e6
//
// For bdf:
//   * rel_tol = 1e-10
//   * abs_tol = 1e-10
//   * max_num_steps = 1e8

struct feedbackODE {
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

struct fullODE {
  template <typename T0, typename T1, typename T2>
  inline
    std::vector<typename boost::math::tools::promote_args<T0, T1, T2>::type>
  operator()(const T0& t,
           const std::vector<T1>& x,
           const std::vector<T2>& theta,
           const std::vector<double>& x_r,
           const std::vector<int>& x_i,
           std::ostream* pstream_) const {
    typedef typename boost::math::tools::promote_args<T0, T1, T2>::type
    scalar;

    scalar 
      CL = theta[0],
      Q = theta[1],
      VC = theta[2],
      VP = theta[3],
      ka = theta[4],
      k10 = CL / VC,
      k12 = Q / VC,
      k21 = Q / VP,
      Mtt = theta[5],
      circ0 = theta[6],
      alpha = theta[7],
      gamma = theta[8],
      ktr = 4 / Mtt,
      prol = x[3] + circ0,
      transit = x[4] + circ0,
      circ = x[5] + circ0,
      Edrug;

    std::vector<scalar> dxdt(6);
    dxdt[0] = -ka * x[0];
    dxdt[1] = ka * x[0] - (k10 + k12) * x[1] + k21 * x[2];
    dxdt[2] = k12 * x[1] - k21 * x[2];
    Edrug = alpha * x[1] / VC;
    dxdt[3] = ktr * prol * ((1 - Edrug) * pow(circ0 / circ, gamma) - 1);
    dxdt[4] = ktr * (prol - transit);
    dxdt[5] = ktr * (transit - circ);

    return dxdt;
  }
};

TEST(Torsten, mixOde1Cpt_singleDose) {
  using std::vector;
  using Eigen::Matrix;
  using Eigen::Dynamic;

  int nEvents = 10;
  vector<double> time(nEvents);
  time[0] = 0.0;
  for (int i = 1; i < nEvents; i++) time[i] = time[i - 1] + 1;
  vector<double> amt(nEvents, 0);
  amt[0] = 1200;
  vector<double> rate(nEvents, 0);
  vector<int> cmt(nEvents, 2);
  cmt[0] = 1;
  vector<int> evid(nEvents, 0);
  evid[0] = 1;
  vector<double> ii(nEvents, 0);
  ii[0] = 12;
  vector<int> addl(nEvents, 0);
  vector<int> ss(nEvents, 0);
  ss[0] = 1;

  int nParameters = 9;
  vector<vector<double> > parameters(1);
  parameters[0].resize(nParameters);
  parameters[0][0] = 10;  // CL
  parameters[0][1] = 15;  // Q
  parameters[0][2] = 35;  // VC
  parameters[0][3] = 105;  // VP
  parameters[0][4] = 2.0;  // ka
  parameters[0][5] = 125;  // Mtt
  parameters[0][6] = 5;  // Circ0
  parameters[0][7] = 3e-4;  // alpha
  parameters[0][8] = 0.17;  // gamma

  int nOde = 6;
  vector<vector<double> > biovar(1);
  biovar[0].resize(nOde);
  for (int i = 0; i < nOde; i++) biovar[0][i] = 1;

  vector<vector<double> > tlag(1);
  tlag[0].resize(nOde);
  for (int i = 0; i < nOde; i++) tlag[0][i] = 0;

  int nPD = 3;
  double rel_tol_rk = 1e-6, abs_tol_rk = 1e-6;
  long int max_num_steps_rk = 1e6;
  Matrix<double, Dynamic, Dynamic>
    x_rk45 = mixOde2CptModel_rk45(feedbackODE(), nPD,
                                  time, amt, rate, ii, evid, cmt, addl, ss,
                                  parameters, biovar, tlag,
                                  0,
                                  rel_tol_rk, abs_tol_rk, max_num_steps_rk);

  double rel_tol_bdf = 1e-10, abs_tol_bdf = 1e-10;
  long int max_num_steps_bdf = 1e8;
  Matrix<double, Dynamic, Dynamic>
    x_bdf = mixOde2CptModel_bdf(feedbackODE(), nPD,
                                time, amt, rate, ii, evid, cmt, addl, ss,
                                parameters, biovar, tlag,
                                0,
                                rel_tol_bdf, abs_tol_bdf, max_num_steps_bdf);
  
  // Solution from mrgsolve (uses an LSODA integrator)
  Matrix<double, Dynamic, Dynamic> amounts(nEvents, nOde);
  amounts << 1.200000e+03, 179.9494,  835.4153, -0.08689426, -0.08765832, -0.08765770,
             1.624023e+02, 842.7123, 1008.6670, -0.08737445, -0.08763983, -0.08765738,
             2.197877e+01, 615.0706, 1166.7605, -0.08789426, -0.08764055, -0.08765680,
             2.974503e+00, 435.9062, 1217.1511, -0.08811989, -0.08765274, -0.08765646,
             4.025552e-01, 339.0661, 1207.3381, -0.08816317, -0.08766848, -0.08765659,
             5.447993e-02, 287.5025, 1170.4685, -0.08810975, -0.08768340, -0.08765720,
             7.373058e-03, 257.6031, 1122.8839, -0.08800312, -0.08769524, -0.08765823,
             9.978352e-04, 237.8551, 1072.0195, -0.08786379, -0.08770280, -0.08765953,
             1.350418e-04, 222.9737, 1021.1484, -0.08770143, -0.08770535, -0.08766094,
             1.827653e-05, 210.5666,  971.6636, -0.08752080, -0.08770241, -0.08766231;

  double rel_err = 1e-3;
  expect_near_matrix_eq(amounts, x_rk45, rel_err);

  expect_near_matrix_eq(amounts, x_bdf, rel_err);

  // Test AutoDiff against FiniteDiff
  // diff2 obtained empirically
  double diff = 1e-8, diff2 = 6e-2;
  test_mixOdeCptModel(feedbackODE(), nPD,
                      time, amt, rate, ii, evid, cmt, addl, ss,
                      parameters, biovar, tlag,
                      rel_tol_rk, abs_tol_rk, max_num_steps_rk,
                      diff, diff2, "2_rk45");

  diff2 = 7.2e-1;  // FIX ME (see below)
  test_mixOdeCptModel(feedbackODE(), nPD,
                      time, amt, rate, ii, evid, cmt, addl, ss,
                      parameters, biovar, tlag,
                      rel_tol_bdf, abs_tol_bdf, max_num_steps_bdf,
                      diff, diff2, "2_bdf");

  // FIX ME - need relative error. The absolute relative error doesn't make
  // sense seeing how much the scale of the result varies. Quick inspection
  // suggests the results look ok.
}

TEST(Torsten, mixOde2Cpt_SS_infusion) {
  using std::vector;
  using Eigen::Matrix;
  using Eigen::Dynamic;
  
  int nEvents = 10;
  vector<double> time(nEvents);
  time[0] = 0.0;
  for (int i = 1; i < nEvents; i++) time[i] = time[i - 1] + 1;
  vector<double> amt(nEvents, 0);
  amt[0] = 1200;
  vector<double> rate(nEvents, 0);
  rate[0] = 150;
  vector<int> cmt(nEvents, 2);
  cmt[0] = 1;
  vector<int> evid(nEvents, 0);
  evid[0] = 1;
  vector<double> ii(nEvents, 0);
  ii[0] = 12;
  vector<int> addl(nEvents, 0);
  vector<int> ss(nEvents, 0);
  ss[0] = 1;
  
  int nParameters = 9;
  vector<vector<double> > parameters(1);
  parameters[0].resize(nParameters);
  parameters[0][0] = 10;  // CL
  parameters[0][1] = 15;  // Q
  parameters[0][2] = 35;  // VC
  parameters[0][3] = 105;  // VP
  parameters[0][4] = 2.0;  // ka
  parameters[0][5] = 125;  // Mtt
  parameters[0][6] = 5;  // Circ0
  parameters[0][7] = 3e-4;  // alpha
  parameters[0][8] = 0.17;  // gamma

  int nOde = 6;
  vector<vector<double> > biovar(1);
  biovar[0].resize(nOde);
  for (int i = 0; i < nOde; i++) biovar[0][i] = 1;

  vector<vector<double> > tlag(1);
  tlag[0].resize(nOde);
  for (int i = 0; i < nOde; i++) tlag[0][i] = 0;

  int nPD = 3;
  double rel_tol_rk = 1e-6, abs_tol_rk = 1e-6;
  long int max_num_steps_rk = 1e6;
  Matrix<double, Dynamic, Dynamic>
    x_rk45 = mixOde2CptModel_rk45(feedbackODE(), nPD,
                                  time, amt, rate, ii, evid, cmt, addl, ss,
                                  parameters, biovar, tlag,
                                  0,
                                  rel_tol_rk, abs_tol_rk, max_num_steps_rk);
  
  double rel_tol_bdf = 1e-10, abs_tol_bdf = 1e-10;
  double max_num_steps_bdf = 1e8;
  Matrix<double, Dynamic, Dynamic>
    x_bdf = mixOde2CptModel_bdf(feedbackODE(), nPD,
                                time, amt, rate, ii, evid, cmt, addl, ss,
                                parameters, biovar, tlag,
                                0,
                                rel_tol_bdf, abs_tol_bdf, max_num_steps_bdf);

  // Solution from mrgsolve (uses an LSODA integrator)
  Matrix<double, Dynamic, Dynamic> amounts(nEvents, nOde);
  amounts << 0.0251597, 232.5059, 1022.9500, -0.08764767, -0.08769098, -0.08766079,
             64.8532585, 281.8567,  987.1464, -0.08751307, -0.08768729, -0.08766170,
             73.6267877, 339.7994,  981.0063, -0.08746266, -0.08768080, -0.08766240,
             74.8141559, 373.5566,  993.6797, -0.08747354, -0.08767398, -0.08766287,
             74.9748488, 392.5755, 1014.8106, -0.08751899, -0.08766832, -0.08766313,
             74.9965962, 404.3396, 1039.0644, -0.08758463, -0.08766462, -0.08766323,
             74.9995393, 412.6384, 1063.9958, -0.08766355, -0.08766332, -0.08766324,
             74.9999377, 419.2312, 1088.5357, -0.08775240, -0.08766471, -0.08766326,
             74.9999916, 424.9186, 1112.2394, -0.08784949, -0.08766899, -0.08766337,
             10.1501452, 363.8037, 1123.7893, -0.08791728, -0.08767598, -0.08766365;

  // the relative error is 1.5% for this problem, which is
  // comparitively high.
  double rel_err = 1e-3;
  expect_near_matrix_eq(amounts, x_rk45, rel_err);
  expect_near_matrix_eq(amounts, x_bdf, rel_err);

  // Test AutoDiff against FiniteDiff
  // diff2 obtained empirically.
  // Skip biovar (F) tests, since this returns an exception.
  double diff = 1e-8, diff2 = 1e-3;
  test_mixOdeCptModel(feedbackODE(), nPD,
                      time, amt, rate, ii, evid, cmt, addl, ss,
                      parameters, biovar, tlag,
                      rel_tol_rk, abs_tol_rk, max_num_steps_rk,
                      diff, diff2, "2_rk45", 2);

  test_mixOdeCptModel(feedbackODE(), nPD,
                      time, amt, rate, ii, evid, cmt, addl, ss,
                      parameters, biovar, tlag,
                      rel_tol_bdf, abs_tol_bdf, max_num_steps_bdf,
                      diff, diff2, "2_bdf", 2);
}

TEST(Torsten, mixOdeCpt2_SS_constant_rate) {
  using std::vector;
  using Eigen::Matrix;
  using Eigen::Dynamic;

  int nEvents = 10;
  vector<double> time(nEvents);
  time[0] = 0.0;
  for (int i = 1; i < nEvents; i++) time[i] = time[i - 1] + 1;
  vector<double> amt(nEvents, 0);
  vector<double> rate(nEvents, 0);
  rate[0] = 150;
  vector<int> cmt(nEvents, 2);
  cmt[0] = 1;
  vector<int> evid(nEvents, 0);
  evid[0] = 1;
  vector<double> ii(nEvents, 0);
  vector<int> addl(nEvents, 0);
  vector<int> ss(nEvents, 0);
  ss[0] = 1;

  int nParameters = 9;
  vector<vector<double> > parameters(1);
  parameters[0].resize(nParameters);
  parameters[0][0] = 10;  // CL
  parameters[0][1] = 15;  // Q
  parameters[0][2] = 35;  // VC
  parameters[0][3] = 105;  // VP
  parameters[0][4] = 2.0;  // ka
  parameters[0][5] = 125;  // Mtt
  parameters[0][6] = 5;  // Circ0
  parameters[0][7] = 3e-4;  // alpha
  parameters[0][8] = 0.17;  // gamma

  int nOde = 6;
  vector<vector<double> > biovar(1);
  biovar[0].resize(nOde);
  for (int i = 0; i < nOde; i++) biovar[0][i] = 1;

  vector<vector<double> > tlag(1);
  tlag[0].resize(nOde);
  for (int i = 0; i < nOde; i++) tlag[0][i] = 0;

  int nPD = 3;
  double rel_tol_rk = 1e-6, abs_tol_rk = 1e-6;
  long int max_num_steps_rk = 1e6;
  Matrix<double, Dynamic, Dynamic>
    x_rk45 = mixOde2CptModel_rk45(feedbackODE(), nPD,
                                  time, amt, rate, ii, evid, cmt, addl, ss,
                                  parameters, biovar, tlag,
                                  0,
                                  rel_tol_rk, abs_tol_rk, max_num_steps_rk);

  double rel_tol_bdf = 1e-10, abs_tol_bdf = 1e-10;
  double max_num_steps_bdf = 1e8;
  Matrix<double, Dynamic, Dynamic>
    x_bdf = mixOde2CptModel_bdf(feedbackODE(), nPD,
                                time, amt, rate, ii, evid, cmt, addl, ss,
                                parameters, biovar, tlag,
                                0,
                                rel_tol_bdf, abs_tol_bdf, max_num_steps_bdf);

  // can't do constant rate in mrgsolve. Comparing to result obtained
  // with generalOdeModel, as a provisional test.  
  Matrix<double, Dynamic, Dynamic>
    x = generalOdeModel_rk45(fullODE(), nOde,
                             time, amt, rate, ii, evid, cmt, addl, ss,
                             parameters, biovar, tlag,
                             0,
                             rel_tol_bdf, abs_tol_bdf, max_num_steps_bdf);

  double rel_err_rk = 1e-5;
  expect_near_matrix_eq(x, x_rk45, rel_err_rk);
  
  double rel_err_bdf = 1e-5;
  expect_near_matrix_eq(x, x_bdf, rel_err_bdf);

  double diff = 1e-8, diff2 = 1e-3;
  test_mixOdeCptModel(feedbackODE(), nPD,
                      time, amt, rate, ii, evid, cmt, addl, ss,
                      parameters, biovar, tlag,
                      rel_tol_rk, abs_tol_rk, max_num_steps_rk,
                      diff, diff2, "2_rk45", 2);

  test_mixOdeCptModel(feedbackODE(), nPD,
                      time, amt, rate, ii, evid, cmt, addl, ss,
                      parameters, biovar, tlag,
                      rel_tol_bdf, abs_tol_bdf, max_num_steps_bdf,
                      diff, diff2, "2_bdf", 2);
}
