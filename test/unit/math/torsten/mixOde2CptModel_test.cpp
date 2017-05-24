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

TEST(Torsten, mixOde1Cpt_singleDose) {
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
  double rel_tol = 1e-6, abs_tol = 1e-6;
  long int max_num_steps = 1e6;
  Matrix<double, Dynamic, Dynamic>
    x_rk45 = mixOde2CptModel_rk45(feedbackODE(), nPD,
                                  time, amt, rate, ii, evid, cmt, addl, ss,
                                  parameters, biovar, tlag,
                                  0,
                                  rel_tol, abs_tol, max_num_steps);

  rel_tol = 1e-10, abs_tol = 1e-10;
  max_num_steps = 1e8;
  Matrix<double, Dynamic, Dynamic>
    x_bdf = mixOde2CptModel_bdf(feedbackODE(), nPD,
                                time, amt, rate, ii, evid, cmt, addl, ss,
                                parameters, biovar, tlag,
                                0,
                                rel_tol, abs_tol, max_num_steps);

  // Solution from mrgsolve (uses an LSODA integrator)
  Matrix<double, Dynamic, Dynamic> amounts(nEvents, nOde);
  amounts << 10000, 0, 0, 0, 0, 0,
    6065.306597, 3579.304,  212.1623, -0.0006874417, -1.933282e-06, -3.990995e-09,
    3678.794412, 5177.749,  678.8210, -0.0022297559, -1.318812e-05, -5.608541e-08,
    2231.301599, 5678.265, 1233.3871, -0.0041121287, -3.824790e-05, -2.508047e-07,
    1353.352829, 5597.489, 1787.0134, -0.0060546255, -7.847821e-05, -7.039447e-07,
    820.849983, 5233.332, 2295.7780, -0.0079139199, -1.335979e-04, -1.533960e-06,
    497.870681, 4753.865, 2741.1870, -0.0096246889, -2.025267e-04, -2.852515e-06,
    301.973832, 4250.712, 3118.6808, -0.0111649478, -2.838628e-04, -4.760265e-06,
    183.156387, 3771.009, 3430.9355, -0.0125357580, -3.761421e-04, -7.345522e-06,
    3.354626, 1601.493, 4374.6747, -0.0192607813, -1.370742e-03, -5.951920e-05;

  double rel_err = 1e-5;
  expect_near_matrix_eq(amounts, x_rk45, rel_err);
  expect_near_matrix_eq(amounts, x_bdf, rel_err);

  // Test AutoDiff against FiniteDiff
  double diff = 1e-8, diff2 = 5e-3;
  test_mixOdeCptModel(feedbackODE(), nPD,
                       time, amt, rate, ii, evid, cmt, addl, ss,
                       parameters, biovar, tlag,
                       rel_tol, abs_tol, max_num_steps, diff, diff2, "2_rk45");

  test_mixOdeCptModel(feedbackODE(), nPD,
                      time, amt, rate, ii, evid, cmt, addl, ss,
                      parameters, biovar, tlag,
                      rel_tol, abs_tol, max_num_steps, diff, diff2, "2_bdf");
}


TEST(Torsten, mixOde2Cpt_singleDose_overload) {
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


  double rel_tol = 1e-6, abs_tol = 1e-6;
  long int max_num_steps = 1e6;
  Matrix<double, Eigen::Dynamic, Eigen::Dynamic> x_rk45_122, x_rk45_112, 
    x_rk45_111, x_rk45_121, x_rk45_212, x_rk45_211, x_rk45_221;

  x_rk45_122 = mixOde2CptModel_rk45(feedbackODE(), nPD,
                                time, amt, rate, ii, evid, cmt, addl, ss,
                                parameters[0], biovar, tlag,
                                0,
                                rel_tol, abs_tol, max_num_steps);

  x_rk45_112 = mixOde2CptModel_rk45(feedbackODE(), nPD,
                                    time, amt, rate, ii, evid, cmt, addl, ss,
                                    parameters[0], biovar[0], tlag,
                                    0,
                                    rel_tol, abs_tol, max_num_steps);

  x_rk45_111 = mixOde2CptModel_rk45(feedbackODE(), nPD,
                                    time, amt, rate, ii, evid, cmt, addl, ss,
                                    parameters[0], biovar[0], tlag[0],
                                    0,
                                    rel_tol, abs_tol, max_num_steps);

  x_rk45_121 = mixOde2CptModel_rk45(feedbackODE(), nPD,
                                    time, amt, rate, ii, evid, cmt, addl, ss,
                                    parameters[0], biovar, tlag[0],
                                    0,
                                    rel_tol, abs_tol, max_num_steps);

  x_rk45_212 = mixOde2CptModel_rk45(feedbackODE(), nPD,
                                    time, amt, rate, ii, evid, cmt, addl, ss,
                                    parameters, biovar[0], tlag,
                                    0,
                                    rel_tol, abs_tol, max_num_steps);

  x_rk45_211 = mixOde2CptModel_rk45(feedbackODE(), nPD,
                                    time, amt, rate, ii, evid, cmt, addl, ss,
                                    parameters, biovar[0], tlag[0],
                                    0,
                                    rel_tol, abs_tol, max_num_steps);

  x_rk45_221 = mixOde2CptModel_rk45(feedbackODE(), nPD,
                                    time, amt, rate, ii, evid, cmt, addl, ss,
                                    parameters, biovar, tlag[0],
                                    0,
                                    rel_tol, abs_tol, max_num_steps);

  double rel_tol_bdf = 1e-10, abs_tol_bdf = 1e-10;
  long int max_num_steps_bdf = 1e8;
  Matrix<double, Eigen::Dynamic, Eigen::Dynamic> x_bdf_122, x_bdf_112, 
    x_bdf_111, x_bdf_121, x_bdf_212, x_bdf_211, x_bdf_221;

  x_bdf_122 = mixOde2CptModel_bdf(feedbackODE(), nPD,
                                time, amt, rate, ii, evid, cmt, addl, ss,
                                parameters[0], biovar, tlag,
                                0,
                                rel_tol_bdf, abs_tol_bdf, max_num_steps_bdf);

  x_bdf_112 = mixOde2CptModel_bdf(feedbackODE(), nPD,
                                    time, amt, rate, ii, evid, cmt, addl, ss,
                                    parameters[0], biovar[0], tlag,
                                    0,
                                    rel_tol_bdf, abs_tol_bdf, max_num_steps_bdf);

  x_bdf_111 = mixOde2CptModel_bdf(feedbackODE(), nPD,
                                    time, amt, rate, ii, evid, cmt, addl, ss,
                                    parameters[0], biovar[0], tlag[0],
                                    0,
                                    rel_tol_bdf, abs_tol_bdf, max_num_steps_bdf);

  x_bdf_121 = mixOde2CptModel_bdf(feedbackODE(), nPD,
                                    time, amt, rate, ii, evid, cmt, addl, ss,
                                    parameters[0], biovar, tlag[0],
                                    0,
                                    rel_tol_bdf, abs_tol_bdf, max_num_steps_bdf);

  x_bdf_212 = mixOde2CptModel_bdf(feedbackODE(), nPD,
                                    time, amt, rate, ii, evid, cmt, addl, ss,
                                    parameters, biovar[0], tlag,
                                    0,
                                    rel_tol_bdf, abs_tol_bdf, max_num_steps_bdf);

  x_bdf_211 = mixOde2CptModel_bdf(feedbackODE(), nPD,
                                    time, amt, rate, ii, evid, cmt, addl, ss,
                                    parameters, biovar[0], tlag[0],
                                    0,
                                    rel_tol_bdf, abs_tol_bdf, max_num_steps_bdf);

  x_bdf_221 = mixOde2CptModel_bdf(feedbackODE(), nPD,
                                    time, amt, rate, ii, evid, cmt, addl, ss,
                                    parameters, biovar, tlag[0],
                                    0,
                                    rel_tol_bdf, abs_tol_bdf, max_num_steps_bdf);

  // Solution from mrgsolve
  Matrix<double, Dynamic, Dynamic> amounts(nEvents, nOde);
  amounts << 10000, 0, 0, 0, 0, 0,
    6065.306597, 3579.304,  212.1623, -0.0006874417, -1.933282e-06, -3.990995e-09,
    3678.794412, 5177.749,  678.8210, -0.0022297559, -1.318812e-05, -5.608541e-08,
    2231.301599, 5678.265, 1233.3871, -0.0041121287, -3.824790e-05, -2.508047e-07,
    1353.352829, 5597.489, 1787.0134, -0.0060546255, -7.847821e-05, -7.039447e-07,
    820.849983, 5233.332, 2295.7780, -0.0079139199, -1.335979e-04, -1.533960e-06,
    497.870681, 4753.865, 2741.1870, -0.0096246889, -2.025267e-04, -2.852515e-06,
    301.973832, 4250.712, 3118.6808, -0.0111649478, -2.838628e-04, -4.760265e-06,
    183.156387, 3771.009, 3430.9355, -0.0125357580, -3.761421e-04, -7.345522e-06,
    3.354626, 1601.493, 4374.6747, -0.0192607813, -1.370742e-03, -5.951920e-05;


  double rel_err = 1e-5;
  expect_near_matrix_eq(amounts, x_rk45_122, rel_err);
  expect_near_matrix_eq(amounts, x_rk45_112, rel_err);
  expect_near_matrix_eq(amounts, x_rk45_111, rel_err);
  expect_near_matrix_eq(amounts, x_rk45_121, rel_err);
  expect_near_matrix_eq(amounts, x_rk45_212, rel_err);
  expect_near_matrix_eq(amounts, x_rk45_211, rel_err);
  expect_near_matrix_eq(amounts, x_rk45_221, rel_err);

  double rel_err_bdf = 1e-5;

  expect_near_matrix_eq(amounts, x_bdf_122, rel_err_bdf);
  expect_near_matrix_eq(amounts, x_bdf_112, rel_err_bdf);
  expect_near_matrix_eq(amounts, x_bdf_111, rel_err_bdf);
  expect_near_matrix_eq(amounts, x_bdf_121, rel_err_bdf);
  expect_near_matrix_eq(amounts, x_bdf_212, rel_err_bdf);
  expect_near_matrix_eq(amounts, x_bdf_211, rel_err_bdf);
  expect_near_matrix_eq(amounts, x_bdf_221, rel_err_bdf);
}
