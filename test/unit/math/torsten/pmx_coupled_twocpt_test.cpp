#include <stan/math/rev/mat.hpp>  // FIX ME - includes should be more specific
#include <test/unit/math/torsten/test_util.hpp>
#include <gtest/gtest.h>
#include <test/unit/math/torsten/util_mixOdeCptModel.hpp>
#include <test/unit/math/torsten/expect_near_matrix_eq.hpp>

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

using std::vector;
using Eigen::Matrix;
using Eigen::Dynamic;
using Eigen::MatrixXd;
using torsten::pmx_solve_twocpt_rk45;
using torsten::pmx_solve_twocpt_bdf;

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

TEST(Torsten, mixOde2Cpt_singleDose) {
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
    x_rk45 = torsten::pmx_solve_twocpt_rk45(feedbackODE(), nPD,
                                  time, amt, rate, ii, evid, cmt, addl, ss,
                                  parameters, biovar, tlag,
                                  0,
                                  rel_tol, abs_tol, max_num_steps);

  rel_tol = 1e-10, abs_tol = 1e-10;
  max_num_steps = 1e8;
  Matrix<double, Dynamic, Dynamic>
    x_bdf = torsten::pmx_solve_twocpt_bdf(feedbackODE(), nPD,
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
  MatrixXd xt = amounts.transpose();
  torsten::test::test_val(x_rk45, xt, 1e-5, 1e-5);
  torsten::test::test_val(x_bdf, xt, 2e-5, 1e-5);

  // Test AutoDiff against FiniteDiff
  biovar[0] = std::vector<double>(nOde, 0.8);
  tlag[0] = std::vector<double>(nOde, 1.9);
  TORSTEN_ODE_GRAD_THETA_TEST(pmx_solve_twocpt_rk45, feedbackODE(), nPD,
                              time, amt, rate, ii, evid, cmt, addl, ss,
                              parameters, biovar, tlag,
                              rel_tol, abs_tol, max_num_steps,
                              2e-5, 1e-6, 1e-3, 1e-5);
  TORSTEN_ODE_GRAD_BIOVAR_TEST(pmx_solve_twocpt_rk45, feedbackODE(), nPD,
                               time, amt, rate, ii, evid, cmt, addl, ss,
                               parameters, biovar, tlag,
                               rel_tol, abs_tol, max_num_steps,
                               2e-5, 1e-6, 1e-3, 1e-5);
  TORSTEN_ODE_GRAD_THETA_TEST(pmx_solve_twocpt_bdf, feedbackODE(), nPD,
                              time, amt, rate, ii, evid, cmt, addl, ss,
                              parameters, biovar, tlag,
                              rel_tol, abs_tol, max_num_steps,
                              2e-5, 1e-6, 1e-3, 1e-5);
  TORSTEN_ODE_GRAD_BIOVAR_TEST(pmx_solve_twocpt_bdf, feedbackODE(), nPD,
                               time, amt, rate, ii, evid, cmt, addl, ss,
                               parameters, biovar, tlag,
                               rel_tol, abs_tol, max_num_steps,
                               2e-5, 1e-6, 1e-3, 1e-5);
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

  x_rk45_122 = torsten::pmx_solve_twocpt_rk45(feedbackODE(), nPD,
                                time, amt, rate, ii, evid, cmt, addl, ss,
                                parameters[0], biovar, tlag,
                                0,
                                rel_tol, abs_tol, max_num_steps);

  x_rk45_112 = torsten::pmx_solve_twocpt_rk45(feedbackODE(), nPD,
                                    time, amt, rate, ii, evid, cmt, addl, ss,
                                    parameters[0], biovar[0], tlag,
                                    0,
                                    rel_tol, abs_tol, max_num_steps);

  x_rk45_111 = torsten::pmx_solve_twocpt_rk45(feedbackODE(), nPD,
                                    time, amt, rate, ii, evid, cmt, addl, ss,
                                    parameters[0], biovar[0], tlag[0],
                                    0,
                                    rel_tol, abs_tol, max_num_steps);

  x_rk45_121 = torsten::pmx_solve_twocpt_rk45(feedbackODE(), nPD,
                                    time, amt, rate, ii, evid, cmt, addl, ss,
                                    parameters[0], biovar, tlag[0],
                                    0,
                                    rel_tol, abs_tol, max_num_steps);

  x_rk45_212 = torsten::pmx_solve_twocpt_rk45(feedbackODE(), nPD,
                                    time, amt, rate, ii, evid, cmt, addl, ss,
                                    parameters, biovar[0], tlag,
                                    0,
                                    rel_tol, abs_tol, max_num_steps);

  x_rk45_211 = torsten::pmx_solve_twocpt_rk45(feedbackODE(), nPD,
                                    time, amt, rate, ii, evid, cmt, addl, ss,
                                    parameters, biovar[0], tlag[0],
                                    0,
                                    rel_tol, abs_tol, max_num_steps);

  x_rk45_221 = torsten::pmx_solve_twocpt_rk45(feedbackODE(), nPD,
                                    time, amt, rate, ii, evid, cmt, addl, ss,
                                    parameters, biovar, tlag[0],
                                    0,
                                    rel_tol, abs_tol, max_num_steps);

  double rel_tol_bdf = 1e-10, abs_tol_bdf = 1e-10;
  long int max_num_steps_bdf = 1e8;
  Matrix<double, Eigen::Dynamic, Eigen::Dynamic> x_bdf_122, x_bdf_112,
    x_bdf_111, x_bdf_121, x_bdf_212, x_bdf_211, x_bdf_221;

  x_bdf_122 = torsten::pmx_solve_twocpt_bdf(feedbackODE(), nPD,
                                time, amt, rate, ii, evid, cmt, addl, ss,
                                parameters[0], biovar, tlag,
                                0,
                                rel_tol_bdf, abs_tol_bdf, max_num_steps_bdf);

  x_bdf_112 = torsten::pmx_solve_twocpt_bdf(feedbackODE(), nPD,
                                    time, amt, rate, ii, evid, cmt, addl, ss,
                                    parameters[0], biovar[0], tlag,
                                    0,
                                    rel_tol_bdf, abs_tol_bdf, max_num_steps_bdf);

  x_bdf_111 = torsten::pmx_solve_twocpt_bdf(feedbackODE(), nPD,
                                    time, amt, rate, ii, evid, cmt, addl, ss,
                                    parameters[0], biovar[0], tlag[0],
                                    0,
                                    rel_tol_bdf, abs_tol_bdf, max_num_steps_bdf);

  x_bdf_121 = torsten::pmx_solve_twocpt_bdf(feedbackODE(), nPD,
                                    time, amt, rate, ii, evid, cmt, addl, ss,
                                    parameters[0], biovar, tlag[0],
                                    0,
                                    rel_tol_bdf, abs_tol_bdf, max_num_steps_bdf);

  x_bdf_212 = torsten::pmx_solve_twocpt_bdf(feedbackODE(), nPD,
                                    time, amt, rate, ii, evid, cmt, addl, ss,
                                    parameters, biovar[0], tlag,
                                    0,
                                    rel_tol_bdf, abs_tol_bdf, max_num_steps_bdf);

  x_bdf_211 = torsten::pmx_solve_twocpt_bdf(feedbackODE(), nPD,
                                    time, amt, rate, ii, evid, cmt, addl, ss,
                                    parameters, biovar[0], tlag[0],
                                    0,
                                    rel_tol_bdf, abs_tol_bdf, max_num_steps_bdf);

  x_bdf_221 = torsten::pmx_solve_twocpt_bdf(feedbackODE(), nPD,
                                    time, amt, rate, ii, evid, cmt, addl, ss,
                                    parameters, biovar, tlag[0],
                                    0,
                                    rel_tol_bdf, abs_tol_bdf, max_num_steps_bdf);

  // Solution from mrgsolve
  MatrixXd amounts(nEvents, nOde);
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
  MatrixXd xt = amounts.transpose();

  torsten::test::test_val(xt, x_rk45_122, 1e-5, 1e-7);
  torsten::test::test_val(xt, x_rk45_112, 1e-5, 1e-7);
  torsten::test::test_val(xt, x_rk45_111, 1e-5, 1e-7);
  torsten::test::test_val(xt, x_rk45_121, 1e-5, 1e-7);
  torsten::test::test_val(xt, x_rk45_212, 1e-5, 1e-7);
  torsten::test::test_val(xt, x_rk45_211, 1e-5, 1e-7);
  torsten::test::test_val(xt, x_rk45_221, 1e-5, 1e-7);

  torsten::test::test_val(xt, x_bdf_122, 1e-5, 1e-7);
  torsten::test::test_val(xt, x_bdf_112, 1e-5, 1e-7);
  torsten::test::test_val(xt, x_bdf_111, 1e-5, 1e-7);
  torsten::test::test_val(xt, x_bdf_121, 1e-5, 1e-7);
  torsten::test::test_val(xt, x_bdf_212, 1e-5, 1e-7);
  torsten::test::test_val(xt, x_bdf_211, 1e-5, 1e-7);
  torsten::test::test_val(xt, x_bdf_221, 1e-5, 1e-7);
}

TEST(Torsten, mixOde2Cpt_rate) {
  using std::vector;
  using Eigen::Matrix;
  using Eigen::Dynamic;

  int nEvents = 10;
  vector<double> time(nEvents);
  time[0] = 0.0;
  for (int i = 1; i < (nEvents - 1); i++) time[i] = time[i - 1] + 0.25;
  time[9] = 4.0;
  vector<double> amt(nEvents, 0);
  amt[0] = 12000;
  vector<double> rate(nEvents, 0);
  rate[0] = 12000;
  vector<int> cmt(nEvents, 2);
  cmt[0] = 1;
  vector<int> evid(nEvents, 0);
  evid[0] = 1;
  vector<double> ii(nEvents, 0);
  ii[0] = 12;
  vector<int> addl(nEvents, 0);
  addl[0] = 14;
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
  MatrixXd
    x_rk45 = torsten::pmx_solve_twocpt_rk45(feedbackODE(), nPD,
                                  time, amt, rate, ii, evid, cmt, addl, ss,
                                  parameters, biovar, tlag,
                                  0,
                                  rel_tol, abs_tol, max_num_steps);

  rel_tol = 1e-10, abs_tol = 1e-10;
  max_num_steps = 1e8;
  MatrixXd
    x_bdf = torsten::pmx_solve_twocpt_bdf(feedbackODE(), nPD,
                                time, amt, rate, ii, evid, cmt, addl, ss,
                                parameters, biovar, tlag,
                                0,
                                rel_tol, abs_tol, max_num_steps);

  // Solution from mrgsolve (uses an LSODA integrator)
  MatrixXd amounts(nEvents, nOde);
  amounts << 0.00000,    0.0000,    0.00000,  0.0000000000,  0.000000,  0.000000,
             2360.81604,  601.5528,   22.49548, -0.0000726505, -1.499109e-07, -2.448992e-10,
             3792.72335, 1951.4716,  152.31877, -0.0004967093, -2.110375e-06, -7.029932e-09,
             4661.21904, 3599.5932,  438.32811, -0.0014439180, -9.454212e-06, -4.809775e-08,
             5187.98830, 5301.0082,  892.11236, -0.0029697950, -2.658409e-05, -1.833674e-07,
             3146.67402, 6328.6148, 1483.47278, -0.0049954464, -5.788643e-05, -5.079766e-07,
             1908.55427, 6478.2514, 2110.88020, -0.0072059086, -1.060159e-04, -1.145673e-06,
             1157.59667, 6180.6685, 2705.53777, -0.0093809310, -1.713345e-04, -2.230652e-06,
             702.11786, 5681.5688, 3235.76051, -0.0114135989, -2.529409e-04, -3.893275e-06,
             12.85974, 2316.1858, 5156.18535, -0.0216237650, -1.312811e-03, -4.956987e-05;
  MatrixXd xt = amounts.transpose();
  torsten::test::test_val(xt, x_rk45, 1e-4, 1e-8);
  torsten::test::test_val(xt, x_bdf,  1e-4, 1e-8);

  // Test AutoDiff against FiniteDiff
  // FIX ME - lag time
  TORSTEN_ODE_GRAD_THETA_TEST(pmx_solve_twocpt_rk45, feedbackODE(), nPD,
                              time, amt, rate, ii, evid, cmt, addl, ss,
                              parameters, biovar, tlag,
                              rel_tol, abs_tol, max_num_steps,
                              2e-5, 1e-6, 1e-3, 1e-5);
  TORSTEN_ODE_GRAD_BIOVAR_TEST(pmx_solve_twocpt_rk45, feedbackODE(), nPD,
                               time, amt, rate, ii, evid, cmt, addl, ss,
                               parameters, biovar, tlag,
                               rel_tol, abs_tol, max_num_steps,
                               2e-5, 1e-6, 1e-3, 1e-5);
  // TORSTEN_ODE_GRAD_THETA_TEST(pmx_solve_twocpt_bdf, feedbackODE(), nPD,
  //                             time, amt, rate, ii, evid, cmt, addl, ss,
  //                             parameters, biovar, tlag,
  //                             rel_tol, abs_tol, max_num_steps,
  //                             2e-5, 1e-6, 1e-3, 1e-5);
  // TORSTEN_ODE_GRAD_BIOVAR_TEST(pmx_solve_twocpt_bdf, feedbackODE(), nPD,
  //                              time, amt, rate, ii, evid, cmt, addl, ss,
  //                              parameters, biovar, tlag,
  //                              rel_tol, abs_tol, max_num_steps,
  //                              2e-5, 1e-6, 1e-3, 1e-5);
}

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

TEST(Torsten, mixOde2Cpt_SS_bolus) {
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
  MatrixXd
    x_rk45 = torsten::pmx_solve_twocpt_rk45(feedbackODE(), nPD,
                                  time, amt, rate, ii, evid, cmt, addl, ss,
                                  parameters, biovar, tlag,
                                  0,
                                  rel_tol_rk, abs_tol_rk, max_num_steps_rk);

  double rel_tol_bdf = 1e-10, abs_tol_bdf = 1e-10;
  long int max_num_steps_bdf = 1e8;
  MatrixXd
    x_bdf = torsten::pmx_solve_twocpt_bdf(feedbackODE(), nPD,
                                time, amt, rate, ii, evid, cmt, addl, ss,
                                parameters, biovar, tlag,
                                0,
                                rel_tol_bdf, abs_tol_bdf, max_num_steps_bdf);

  // Solution from mrgsolve (uses an LSODA integrator)
  MatrixXd amounts(nEvents, nOde);
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
  MatrixXd xt = amounts.transpose();
  torsten::test::test_val(xt, x_rk45, 1.e-3, 1.e-8);
  torsten::test::test_val(xt, x_bdf,  1.e-3, 1.e-8);

  // Test AutoDiff against FiniteDiff
  TORSTEN_ODE_GRAD_THETA_TEST(pmx_solve_twocpt_rk45, feedbackODE(), nPD,
                              time, amt, rate, ii, evid, cmt, addl, ss,
                              parameters, biovar, tlag,
                              rel_tol_rk, abs_tol_rk, max_num_steps_rk,
                              1e-5, 1e-6, 1e-1, 1e-3);
  TORSTEN_ODE_GRAD_BIOVAR_TEST(pmx_solve_twocpt_rk45, feedbackODE(), nPD,
                               time, amt, rate, ii, evid, cmt, addl, ss,
                               parameters, biovar, tlag,
                               rel_tol_rk, abs_tol_rk, max_num_steps_rk,
                               1e-5, 1e-6, 1e-1, 1e-3);
  TORSTEN_ODE_GRAD_THETA_TEST(pmx_solve_twocpt_bdf, feedbackODE(), nPD,
                              time, amt, rate, ii, evid, cmt, addl, ss,
                              parameters, biovar, tlag,
                              rel_tol_bdf, abs_tol_bdf, max_num_steps_bdf,
                              1e-6, 1e-6, 8e-1, 1e-3);
  TORSTEN_ODE_GRAD_BIOVAR_TEST(pmx_solve_twocpt_bdf, feedbackODE(), nPD,
                               time, amt, rate, ii, evid, cmt, addl, ss,
                               parameters, biovar, tlag,
                               rel_tol_bdf, abs_tol_bdf, max_num_steps_bdf,
                               1e-6, 1e-6, 1e-1, 1e-3);
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
  MatrixXd
    x_rk45 = torsten::pmx_solve_twocpt_rk45(feedbackODE(), nPD,
                                  time, amt, rate, ii, evid, cmt, addl, ss,
                                  parameters, biovar, tlag,
                                  0,
                                  rel_tol_rk, abs_tol_rk, max_num_steps_rk);

  //std::cout << x_rk45 << std::endl << std::endl;

  double rel_tol_bdf = 1e-10, abs_tol_bdf = 1e-10;
  double max_num_steps_bdf = 1e8;
  MatrixXd
    x_bdf = torsten::pmx_solve_twocpt_bdf(feedbackODE(), nPD,
                                time, amt, rate, ii, evid, cmt, addl, ss,
                                parameters, biovar, tlag,
                                0,
                                rel_tol_bdf, abs_tol_bdf, max_num_steps_bdf);

  // Solution from mrgsolve (uses an LSODA integrator)
  MatrixXd amounts(nEvents, nOde);
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
  MatrixXd xt = amounts.transpose();
  torsten::test::test_val(xt, x_rk45, 1.e-3, 1.e-8);
  torsten::test::test_val(xt, x_bdf,  1.e-3, 1.e-8);

  // Test AutoDiff against FiniteDiff
  // diff2 obtained empirically.
  TORSTEN_ODE_GRAD_THETA_TEST(pmx_solve_twocpt_rk45, feedbackODE(), nPD,
                              time, amt, rate, ii, evid, cmt, addl, ss,
                              parameters, biovar, tlag,
                              rel_tol_rk, abs_tol_rk, max_num_steps_rk,
                              2e-5, 1e-6, 1e-5, 1e-5);
  // FIXME: steady state exception
  // TORSTEN_ODE_GRAD_BIOVAR_TEST(pmx_solve_onecpt_rk45, feedbackODE(), nPD,
  //                              time, amt, rate, ii, evid, cmt, addl, ss,
  //                              parameters, biovar, tlag,
  //                              rel_tol, abs_tol, max_num_steps,
  //                              2e-5, 1e-6, 1e-5, 1e-5)
  TORSTEN_ODE_GRAD_THETA_TEST(pmx_solve_twocpt_bdf, feedbackODE(), nPD,
                              time, amt, rate, ii, evid, cmt, addl, ss,
                              parameters, biovar, tlag,
                              rel_tol_bdf, abs_tol_bdf, max_num_steps_bdf,
                              2e-5, 1e-6, 1e-3, 1e-5);
  // FIXME: steady state exception
  // TORSTEN_ODE_GRAD_BIOVAR_TEST(pmx_solve_onecpt_bdf, feedbackODE(), nPD,
  //                              time, amt, rate, ii, evid, cmt, addl, ss,
  //                              parameters, biovar, tlag,
  //                              rel_tol, abs_tol, max_num_steps,
  //                              2e-5, 1e-6, 1e-3, 1e-5)
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
  MatrixXd
    x_rk45 = torsten::pmx_solve_twocpt_rk45(feedbackODE(), nPD,
                                  time, amt, rate, ii, evid, cmt, addl, ss,
                                  parameters, biovar, tlag,
                                  0,
                                  rel_tol_rk, abs_tol_rk, max_num_steps_rk);

  double rel_tol_bdf = 1e-10, abs_tol_bdf = 1e-10;
  double max_num_steps_bdf = 1e8;
  MatrixXd
    x_bdf = torsten::pmx_solve_twocpt_bdf(feedbackODE(), nPD,
                                time, amt, rate, ii, evid, cmt, addl, ss,
                                parameters, biovar, tlag,
                                0,
                                rel_tol_bdf, abs_tol_bdf, max_num_steps_bdf);

  // can't do constant rate in mrgsolve. Comparing to result obtained
  // with generalOdeModel, as a provisional test.
  MatrixXd
    x = torsten::pmx_solve_rk45(fullODE(), nOde,
                             time, amt, rate, ii, evid, cmt, addl, ss,
                             parameters, biovar, tlag,
                             0,
                             rel_tol_bdf, abs_tol_bdf, max_num_steps_bdf);

  torsten::test::test_val(x, x_rk45, 1.e-5, 1.e-8);
  torsten::test::test_val(x, x_bdf,  1.e-5, 1.e-8);

  TORSTEN_ODE_GRAD_THETA_TEST(pmx_solve_twocpt_rk45, feedbackODE(), nPD,
                              time, amt, rate, ii, evid, cmt, addl, ss,
                              parameters, biovar, tlag,
                              rel_tol_rk, abs_tol_rk, max_num_steps_rk,
                              2e-5, 1e-6, 1e-5, 1e-5);
  // FIXME: steady state exception
  // TORSTEN_ODE_GRAD_BIOVAR_TEST(pmx_solve_onecpt_rk45, feedbackODE(), nPD,
  //                              time, amt, rate, ii, evid, cmt, addl, ss,
  //                              parameters, biovar, tlag,
  //                              rel_tol, abs_tol, max_num_steps,
  //                              2e-5, 1e-6, 1e-5, 1e-5)
  TORSTEN_ODE_GRAD_THETA_TEST(pmx_solve_twocpt_bdf, feedbackODE(), nPD,
                              time, amt, rate, ii, evid, cmt, addl, ss,
                              parameters, biovar, tlag,
                              rel_tol_bdf, abs_tol_bdf, max_num_steps_bdf,
                              2e-5, 1e-6, 1e-3, 1e-5);
  // FIXME: steady state exception
  // TORSTEN_ODE_GRAD_BIOVAR_TEST(pmx_solve_onecpt_bdf, feedbackODE(), nPD,
  //                              time, amt, rate, ii, evid, cmt, addl, ss,
  //                              parameters, biovar, tlag,
  //                              rel_tol, abs_tol, max_num_steps,
  //                              2e-5, 1e-6, 1e-3, 1e-5)
}
