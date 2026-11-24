#include <stan/math/rev/core.hpp>
#include <stan/math/prim/fun/Eigen.hpp>
#include <stan/math/torsten/pmx_solve_onecpt_rk45.hpp>
#include <stan/math/torsten/pmx_solve_rk45.hpp>
#include <stan/math/torsten/pmx_solve_onecpt_bdf.hpp>
#include <stan/math/torsten/test/unit/test_macros.hpp>
#include <stan/math/torsten/test/unit/util_mixOdeCptModel.hpp>
#include <stan/math/torsten/test/unit/pmx_coupled_model_fixture.hpp>
#include <stan/math/torsten/test/unit/test_macros.hpp>
#include <gtest/gtest.h>

using std::vector;
using Eigen::Matrix;
using Eigen::Dynamic;
using Eigen::MatrixXd;
using torsten::pmx_solve_onecpt_rk45;
using torsten::pmx_solve_rk45;
using torsten::pmx_solve_onecpt_bdf;

TEST_F(TorstenCoupledOneCptTest, single_bolus) {
  double rel_tol = 1e-6, abs_tol = 1e-6;
  long int max_num_steps = 1e6;
  CoupledOneCptODE f;

  MatrixXd x_rk45 = pmx_solve_onecpt_rk45(f, nPD,
                                         time, amt, rate, ii, evid, cmt, addl, ss,
                                         parameters, biovar, tlag,
                                         rel_tol, abs_tol, max_num_steps);

  rel_tol = 1e-10, abs_tol = 1e-10;
  max_num_steps = 1e8;
  MatrixXd x_bdf = pmx_solve_onecpt_bdf(f, nPD,
                                       time, amt, rate, ii, evid, cmt, addl, ss,
                                       parameters, biovar, tlag,
                                       rel_tol, abs_tol, max_num_steps);

  // Solution from mrgsolve (uses an LSODA integrator)
  MatrixXd amounts(nt, nOde);
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
  MatrixXd xt = amounts.transpose();

  torsten::test::test_val(x_rk45, xt, 1e-5, 1e-5);
  torsten::test::test_val(x_bdf, xt, 2e-5, 1e-5);

  // Test AD against finite diff
  biovar[0] = std::vector<double>(nOde, 0.8);
  tlag[0] = std::vector<double>(nOde, 1.9);
  TORSTEN_ODE_GRAD_THETA_TEST(pmx_solve_onecpt_rk45, f, nPD,
                              time, amt, rate, ii, evid, cmt, addl, ss,
                              parameters, biovar, tlag,
                              rel_tol, abs_tol, max_num_steps,
                              2e-5, 1e-6, 1e-5, 1e-5);
  TORSTEN_ODE_GRAD_BIOVAR_TEST(pmx_solve_onecpt_rk45, f, nPD,
                               time, amt, rate, ii, evid, cmt, addl, ss,
                               parameters, biovar, tlag,
                               rel_tol, abs_tol, max_num_steps,
                               2e-5, 1e-6, 1e-5, 1e-5);
  TORSTEN_ODE_GRAD_THETA_TEST(pmx_solve_onecpt_bdf, f, nPD,
                              time, amt, rate, ii, evid, cmt, addl, ss,
                              parameters, biovar, tlag,
                              rel_tol, abs_tol, max_num_steps,
                              2e-5, 1e-6, 1e-3, 1e-5);
  TORSTEN_ODE_GRAD_BIOVAR_TEST(pmx_solve_onecpt_bdf, f, nPD,
                               time, amt, rate, ii, evid, cmt, addl, ss,
                               parameters, biovar, tlag,
                               rel_tol, abs_tol, max_num_steps,
                               2e-5, 1e-6, 1e-3, 1e-5);
    // FIXME: lag time sensitivity
}

// TEST_F(TorstenCoupledOneCptTest, single_dose_overload) {
//   CoupledOneCptODE f;

//   TORSTEN_ODE_PARAM_OVERLOAD_TEST(pmx_solve_onecpt_rk45, f,
//                                   nPD, time, amt, rate, ii, evid, cmt, addl, ss, parameters,
//                                   biovar, tlag, 1e-10, 1e-10);
//   TORSTEN_ODE_PARAM_OVERLOAD_TEST(pmx_solve_onecpt_bdf, f,
//                                   nPD, time, amt, rate, ii, evid, cmt, addl, ss, parameters,
//                                   biovar, tlag, 1e-10, 1e-10);
// }

// TEST_F(TorstenCoupledOneCptTest, trucated_infusion) {
//   amt[0] = 12000;
//   rate[0] = 12000;
//   ii[0] = 12;
//   addl[0] = 14;

//   CoupledOneCptODE f;

//   double rel_tol = 1e-6, abs_tol = 1e-6;
//   long int max_num_steps = 1e6;
//   MatrixXd x_rk45 = pmx_solve_onecpt_rk45(f, nPD,
//                                           time, amt, rate, ii, evid, cmt, addl, ss,
//                                           parameters, biovar, tlag,
//                                           rel_tol, abs_tol, max_num_steps);

//   rel_tol = 1e-10, abs_tol = 1e-10;
//   max_num_steps = 1e8;
//   MatrixXd x_bdf = pmx_solve_onecpt_bdf(f, nPD,
//                                         time, amt, rate, ii, evid, cmt, addl, ss,
//                                         parameters, biovar, tlag,
//                                         rel_tol, abs_tol, max_num_steps);

//   // Solution from mrgsolve (uses an LSODA integrator)
//   MatrixXd amounts(nt, nOde);
//   amounts << 0.00000,    0.0000,  0.000000e+00,  0.000000e+00,  0.000000e+00,
//              2360.81604,  623.6384, -7.461806e-05, -1.531363e-07, -2.492751e-10,
//              3792.72335, 2098.1390, -5.238333e-04, -2.201387e-06, -7.280973e-09,
//              4661.21904, 4013.1415, -1.562825e-03, -1.006606e-05, -5.066933e-08,
//              5187.98830, 6124.9596, -3.296762e-03, -2.887519e-05, -1.964005e-07,
//              3146.67405, 7667.0022, -5.691112e-03, -6.411818e-05, -5.529458e-07,
//              1908.55429, 8329.8566, -8.448480e-03, -1.198063e-04, -1.267227e-06,
//              1157.59668, 8478.2378, -1.133512e-02, -1.976540e-04, -2.507387e-06,
//              702.11787, 8332.0619, -1.421568e-02, -2.979247e-04, -4.447571e-06,
//              12.85974, 5152.8451, -3.269696e-02, -1.786827e-03, -6.367482e-05;
//   MatrixXd xt = amounts.transpose();
//   torsten::test::test_val(xt, x_rk45, 1.0e-5, 1.0e-10);
//   torsten::test::test_val(xt, x_bdf, 2.0e-5, 1.0e-10);

//   // FIX ME - lag time
//   TORSTEN_ODE_GRAD_THETA_TEST(pmx_solve_onecpt_rk45, f, nPD,
//                               time, amt, rate, ii, evid, cmt, addl, ss,
//                               parameters, biovar, tlag,
//                               rel_tol, abs_tol, max_num_steps,
//                               2e-5, 1e-6, 1e-5, 1e-5);
//   TORSTEN_ODE_GRAD_BIOVAR_TEST(pmx_solve_onecpt_rk45, f, nPD,
//                                time, amt, rate, ii, evid, cmt, addl, ss,
//                                parameters, biovar, tlag,
//                                rel_tol, abs_tol, max_num_steps,
//                                2e-5, 1e-6, 1e-5, 1e-5)
//   TORSTEN_ODE_GRAD_THETA_TEST(pmx_solve_onecpt_bdf, f, nPD,
//                               time, amt, rate, ii, evid, cmt, addl, ss,
//                               parameters, biovar, tlag,
//                               rel_tol, abs_tol, max_num_steps,
//                               2e-5, 1e-6, 1e-3, 1e-5);
//   TORSTEN_ODE_GRAD_BIOVAR_TEST(pmx_solve_onecpt_bdf, f, nPD,
//                                time, amt, rate, ii, evid, cmt, addl, ss,
//                                parameters, biovar, tlag,
//                                rel_tol, abs_tol, max_num_steps,
//                                2e-5, 1e-6, 1e-3, 1e-5)
// }

// TEST_F(TorstenCoupledOneCptTest, ss_bolus) {
//   for (int i = 0; i < nt; i++) time[i] = i * 1.0;
//   amt[0] = 1200;
//   ii[0] = 12;
//   ss[0] = 1;
  
//   CoupledOneCptODE f;

//   double rel_tol = 1e-6;
//   double abs_tol = 1e-6;
//   double max_num_steps = 1e6;
//   MatrixXd x_rk45 = torsten::pmx_solve_onecpt_rk45(f, nPD,
//                                   time, amt, rate, ii, evid, cmt, addl, ss,
//                                   parameters, biovar, tlag,
//                                   rel_tol, abs_tol, max_num_steps);
  
//   rel_tol = 1e-10;
//   abs_tol = 1e-10;
//   max_num_steps = 1e8;
//   MatrixXd x_bdf = torsten::pmx_solve_onecpt_bdf(f, nPD,
//                                 time, amt, rate, ii, evid, cmt, addl, ss,
//                                 parameters, biovar, tlag,
//                                 rel_tol, abs_tol, max_num_steps);
  
//   // Solution from mrgsolve (uses an LSODA integrator)
//   MatrixXd amounts(nt, nOde);
//   amounts << 1.200000e+03,  46.92858, -0.08650993, -0.08766715, -0.08766290,
//              1.624023e+02, 897.86458, -0.08691954, -0.08763443, -0.08766247,
//              2.197877e+01, 791.46490, -0.08761369, -0.08762332, -0.08766135,
//              2.974503e+00, 610.56695, -0.08808584, -0.08763114, -0.08766024,
//              4.025552e-01, 460.96537, -0.08833222, -0.08764989, -0.08765960,
//              5.447993e-02, 346.69437, -0.08840099, -0.08767287, -0.08765965,
//              7.373058e-03, 260.57211, -0.08833521, -0.08769507, -0.08766043,
//              9.978351e-04, 195.81933, -0.08816816, -0.08771281, -0.08766182,
//              1.350418e-04, 147.15449, -0.08792498, -0.08772347, -0.08766361,
//              1.827653e-05, 110.58336, -0.08762457, -0.08772519, -0.08766555;
//   MatrixXd xt = amounts.transpose();
//   torsten::test::test_val(xt, x_rk45, 1.5e-3, 1.0e-10);
//   torsten::test::test_val(xt, x_bdf, 1.5e-3, 1.0e-10);
  
//   TORSTEN_ODE_GRAD_THETA_TEST(pmx_solve_onecpt_rk45, f, nPD,
//                               time, amt, rate, ii, evid, cmt, addl, ss,
//                               parameters, biovar, tlag,
//                               rel_tol, abs_tol, max_num_steps,
//                               2e-5, 1e-6, 1e-2, 1e-5);
//   TORSTEN_ODE_GRAD_BIOVAR_TEST(pmx_solve_onecpt_rk45, f, nPD,
//                                time, amt, rate, ii, evid, cmt, addl, ss,
//                                parameters, biovar, tlag,
//                                rel_tol, abs_tol, max_num_steps,
//                                2e-5, 1e-6, 1e-2, 1e-5)
//   TORSTEN_ODE_GRAD_THETA_TEST(pmx_solve_onecpt_bdf, f, nPD,
//                               time, amt, rate, ii, evid, cmt, addl, ss,
//                               parameters, biovar, tlag,
//                               rel_tol, abs_tol, max_num_steps,
//                               2e-5, 1e-6, 1e-2, 1e-5);
//   TORSTEN_ODE_GRAD_BIOVAR_TEST(pmx_solve_onecpt_bdf, f, nPD,
//                                time, amt, rate, ii, evid, cmt, addl, ss,
//                                parameters, biovar, tlag,
//                                rel_tol, abs_tol, max_num_steps,
//                                2e-5, 1e-6, 1e-2, 1e-5)
// }

// TEST_F(TorstenCoupledOneCptTest, ss_truncated_infusion) {
//   for (int i = 0; i < nt; i++) time[i] = i * 1.0;
//   amt[0] = 1200;
//   rate[0] = 150;
//   ii[0] = 12;
//   ss[0] = 1;
  
//   CoupledOneCptODE f;

//   double rel_tol = 1e-6, abs_tol = 1e-6;
//   long int max_num_steps = 1e6;
//   MatrixXd x_rk45 = torsten::pmx_solve_onecpt_rk45(f, nPD,
//                                   time, amt, rate, ii, evid, cmt, addl, ss,
//                                   parameters, biovar, tlag,
//                                   rel_tol, abs_tol, max_num_steps);
  
//   rel_tol = 1e-10, abs_tol = 1e-10;
//   max_num_steps = 1e8;
//   MatrixXd x_bdf = torsten::pmx_solve_onecpt_bdf(f, nPD,
//                                 time, amt, rate, ii, evid, cmt, addl, ss,
//                                 parameters, biovar, tlag,
//                                 rel_tol, abs_tol, max_num_steps);
  
//   // Solution from mrgsolve (uses an LSODA integrator)
//   MatrixXd amounts(nt, nOde);
//   amounts << 0.0251597, 181.3172, -0.08776322, -0.08770075, -0.08766391,
//              64.8532585, 212.8358, -0.08754360, -0.08769911, -0.08766506,
//              73.6267877, 283.1219, -0.08740574, -0.08769178, -0.08766603,
//              74.8141559, 342.2470, -0.08735655, -0.08768178, -0.08766669,
//              74.9748488, 387.5317, -0.08737769, -0.08767172, -0.08766700,
//              74.9965962, 421.6776, -0.08745217, -0.08766351, -0.08766701,
//              74.9995393, 447.3531, -0.08756680, -0.08765858, -0.08766681,
//              74.9999377, 466.6498, -0.08771162, -0.08765792, -0.08766653,
//              74.9999915, 481.1511, -0.08787912, -0.08766221, -0.08766631,
//              10.1501451, 415.4866, -0.08802304, -0.08767157, -0.08766632;
//   MatrixXd xt = amounts.transpose();
//   torsten::test::test_val(xt, x_rk45, 1e-3, 1.0e-10);
//   torsten::test::test_val(xt, x_bdf, 1e-3, 1.0e-10);
  
//   TORSTEN_ODE_GRAD_THETA_TEST(pmx_solve_onecpt_rk45, f, nPD,
//                               time, amt, rate, ii, evid, cmt, addl, ss,
//                               parameters, biovar, tlag,
//                               rel_tol, abs_tol, max_num_steps,
//                               2e-5, 1e-6, 1e-5, 1e-5);
//   // FIXME: steady state exception
//   // TORSTEN_ODE_GRAD_BIOVAR_TEST(pmx_solve_onecpt_rk45, f, nPD,
//   //                              time, amt, rate, ii, evid, cmt, addl, ss,
//   //                              parameters, biovar, tlag,
//   //                              rel_tol, abs_tol, max_num_steps,
//   //                              2e-5, 1e-6, 1e-5, 1e-5)
//   TORSTEN_ODE_GRAD_THETA_TEST(pmx_solve_onecpt_bdf, f, nPD,
//                               time, amt, rate, ii, evid, cmt, addl, ss,
//                               parameters, biovar, tlag,
//                               rel_tol, abs_tol, max_num_steps,
//                               2e-5, 1e-6, 1e-3, 1e-5);
//   // FIXME: steady state exception
//   // TORSTEN_ODE_GRAD_BIOVAR_TEST(pmx_solve_onecpt_bdf, f, nPD,
//   //                              time, amt, rate, ii, evid, cmt, addl, ss,
//   //                              parameters, biovar, tlag,
//   //                              rel_tol, abs_tol, max_num_steps,
//   //                              2e-5, 1e-6, 1e-3, 1e-5)
// }

// TEST_F(TorstenCoupledOneCptTest, ss_const_infusion) {
//   for (int i = 0; i < nt; i++) time[i] = i * 1.0;
//   amt[0] = 0.0;
//   rate[0] = 150;
//   ss[0] = 1;
  
//   CoupledOneCptODE f;

//   double rel_tol_rk = 1e-6, abs_tol_rk = 1e-6;
//   long int max_num_steps_rk = 1e6;
//   MatrixXd x_rk45 = torsten::pmx_solve_onecpt_rk45(f, nPD,
//                                   time, amt, rate, ii, evid, cmt, addl, ss,
//                                   parameters, biovar, tlag,
//                                   rel_tol_rk, abs_tol_rk, max_num_steps_rk);
  
//   double rel_tol_bdf = 1e-10, abs_tol_bdf = 1e-10;
//   long int max_num_steps_bdf = 1e8;
//   MatrixXd x_bdf = torsten::pmx_solve_onecpt_bdf(f, nPD,
//                                 time, amt, rate, ii, evid, cmt, addl, ss,
//                                 parameters, biovar, tlag,
//                                 rel_tol_bdf, abs_tol_bdf, max_num_steps_bdf);
  
//   // can't do constant rate in mrgsolve. Comparing to result obtained
//   // with generalOdeModel, as a provisional test.
//   MatrixXd x = torsten::pmx_solve_rk45(f, nOde,
//                                        time, amt, rate, ii, evid, cmt, addl, ss,
//                                        parameters, biovar, tlag,
//                                        rel_tol_rk, abs_tol_rk, max_num_steps_rk, nullptr);
//   torsten::test::test_val(x, x_rk45, 1.5e-2, 1.0e-5);
//   torsten::test::test_val(x, x_bdf, 1.5e-2, 1.0e-5);
  
//   // Test AutoDiff against FiniteDiff
//   // diff2 obtained empirically.
//   TORSTEN_ODE_GRAD_THETA_TEST(pmx_solve_onecpt_rk45, f, nPD,
//                               time, amt, rate, ii, evid, cmt, addl, ss,
//                               parameters, biovar, tlag,
//                               rel_tol_rk, abs_tol_rk, max_num_steps_rk,
//                               2e-5, 1e-6, 1e-5, 1e-5);
//   // FIXME: steady state exception
//   // TORSTEN_ODE_GRAD_BIOVAR_TEST(pmx_solve_onecpt_rk45, f, nPD,
//   //                              time, amt, rate, ii, evid, cmt, addl, ss,
//   //                              parameters, biovar, tlag,
//   //                              rel_tol, abs_tol, max_num_steps,
//   //                              2e-5, 1e-6, 1e-5, 1e-5)
//   TORSTEN_ODE_GRAD_THETA_TEST(pmx_solve_onecpt_bdf, f, nPD,
//                               time, amt, rate, ii, evid, cmt, addl, ss,
//                               parameters, biovar, tlag,
//                               rel_tol_bdf, abs_tol_bdf, max_num_steps_bdf,
//                               2e-5, 1e-6, 1e-3, 1e-5);
//   // FIXME: steady state exception
//   // TORSTEN_ODE_GRAD_BIOVAR_TEST(pmx_solve_onecpt_bdf, f, nPD,
//   //                              time, amt, rate, ii, evid, cmt, addl, ss,
//   //                              parameters, biovar, tlag,
//   //                              rel_tol, abs_tol, max_num_steps,
//   //                              2e-5, 1e-6, 1e-3, 1e-5)
// }
