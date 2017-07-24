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
  double rel_tol = 1e-6, abs_tol = 1e-6;
  long int max_num_steps = 1e6;
  Matrix<double, Dynamic, Dynamic>
    x_rk45 = mixOde1CptModel_rk45(feedbackODE(), nPD,
                                  time, amt, rate, ii, evid, cmt, addl, ss,
                                  parameters, biovar, tlag,
                                  0,
                                  rel_tol, abs_tol, max_num_steps);

  rel_tol = 1e-10, abs_tol = 1e-10;
  max_num_steps = 1e8;
  Matrix<double, Dynamic, Dynamic>
    x_bdf = mixOde1CptModel_bdf(feedbackODE(), nPD,
                                time, amt, rate, ii, evid, cmt, addl, ss,
                                parameters, biovar, tlag,
                                0,
                                rel_tol, abs_tol, max_num_steps);

  // Solution from mrgsolve (uses an LSODA integrator)
  Matrix<double, Dynamic, Dynamic> amounts(nEvents, nOde);
  amounts << 1.200000e+03,  46.92858, -0.08650993, -0.08766715, -0.08766290,
             1.624023e+02, 897.86458, -0.08691954, -0.08763443, -0.08766247,
             2.197877e+01, 791.46490, -0.08761369, -0.08762332, -0.08766135,
             2.974503e+00, 610.56695, -0.08808584, -0.08763114, -0.08766024,
             4.025552e-01, 460.96537, -0.08833222, -0.08764989, -0.08765960,
             5.447993e-02, 346.69437, -0.08840099, -0.08767287, -0.08765965,
             7.373058e-03, 260.57211, -0.08833521, -0.08769507, -0.08766043,
             9.978351e-04, 195.81933, -0.08816816, -0.08771281, -0.08766182,
             1.350418e-04, 147.15449, -0.08792498, -0.08772347, -0.08766361,
             1.827653e-05, 110.58336, -0.08762457, -0.08772519, -0.08766555;

  double rel_err = 1.1e-3;
  expect_near_matrix_eq(amounts, x_rk45, rel_err);

  expect_near_matrix_eq(amounts, x_bdf, rel_err);


  // Test AutoDiff against FiniteDiff
  // diff2 obtained empirically
  double diff = 1e-8, diff2 = 6e-2;
  test_mixOdeCptModel(feedbackODE(), nPD,
                      time, amt, rate, ii, evid, cmt, addl, ss,
                      parameters, biovar, tlag,
                      rel_tol, abs_tol, max_num_steps, diff, diff2, "1_rk45");

  test_mixOdeCptModel(feedbackODE(), nPD,
                      time, amt, rate, ii, evid, cmt, addl, ss,
                      parameters, biovar, tlag,
                      rel_tol, abs_tol, max_num_steps, diff, diff2, "1_bdf");
}

TEST(Torsten, mixOde1Cpt_SS_infusion) {
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
  double rel_tol = 1e-6, abs_tol = 1e-6;
  long int max_num_steps = 1e6;
  Matrix<double, Dynamic, Dynamic>
    x_rk45 = mixOde1CptModel_rk45(feedbackODE(), nPD,
                                  time, amt, rate, ii, evid, cmt, addl, ss,
                                  parameters, biovar, tlag,
                                  0,
                                  rel_tol, abs_tol, max_num_steps);
  
  // std::cout << x_rk45 << std::endl;
  
  rel_tol = 1e-10, abs_tol = 1e-10;
  max_num_steps = 1e8;
  Matrix<double, Dynamic, Dynamic>
    x_bdf = mixOde1CptModel_bdf(feedbackODE(), nPD,
                                time, amt, rate, ii, evid, cmt, addl, ss,
                                parameters, biovar, tlag,
                                0,
                                rel_tol, abs_tol, max_num_steps);
  
  // Solution from mrgsolve (uses an LSODA integrator)
  Matrix<double, Dynamic, Dynamic> amounts(nEvents, nOde);
  amounts << 0.0251597, 181.3172, -0.08776322, -0.08770075, -0.08766391,
             64.8532585, 212.8358, -0.08754360, -0.08769911, -0.08766506,
             73.6267877, 283.1219, -0.08740574, -0.08769178, -0.08766603,
             74.8141559, 342.2470, -0.08735655, -0.08768178, -0.08766669,
             74.9748488, 387.5317, -0.08737769, -0.08767172, -0.08766700,
             74.9965962, 421.6776, -0.08745217, -0.08766351, -0.08766701,
             74.9995393, 447.3531, -0.08756680, -0.08765858, -0.08766681,
             74.9999377, 466.6498, -0.08771162, -0.08765792, -0.08766653,
             74.9999915, 481.1511, -0.08787912, -0.08766221, -0.08766631,
             10.1501451, 415.4866, -0.08802304, -0.08767157, -0.08766632;
  

  // the relative error is 1.5% for this problem, which is
  // comparitively high.
  double rel_err = 1.5e-2;
  expect_near_matrix_eq(amounts, x_rk45, rel_err);
  expect_near_matrix_eq(amounts, x_bdf, rel_err);
  
  // Test AutoDiff against FiniteDiff
  // diff2 obtained empirically
  // double diff = 1e-8, diff2 = 6e-2;
  // test_mixOdeCptModel(feedbackODE(), nPD,
  //                     time, amt, rate, ii, evid, cmt, addl, ss,
  //                     parameters, biovar, tlag,
  //                     rel_tol, abs_tol, max_num_steps, diff, diff2, "1_rk45");
  // 
  // test_mixOdeCptModel(feedbackODE(), nPD,
  //                     time, amt, rate, ii, evid, cmt, addl, ss,
  //                     parameters, biovar, tlag,
  //                     rel_tol, abs_tol, max_num_steps, diff, diff2, "1_bdf");
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
      VC = theta[1],
      ka = theta[2],
      Mtt = theta[3],
      circ0 = theta[4],
      alpha = theta[5],
      gamma = theta[6],
      ktr = 4 / Mtt,
      prol = x[2] + circ0,
      transit = x[3] + circ0,
      circ = x[4] + circ0,
      Edrug;

    std::vector<scalar> dxdt(5);
    dxdt[0] = - ka * x[0];
    dxdt[1] = ka * x[0] - CL / VC * x[1];

    Edrug = alpha * x[1] / VC;

    dxdt[2] = ktr * prol * ((1 - Edrug) * pow(circ0 / circ, gamma) - 1);
    dxdt[3] = ktr * (prol - transit);
    dxdt[4] = ktr * (transit - circ);

    return dxdt;
  }
};


TEST(Torsten, generalOdeCpt_SS_infusion) {
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

double rel_tol = 1e-6, abs_tol = 1e-6;
long int max_num_steps = 1e6;
Matrix<double, Dynamic, Dynamic>
x_rk45 = generalOdeModel_rk45(fullODE(), nOde,
                              time, amt, rate, ii, evid, cmt, addl, ss,
parameters, biovar, tlag,
0,
rel_tol, abs_tol, max_num_steps);

// std::cout << x_rk45 << std::endl;
}
