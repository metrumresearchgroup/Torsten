#include <stan/math/rev/mat.hpp>  // FIX ME - more specific
#include <stan/math/torsten/PKModel/Pred/fTwoCpt.hpp>
#include <gtest/gtest.h>
#include <ostream>
#include <fstream>
#include <vector>


struct FK_functor {
  /**
   * parms contains both the PK and the PD parameters.
   * x contains both the PK and the PD states.
   */
  template <typename T0, typename T1, typename T2, typename T3>
  inline
  std::vector<typename boost::math::tools::promote_args<T0, T1, T2, T3>::type>
  operator()(const T0& t,
             const std::vector<T1>& x,
             const std::vector<T2>& parms,
             const std::vector<T3>& x_r,
             const std::vector<int>& x_i,
             std::ostream* pstream__) const {
    using Eigen::Matrix;
    using Eigen::Dynamic;
    typedef typename boost::math::tools::promote_args<T0, T1, T2, T3>::type scalar;

    // PK variables
    scalar
      CL = parms[0],
      Q = parms[1],
      VC = parms[2],
      VP = parms[3],
      ka = parms[4],
      k10 = CL / VC,
      k12 = Q / VC,
      k21 = Q / VP;

    // PD variables
    scalar
      MTT = parms[5],
      circ0 = parms[6],
      alpha = parms[7],
      gamma = parms[8],
      ktr = 4 / MTT,
      prol = x[3] + circ0,
      transit1 = x[4] + circ0,
      transit2 = x[5] + circ0,
      transit3 = x[6] + circ0,
      circ = x[7] + circ0;

    std::vector<scalar> dxdt(8);
    dxdt[0] = -ka * x[0];
    dxdt[1] = ka * x[0] - (k10 + k12) * x[1] + k21 * x[2];
    dxdt[2] = k12 * x[1] - k21 * x[2];

    scalar conc = x[1] / VC;
    scalar Edrug = alpha * conc;

    dxdt[3] = ktr * prol * ((1 - Edrug) * pow((circ0 / circ), gamma) - 1);
    dxdt[4] = ktr * (prol - transit1);
    dxdt[5] = ktr * (transit1 - transit2);
    dxdt[6] = ktr * (transit2 - transit3);
    dxdt[7] = ktr * (transit3 - circ);

    return dxdt;
  }
};

struct FK_mix_functor {
  /**
   * parms contains both the PK and the PD parameters, and the
   * PK states.
   * x contains the PD states.
   */
  template <typename T0, typename T1, typename T2, typename T3>
  inline
  std::vector<typename boost::math::tools::promote_args<T0, T1, T2, T3>::type>
  operator()(const T0& t,
             const std::vector<T1>& x,
             const std::vector<T2>& parms,
             const std::vector<T3>& rate,
             const std::vector<int>& dummy, std::ostream* pstream__) const {
    using std::vector;
    using Eigen::Matrix;
    using Eigen::Dynamic;
    typedef typename boost::math::tools::promote_args<T0, T1, T2, T3>::type scalar;
    // PK variables
    scalar VC = parms[2];

    // PD variables
    scalar
      MTT = parms[5],
      circ0 = parms[6],
      alpha = parms[7],
      gamma = parms[8],
      ktr = 4 / MTT,
      prol = x[0] + circ0,
      transit1 = x[1] + circ0,
      transit2 = x[2] + circ0,
      transit3 = x[3] + circ0,
      circ = x[4] + circ0;

    vector<T2> initPK(3);
    // Matrix<T2, 1, Dynamic> initPK(3);
    initPK[0] = parms[9];
    initPK[1] = parms[10];
    initPK[2] = parms[11];

    vector<scalar> parmsPK(5);
    // Matrix<T2, 1, Dynamic> parmsPK(5);
    for (size_t i = 0; i < parmsPK.size(); i++)
      parmsPK[i] = parms[i];

    vector<scalar> predPK = fTwoCpt(t, parmsPK, initPK, rate);
    // Matrix<T2, 1, Dynamic> predPK = predTwoCpt(t, parmsPK, initPK);

    scalar conc = predPK[1] / VC;
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

TEST(mixed_solver, dbl) {
  /**
   * TEST 1: -- everything is passed as doubles
   * No sensitivities need to be calculated.
   */
  using std::vector;
  using stan::math::var;
  using stan::math::integrate_ode_rk45;
  using stan::math::to_vector;

  typedef double scalar;

  vector<double> init(8, 0);
  init[0] = 80000;  // bolus dose in the gut
  double t0 = 0;
  vector<double> t(1);
  t[0] = 1;

  vector<scalar> parms(9);
  parms[0] = 10;  // CL
  parms[1] = 15;  // Q
  parms[2] = 35;  // VC
  parms[3] = 105;  // VP
  parms[4] = 2;  // ka
  parms[5] = 125;  // MTT
  parms[6] = 5;  // Circ0
  parms[7] = 3e-4;  // alpha
  parms[8] = 0.17;  // gamma

  vector<double> x_r(1, 0);
  vector<int> x_i(1, 0);

   vector<double> rate(8, 0);  // required for fTwoCpt

  // need to make sure we pass dt (for fTwoCpt)
  vector<double> dt(1);
  dt[0] = t[0] - t0;

  // results from mrgsolve (use to make sure the solver works)
  vector<double> mrgSolution(8);
  mrgSolution[0] = 10826.82;
  mrgSolution[1] = 44779.91;
  mrgSolution[2] = 14296.11;
  mrgSolution[3] = -0.04823222;
  mrgSolution[4] = -0.0006261043;
  mrgSolution[5] = -5.620476e-06;
  mrgSolution[6] = -3.883612e-08;
  mrgSolution[7] = -2.188286e-10;

  // Save CPU times inside output file
  std::ofstream myfile;
  myfile.open ("test/unit/math/torsten/mixSolver/mixSolverResult_dd.csv");

  int N = 1000;  // number of simulations
  for (int i = 0; i < N; i++) myfile << i << ", ";
  myfile << "0\n";

  // Use full numerical method
  for (int i = 0; i < N; i++) {
    clock_t start = clock();

    vector<vector<scalar> >
      temp = integrate_ode_rk45(FK_functor(), init, 0, dt, parms, x_r, x_i);

    clock_t end = clock();

    myfile << (float)(end - start) / CLOCKS_PER_SEC << ", ";

    // Check accuracy of result
    for (int i = 0; i < 8; i++) {
      // FIX ME - what value should I use for the relative error?
      double rel_err  = std::max(std::abs(1e-3 * mrgSolution[i]), 5e-13);
      EXPECT_NEAR(mrgSolution[i], temp[0][i], rel_err);
    }
  }
  myfile << "0\n";

  // Use mix solver
  for (int i = 0; i < N; i++) {
    clock_t start2 = clock();

    parms.push_back(init[0]);
    parms.push_back(init[1]);
    parms.push_back(init[2]);

    int nStates = 5;
    vector<double> init_mix(nStates);
    for (int i = 0; i < nStates; i++) init_mix[i] = init[3 + i];

    vector<vector<scalar> >
      temp_PD = integrate_ode_rk45(FK_mix_functor(), init_mix, 0, dt, parms,
                                   rate, x_i);

    int nParmsPK = 5;
    vector<scalar> parmsPK(nParmsPK);
    for (int i = 0; i < nParmsPK; i++) parmsPK[i] = parms[i];

    int nPK = 3;
    vector<scalar> initPK(nPK);
    for (int i = 0; i < nPK; i++) initPK[i] = init[i];

    vector<scalar> tempPK = fTwoCpt(dt[0], parmsPK, initPK, rate);

    vector<vector<scalar> > temp_mix(1);
    temp_mix[0].resize(8);
    for (int i = 0 ; i < 3; i++) temp_mix[0][i] = tempPK[i];
    for (int i = 0 ; i < 5; i++) temp_mix[0][3 + i] = temp_PD[0][i];

    clock_t end2 = clock();
    myfile << (float)(end2 - start2) / CLOCKS_PER_SEC << ", ";

    for (int i = 0; i < 8; i++) {
      // FIX ME - what value should I use for the relative error?
      double rel_err  = std::max(std::abs(1e-3 * mrgSolution[i]), 5e-13);
      EXPECT_NEAR(mrgSolution[i], temp_mix[0][i], rel_err);
    }
  }

  myfile << "0\n";
  myfile.close();
}

TEST(mixed_solver, var) {
  /**
   * TEST 1: -- everything is passed as var
   * Sensitivities with respect to parameters and
   * initial estimated need  to be calculated.
   */
  using std::vector;
  using stan::math::var;
  using stan::math::integrate_ode_rk45;
  using stan::math::to_vector;

  typedef var scalar;

  vector<scalar> init(8, 0);
  init[0] = 80000;  // bolus dose in the gut
  double t0 = 0;
  vector<double> t(1);
  t[0] = 1;

  vector<scalar> parms(9);
  parms[0] = 10;  // CL
  parms[1] = 15;  // Q
  parms[2] = 35;  // VC
  parms[3] = 105;  // VP
  parms[4] = 2;  // ka
  parms[5] = 125;  // MTT
  parms[6] = 5;  // Circ0
  parms[7] = 3e-4;  // alpha
  parms[8] = 0.17;  // gamma

  vector<double> x_r(1, 0);
  vector<int> x_i(1, 0);

  // need to make sure we pass dt (for fTwoCpt)
  vector<double> dt(1);
  dt[0] = t[0] - t0;

  // results from mrgsolve (use to make sure the solver works)
  vector<double> mrgSolution(8);
  mrgSolution[0] = 10826.82;
  mrgSolution[1] = 44779.91;
  mrgSolution[2] = 14296.11;
  mrgSolution[3] = -0.04823222;
  mrgSolution[4] = -0.0006261043;
  mrgSolution[5] = -5.620476e-06;
  mrgSolution[6] = -3.883612e-08;
  mrgSolution[7] = -2.188286e-10;

  // Save CPU times inside output file
  std::ofstream myfile;
  myfile.open ("test/unit/math/torsten/mixSolver/mixSolverResult_vv.csv");

  int N = 1000;  // number of simulations
  for (int i = 0; i < N; i++) myfile << i << ", ";
  myfile << "0\n";

  // Use full numerical method
  for (int i = 0; i < N; i++) {
    clock_t start = clock();

    vector<vector<scalar> >
      temp = integrate_ode_rk45(FK_functor(), init, 0, dt, parms, x_r, x_i);

    clock_t end = clock();

    myfile << (float)(end - start) / CLOCKS_PER_SEC << ", ";

    for (int i = 0; i < 8; i++) {
      // FIX ME - what value should I use for the relative error?
      double rel_err  = std::max(std::abs(1e-3 * mrgSolution[i]), 5e-13);
      EXPECT_NEAR(mrgSolution[i], temp[0][i].val(), rel_err);
    }
  }
  myfile << "0\n";

  // Use mix solver
  for (int i = 0; i < N; i++) {
    clock_t start2 = clock();

    parms.push_back(init[0]);
    parms.push_back(init[1]);
    parms.push_back(init[2]);

    vector<double> rate(8, 0);  // required for fTwoCpt

    int nStates = 5;
    vector<scalar> init_mix(nStates);
    for (int i = 0; i < nStates; i++) init_mix[i] = init[3 + i];

    vector<vector<scalar> >
      temp_PD = integrate_ode_rk45(FK_mix_functor(), init_mix, 0, dt, parms,
                                 rate, x_i);

    int nParmsPK = 5;
    vector<scalar> parmsPK(nParmsPK);
    for (int i = 0; i < nParmsPK; i++) parmsPK[i] = parms[i];

    int nPK = 3;
    vector<scalar> initPK(nPK);
    for (int i = 0; i < nPK; i++) initPK[i] = init[i];

    vector<scalar> tempPK = fTwoCpt(dt[0], parmsPK, initPK, rate);

    vector<vector<scalar> > temp_mix(1);
    temp_mix[0].resize(8);
    for (int i = 0 ; i < 3; i++) temp_mix[0][i] = tempPK[i];
    for (int i = 0 ; i < 5; i++) temp_mix[0][3 + i] = temp_PD[0][i];

    clock_t end2 = clock();
    myfile << (float)(end2 - start2) / CLOCKS_PER_SEC << ", ";

    for (int i = 0; i < 8; i++) {
      // FIX ME - what value should I use for the relative error?
      double rel_err  = std::max(std::abs(1e-3 * mrgSolution[i]), 5e-13);
      EXPECT_NEAR(mrgSolution[i], temp_mix[0][i].val(), rel_err);
    }
  }

  myfile << "0\n";
  myfile.close();
}
