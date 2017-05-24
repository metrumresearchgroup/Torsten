#include <stan/math/rev/mat.hpp>  // FIX ME - more specific
#include <stan/math/torsten/PKModel/Pred/fTwoCpt.hpp>
#include <gtest/gtest.h>
#include <ostream>
#include <fstream>
#include <vector>


struct OneCpt_functor {
  /**
   * parms contains all parameters.
   * x contains all initial states.
   */
  template <typename T0, typename T1, typename T2, typename T3>
  inline
  std::vector<typename boost::math::tools::promote_args<T0, T1, T2, T3>::type>
  operator()(const T0& t,
             const std::vector<T1>& y,
             const std::vector<T2>& parms,
             const std::vector<T3>& x_r,
             const std::vector<int>& x_i,
             std::ostream* pstream__) const {
    using Eigen::Matrix;
    using Eigen::Dynamic;
    typedef typename boost::math::tools::promote_args<T0, T1, T2, T3>::type scalar;
    T2
      CL = parms[0],
      VC = parms[1],
      ka = parms[2];


    std::vector<scalar> dydt(2);
    dydt[0] = -ka * y[0];
    dydt[1] = ka * y[0] - CL / VC * y[1];

    return dydt;
  }
};

struct OneCpt_mix_functor {
  /**
   * parms contains all parameters and the y1 initial state.
   * x contains the y2 initial state.
   */
  template <typename T0, typename T1, typename T2, typename T3>
  inline
  std::vector<typename boost::math::tools::promote_args<T0, T1, T2, T3>::type>
  operator()(const T0& t,
             const std::vector<T1>& y,
             const std::vector<T2>& parms,
             const std::vector<T3>& x_r,
             const std::vector<int>& x_i,
             std::ostream* pstream__) const {
    using Eigen::Matrix;
    using Eigen::Dynamic;
    typedef typename boost::math::tools::promote_args<T0, T1, T2, T3>::type scalar;
    T2
      CL = parms[0],
      VC = parms[1],
      ka = parms[2],
      y1_0 = parms[3];

    scalar y1 = y1_0 * exp(-ka * t);

    std::vector<scalar> dydt(1);
    dydt[0] = ka * y1 - CL / VC * y[0];

    return dydt;
  }
};



TEST(mix_solver, dbl) {
  /**
   * TEST 1: -- everything is passed as doubles
   * No sensitivities need to be calculated.
   */
  using std::vector;
  using stan::math::integrate_ode_rk45;

  typedef double scalar;

  vector<double> init(2, 0);
  init[0] = 1000;  // bolus dose in the gut
  double t0 = 0;
  vector<double> t(1);
  t[0] = 0.25;

  vector<scalar> parms(3);
  parms[0] = 10;  // CL
  parms[1] = 80;  // VC
  parms[2] = 1.2;  // ka

  vector<double> x_r(1, 0);
  vector<int> x_i(1, 0);

  // need to make sure we pass dt
  vector<double> dt(1);
  dt[0] = t[0] - t0;

  // results from mrgsolve (use to make sure the solver works)
  vector<double> mrgSolution(2);
  mrgSolution[0] = 740.8182;
  mrgSolution[1] = 254.97490;

  // Save CPU times inside output file
  std::ofstream myfile;
  myfile.open ("test/unit/math/torsten/mixSolver/msSimple_dd.csv");

  int N = 1000;  // number of simulations
  for (int i = 0; i < N; i++) myfile << i << ", ";
  myfile << "0\n";

  // Use full numerical method
  for (int i = 0; i < N; i++) {
    clock_t start = clock();

    vector<vector<scalar> >
      temp = integrate_ode_rk45(OneCpt_functor(), init, 0, dt, parms, x_r, x_i);

    clock_t end = clock();

    myfile << (float)(end - start) / CLOCKS_PER_SEC << ", ";

    // Check accuracy of result
    for (int i = 0; i < 2; i++) {
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

    int nStates = 1;
    int nBase = 1;
    vector<double> init_mix(nStates);
    for (int i = 0; i < nStates; i++) init_mix[i] = init[nBase + i];

    // clock_t start2 = clock();

    vector<vector<scalar> >
      temp_PD = integrate_ode_rk45(OneCpt_mix_functor(), init_mix, 0, dt, parms,
                                   x_r, x_i);

    // clock_t end2 = clock();

    scalar tempBase = init[0] * exp(-parms[2] * dt[0]);

    vector<vector<scalar> > temp_mix(1);
    temp_mix[0].resize(nBase + nStates);
    for (int i = 0 ; i < nBase; i++) temp_mix[0][i] = tempBase;
    for (int i = 0 ; i < nStates; i++) temp_mix[0][nBase + i] = temp_PD[0][i];

    clock_t end2 = clock();
    myfile << (float)(end2 - start2) / CLOCKS_PER_SEC << ", ";

    for (int i = 0; i < 2; i++) {
      // FIX ME - what value should I use for the relative error?
      double rel_err  = std::max(std::abs(1e-3 * mrgSolution[i]), 5e-13);
      EXPECT_NEAR(mrgSolution[i], temp_mix[0][i], rel_err);
    }
  }

  myfile << "0\n";
  myfile.close();
}

TEST(mix_solver,var) {
  /**
   * TEST 1: -- everything is passed as doubles
   * No sensitivities need to be calculated.
   */
  using std::vector;
  using stan::math::integrate_ode_rk45;
  using stan::math::var;

  typedef var scalar;

  vector<scalar> init(2, 0);
  init[0] = 1000;  // bolus dose in the gut
  double t0 = 0;
  vector<double> t(1);
  t[0] = 0.25;

  vector<scalar> parms(3);
  parms[0] = 10;  // CL
  parms[1] = 80;  // VC
  parms[2] = 1.2;  // ka

  vector<double> x_r(1, 0);
  vector<int> x_i(1, 0);

  // need to make sure we pass dt
  vector<double> dt(1);
  dt[0] = t[0] - t0;

  // results from mrgsolve (use to make sure the solver works)
  vector<double> mrgSolution(2);
  mrgSolution[0] = 740.8182;
  mrgSolution[1] = 254.97490;

  // Save CPU times inside output file
  std::ofstream myfile;
  myfile.open ("test/unit/math/torsten/mixSolver/msSimple_vv.csv");

  int N = 1000;  // number of simulations
  for (int i = 0; i < N; i++) myfile << i << ", ";
  myfile << "0\n";

  // Use full numerical method
  for (int i = 0; i < N; i++) {
    clock_t start = clock();

    vector<vector<scalar> >
      temp = integrate_ode_rk45(OneCpt_functor(), init, 0, dt, parms, x_r, x_i);

    clock_t end = clock();

    myfile << (float)(end - start) / CLOCKS_PER_SEC << ", ";

    // Check accuracy of result
    for (int i = 0; i < 2; i++) {
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

    int nStates = 1;
    int nBase = 1;
    vector<scalar> init_mix(nStates);
    for (int i = 0; i < nStates; i++) init_mix[i] = init[nBase + i];

    // clock_t start2 = clock();

    vector<vector<scalar> >
      temp_PD = integrate_ode_rk45(OneCpt_mix_functor(), init_mix, 0, dt, parms,
                                   x_r, x_i);

    // clock_t end2 = clock();

    scalar tempBase = init[0] * exp(-parms[2] * dt[0]);

    vector<vector<scalar> > temp_mix(1);
    temp_mix[0].resize(nBase + nStates);
    for (int i = 0 ; i < nBase; i++) temp_mix[0][i] = tempBase;
    for (int i = 0 ; i < nStates; i++) temp_mix[0][nBase + i] = temp_PD[0][i];

    clock_t end2 = clock();
    myfile << (float)(end2 - start2) / CLOCKS_PER_SEC << ", ";

    for (int i = 0; i < 2; i++) {
      // FIX ME - what value should I use for the relative error?
      double rel_err  = std::max(std::abs(1e-3 * mrgSolution[i]), 5e-13);
      EXPECT_NEAR(mrgSolution[i], temp_mix[0][i].val(), rel_err);
    }
  }

  myfile << "0\n";
  myfile.close();
}
