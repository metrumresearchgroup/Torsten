#include <stan/math/rev/mat.hpp>  // FIX ME - more specific
#include <stan/math/torsten/fTwoCpt.hpp>
#include <test/unit/math/torsten/mixSolver/util.hpp>
// #include <test/unit/math/rev/arr/functor/util.hpp>
#include <gtest/gtest.h>


struct FK_functor {
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

struct FK_mixed_functor {
    template <typename T0, typename T1, typename T2, typename T3>
  inline
  std::vector<typename boost::math::tools::promote_args<T0, T1, T2, T3>::type>
  operator()(const T0& t,
             const std::vector<T1>& x,
             const std::vector<T2>& parms,
             const std::vector<T3>& rate,
             const std::vector<int>& dummy, std::ostream* pstream__) const {
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

    Matrix<scalar, Dynamic, 1> initPK(3);
    initPK(0) = parms[9];
    initPK(1) = parms[10];
    initPK(2) = parms[11];

    Matrix<scalar, Dynamic, 1> parmsPK(5);
    for (int i = 0; i < parmsPK.size(); i++)
      parmsPK(i) = parms[i];


    Matrix<scalar, Dynamic, 1> predPK = fTwoCpt(t, parmsPK, initPK, rate);

    scalar conc = predPK(1) / VC;
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
  using std::vector;
  using stan::math::var;
  using stan::math::integrate_ode_rk45;
  using stan::math::to_vector;

  vector<double> init(8, 0);
  init[0] = 80000;  // bolus dose in the gut
  double t0 = 0;
  vector<double> t(1);
  t[0] = 1;

  vector<double> parms(9);
  parms[0] = 10;
  parms[1] = 15;
  parms[2] = 35;
  parms[3] = 105;
  parms[4] = 2;
  parms[5] = 125;
  parms[6] = 5;
  parms[7] = 3e-4;
  parms[8] = 0.17;

  vector<double> x_r(1, 0);
  vector<int> x_i(1, 0);

  // need to make sure we pass dt (for fTwoCpt)
  vector<double> dt(1);
  dt[0] = t[0] - t0;

  // Use full numerical method
  clock_t start = clock();

  vector<vector<double> >
    temp = integrate_ode_rk45(FK_functor(), init, 0, dt, parms, x_r, x_i);

  clock_t end = clock();
  std::cout << "(dbl) CPU time for num solver: "
            <<  (float)(end - start) / CLOCKS_PER_SEC
            << " s"  << std::endl;

  // Use mix solver
  clock_t start2 = clock();

  parms.push_back(init[0]);
  parms.push_back(init[1]);
  parms.push_back(init[2]);

  vector<double> rate(8, 0);  // required for fTwoCpt

  int nStates = 5;
  vector<double> init_mix(nStates);
  for (int i = 0; i < nStates; i++) init_mix[i] = init[3 + i];

  vector<vector<double> >
    temp_PD = integrate_ode_rk45(FK_mixed_functor(), init_mix, 0, dt, parms,
                                 rate, x_i);

  Eigen::Matrix<double, Eigen::Dynamic, 1>
    tempPK = fTwoCpt(dt[0],
                     to_vector(parms),
                     to_vector(init),
                     rate);

  vector<vector<double> > temp_mix(1);
  temp_mix[0].resize(8);
  for (int i = 0 ; i < 3; i++) temp_mix[0][i] = tempPK(i);
  for (int i = 0 ; i < 5; i++) temp_mix[0][3 + i] = temp_PD[0][i];

  clock_t end2 = clock();
  std::cout << "(dbl) CPU time for mix solver: "
            <<  (float)(end2 - start2) / CLOCKS_PER_SEC
            << " s"  << std::endl;

  // Compare output of solver to what is produced by mrgsolve
  vector<double> mrgSolution(8);
  mrgSolution[0] = 10826.82;
  mrgSolution[1] = 44779.91;
  mrgSolution[2] = 14296.11;
  mrgSolution[3] = -0.04823222;
  mrgSolution[4] = -0.0006261043;
  mrgSolution[5] = -5.620476e-06;
  mrgSolution[6] = -3.883612e-08;
  mrgSolution[7] = -2.188286e-10;

  for (int i = 0; i < 8; i++) {
    // FIX ME - what value should I use for the relative error?
    double rel_err  = std::max(std::abs(1e-3 * mrgSolution[i]), 5e-13);
    EXPECT_NEAR(mrgSolution[i], temp[0][i], rel_err);
    EXPECT_NEAR(mrgSolution[i], temp_mix[0][i], rel_err);
  }
}


TEST(mixed_solver, var) {
  using std::vector;
  using stan::math::var;
  using stan::math::integrate_ode_rk45;
  using stan::math::to_vector;

  vector<var> init(8, 0);
  init[0] = 80000;  // bolus dose in the gut
  double t0 = 0;
  vector<double> t(1);
  t[0] = 1;

  vector<var> parms(9);
  parms[0] = 10;
  parms[1] = 15;
  parms[2] = 35;
  parms[3] = 105;
  parms[4] = 2;
  parms[5] = 125;
  parms[6] = 5;
  parms[7] = 3e-4;
  parms[8] = 0.17;

  vector<double> x_r(1, 0);
  vector<int> x_i(1, 0);

  // need to make sure we pass dt (for fTwoCpt)
  vector<double> dt(1);
  dt[0] = t[0] - t0;

  // Use full numerical method
  clock_t start = clock();

  vector<vector<var> >
    temp = integrate_ode_rk45(FK_functor(), init, 0, dt, parms, x_r, x_i);

  clock_t end = clock();
  std::cout << "(var) CPU time for num solver: "
            <<  (float)(end - start) / CLOCKS_PER_SEC
            << " s"  << std::endl;

  // Use mix solver
  clock_t start2 = clock();

  parms.push_back(init[0]);
  parms.push_back(init[1]);
  parms.push_back(init[2]);

  vector<double> rate(8, 0);  // required for fTwoCpt

  int nStates = 5;
  vector<var> init_mix(nStates);
  for (int i = 0; i < nStates; i++) init_mix[i] = init[3 + i];

  vector<vector<var> >
    temp_PD = integrate_ode_rk45(FK_mixed_functor(), init_mix, 0, dt, parms,
                                 rate, x_i);

  Eigen::Matrix<var, Eigen::Dynamic, 1>
    tempPK = fTwoCpt(dt[0],
                     to_vector(parms),
                     to_vector(init),
                     rate);

  vector<vector<var> > temp_mix(1);
  temp_mix[0].resize(8);
  for (int i = 0 ; i < 3; i++) temp_mix[0][i] = tempPK(i);
  for (int i = 0 ; i < 5; i++) temp_mix[0][3 + i] = temp_PD[0][i];

  clock_t end2 = clock();
  std::cout << "(var) CPU time for mix solver: "
            <<  (float)(end2 - start2) / CLOCKS_PER_SEC
            << " s"  << std::endl;

  // Compare output of solver to what is produced by mrgsolve
  vector<double> mrgSolution(8);
  mrgSolution[0] = 10826.82;
  mrgSolution[1] = 44779.91;
  mrgSolution[2] = 14296.11;
  mrgSolution[3] = -0.04823222;
  mrgSolution[4] = -0.0006261043;
  mrgSolution[5] = -5.620476e-06;
  mrgSolution[6] = -3.883612e-08;
  mrgSolution[7] = -2.188286e-10;

  for (int i = 0; i < 8; i++) {
    // double rel_err  = std::max(std::abs(1e-6 * mrgSolution[i]), 1e-14);
    double rel_err = std::max(std::abs(1e-3 * mrgSolution[i]), 1e-14);
    EXPECT_NEAR(mrgSolution[i], temp[0][i].val(), rel_err);
    EXPECT_NEAR(mrgSolution[i], temp_mix[0][i].val(), rel_err);
  }

  // Test autodiff against finite diff (note: only treating the case
  // where both the initial state and parameters are var)
  double diff = 1e-8, diff2 = 1e-3;
  test_ode_finite_diff_vv_rel(FK_functor(), 0, dt,
                              stan::math::value_of(init),
                              stan::math::value_of(parms),
                              x_r, x_i, diff, diff2);
}
