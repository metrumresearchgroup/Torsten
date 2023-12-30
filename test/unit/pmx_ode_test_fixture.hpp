#ifndef TEST_UNIT_TORSTEN_PMX_ODE_TEST_FIXTURE
#define TEST_UNIT_TORSTEN_PMX_ODE_TEST_FIXTURE

// #include <stan/math/rev/mat.hpp>
#include <gtest/gtest.h>
#include <stan/math/rev/fun/fmax.hpp>
#include <boost/numeric/odeint.hpp>
#include <stan/math/torsten/test/unit/test_macros.hpp>
#include <test/unit/math/prim/functor/harmonic_oscillator.hpp>
#include <test/unit/math/prim/functor/lorenz.hpp>
#include <nvector/nvector_serial.h>
#include <sundials/sundials_context.h>
#include <iostream>
#include <sstream>
#include <vector>
#include <string>

/* 
 * The problem is from chemical kinetics, from CVODES examples
 */ 
struct chemical_kinetics {
  template <typename T0, typename T1, typename T2>
  inline std::vector<typename stan::return_type<T1, T2>::type>
  operator()(const T0& t_in,
             const std::vector<T1>& y,
             const std::vector<T2>& theta, const std::vector<double>& x_r,
             const std::vector<int>& x_i, std::ostream* msgs) const {
    if (y.size() != 3)
      throw std::domain_error(
          "the ODE RHS function was called with inconsistent state");

    std::vector<typename stan::return_type<T1, T2>::type> rhs(3);
    rhs[0] = -theta[0] * y[0] + theta[1] * y[1] * y[2];
    rhs[1] =  theta[0] * y[0] - theta[1] * y[1] * y[2] - theta[2] * y[1] * y[1];
    rhs[2] =  theta[2] * y[1] * y[1];

    return rhs;
  }
};

struct chemical_kinetics_eigen {
  template <typename T0, typename T1, typename T2>
  inline Eigen::Matrix<typename stan::return_type<T1, T2>::type,-1,1>
  operator()(const T0& t_in,
             const Eigen::Matrix<T1, -1, 1>& y, std::ostream* msgs,
             const std::vector<T2>& theta, const std::vector<double>& x_r,
             const std::vector<int>& x_i) const {
    if (y.size() != 3)
      throw std::domain_error(
          "the ODE RHS function was called with inconsistent state");

    Eigen::Matrix<typename stan::return_type<T1, T2>::type, -1, 1> rhs(3);
    rhs[0] = -theta[0] * y[0] + theta[1] * y[1] * y[2];
    rhs[1] =  theta[0] * y[0] - theta[1] * y[1] * y[2] - theta[2] * y[1] * y[1];
    rhs[2] =  theta[2] * y[1] * y[1];

    return rhs;
  }
};

struct lorenz_ode_eigen_fun {
  template <typename T0, typename T1, typename T2>
  inline Eigen::Matrix<stan::return_type_t<T1, T2>, -1,1>
  operator()(const T0& t_in, const Eigen::Matrix<T1, -1, 1>& y,
             std::ostream* msgs,
             const std::vector<T2>& theta, const std::vector<double>& x,
             const std::vector<int>& x_int) const {
    Eigen::Matrix<stan::return_type_t<T1, T2>, -1, 1> res(3);
    res[0] = (theta.at(0) * (y[1] - y[0]));
    res[1] = (theta.at(1) * y[0] - y[1] - y[0] * y[2]);
    res[2] = (-theta.at(2) * y[2] + y[0] * y[1]);
    return res;
  }
};

struct TorstenOdeTest : public testing::Test {
  // for events generation
  const int nt;
  std::vector<double> time;
  std::vector<double> amt;
  std::vector<double> rate;
  std::vector<int> cmt;
  std::vector<int> evid;
  std::vector<double> ii;
  std::vector<int> addl;
  std::vector<int> ss;
  std::vector<std::vector<double> > pMatrix;
  std::vector<std::vector<double> > biovar;
  std::vector<std::vector<double> > tlag;

  // for ODE integrator
  double t0;
  std::vector<double> x_r;
  std::vector<int> x_i;
  double rtol;
  double atol;
  int max_num_steps;
  std::ostream* msgs;

  void SetUp() {
    // make sure memory's clean before starting each test
    stan::math::recover_memory();
  }
  TorstenOdeTest() :
    nt(54),
    time{0, 0.083, 0.167, 0.25, 0.5, 0.75, 1, 1.5, 2, 3, 4, 6, 8, 12,
      12.083, 12.167, 12.25, 12.5, 12.75, 13, 13.5, 14, 15, 16, 18, 20, 24, 36,
      48, 60, 72, 84, 96, 108, 120, 132, 144, 156, 168, 168.083, 168.167, 168.25,
      168.5, 168.75, 169, 169.5, 170, 171, 172, 174, 176, 180, 186, 192},
    amt(nt, 0),
    rate(nt, 0),
    cmt(nt, 2),
    evid(nt, 0),
    ii(nt, 0),
    addl(nt, 0),
    ss(nt, 0),
    pMatrix{ {7.4, 1.1, 28, 78, 68 } },
    biovar{ { 1, 1, 1 } },
    tlag{ { 0, 0, 0 } },
    t0(0.0),
    rtol             {1.E-10},
    atol             {1.E-10},
    max_num_steps    {100000},
    msgs             {nullptr} {
      amt[0]  = 80000;
      cmt[0]  = 1;
      evid[0] = 1;
      ii[0]   = 12;
      addl[0] = 14;
      SetUp();
  }

  template <typename Ode>
  void test_cvodes_system(Ode& ode,
                          std::vector<double>& y0,
                          const std::vector<double>& theta,
                          const std::vector<double>& ts) {
    EXPECT_EQ(ode.N, y0.size());
    EXPECT_EQ(ode.M, theta.size());

    sundials::Context sundials_context_;
    N_Vector ydot = N_VNew_Serial(ode.N, sundials_context_);
    N_Vector y = N_VMake_Serial(ode.N, y0.data(), sundials_context_);
    
    EXPECT_EQ(ode.cvodes_rhs(ts.back(), y, ydot, static_cast<void*>(&ode)), 0);
    auto dydt = ode.f_(ts.back(), y0, theta, x_r, x_i, msgs);
    for (size_t i = 0; i < ode.N; ++i)
      EXPECT_EQ(NV_Ith_S(ydot, i), dydt[i]);
  }
};

struct TorstenOdeTest_sho : public TorstenOdeTest {
  const harm_osc_ode_fun f;
  const harm_osc_ode_fun_eigen f_eigen;
  std::vector<double> ts;
  std::vector<double> theta;
  std::vector<double> y0;
  Eigen::VectorXd y0_vec;

  using F = harm_osc_ode_fun;
  using F_eigen = harm_osc_ode_fun_eigen;

  TorstenOdeTest_sho() :
    f(),
    f_eigen(),
    ts           {0.1, 0.2, 0.3, 10},
    theta        {0.15},
    y0           {1.0, 0.0},
    y0_vec(2)
  {
    y0_vec[0] = y0[0];
    y0_vec[1] = y0[1];
  }
};

struct TorstenOdeTest_chem : public TorstenOdeTest {
  const chemical_kinetics f;
  const chemical_kinetics_eigen f_eigen;
  std::vector<double> ts;
  std::vector<double> theta;
  std::vector<double> y0;
  Eigen::VectorXd y0_vec;

  using F = chemical_kinetics;
  using F_eigen = chemical_kinetics_eigen;

  TorstenOdeTest_chem() :
    f(),
    f_eigen(),
    ts    {0.4, 4.0, 40.0},
    theta {0.04, 1.E4, 3.E7},
    y0    {1.0, 0.0, 0.0},
    y0_vec(3)
  {
    y0_vec << 1.0, 0.0, 0.0;
  }
};

struct TorstenOdeTest_lorenz : public TorstenOdeTest {
  const lorenz_ode_fun f;
  const lorenz_ode_eigen_fun f_eigen;
  std::vector<double> ts;
  std::vector<double> theta;
  std::vector<double> y0;
  Eigen::VectorXd y0_vec;

  using F = lorenz_ode_fun;
  using F_eigen = lorenz_ode_eigen_fun;

  TorstenOdeTest_lorenz() :
    f(),
    f_eigen(),
    ts        {0.1, 1.0, 10.0},
    theta     {10.0, 28.0, 8.0/3.0},
    y0        {10.0, 1.0, 1.0},
    y0_vec    (3)
  {
    y0_vec << 10.0, 1.0, 1.0;
  }
};

struct TwoCptNeutModelODE {
  template <typename T0, typename T1, typename T2>
  inline std::vector<typename stan::return_type<T1, T2>::type>
  operator()(const T0& t_in,
             const std::vector<T1>& x,
             const std::vector<T2>& theta,
             const std::vector<double>& x_r,
             const std::vector<int>& x_i,
             std::ostream* msgs) const {
    using scalar_type = typename stan::return_type<T1, T2>::type;

    const T2& CL    = theta[0];
    const T2& Q     = theta[1];
    const T2& V1    = theta[2];
    const T2& V2    = theta[3];
    const T2& ka    = theta[4];
    const T2& mtt   = theta[5];
    const T2& circ0 = theta[6];
    const T2& gamma = theta[7];
    const T2& alpha = theta[8];

    T2 k10 = CL / V1;
    T2 k12 = Q / V1;
    T2 k21 = Q / V2;
    T2 ktr = 4 / mtt;
    
    std::vector<scalar_type> dxdt(8);
    dxdt[0] = -ka * x[0];
    dxdt[1] = ka * x[0] - (k10 + k12) * x[1] + k21 * x[2];
    dxdt[2] = k12 * x[1] - k21 * x[2];    
    scalar_type conc = x[1]/V1;
    scalar_type EDrug = alpha * conc;

    // x[4], x[5], x[6], x[7] and x[8] are differences from circ0.
    scalar_type prol = x[3] + circ0;
    scalar_type transit1 = x[4] + circ0;
    scalar_type transit2 = x[5] + circ0;
    scalar_type transit3 = x[6] + circ0;
    scalar_type circ = stan::math::fmax(stan::math::machine_precision(), x[7] + circ0);

    dxdt[3] = ktr * prol * ((1 - EDrug) * (pow(circ0 / circ, gamma) - 1));
    dxdt[4] = ktr * (prol - transit1);
    dxdt[5] = ktr * (transit1 - transit2);
    dxdt[6] = ktr * (transit2 - transit3);
    dxdt[7] = ktr * (transit3 - circ);

    return dxdt;
  }
};

struct TwoCptNeutModelODE_eigen {
  template <typename T0, typename T1, typename T2>
  inline Eigen::Matrix<typename stan::return_type<T1, T2>::type, -1, 1>
  operator()(const T0& t_in,
             const Eigen::Matrix<T1, -1, 1>& x,
             std::ostream* msgs,
             const std::vector<T2>& theta,
             const std::vector<double>& x_r,
             const std::vector<int>& x_i) const {
    using scalar_type = typename stan::return_type<T1, T2>::type;

    const T2& CL    = theta[0];
    const T2& Q     = theta[1];
    const T2& V1    = theta[2];
    const T2& V2    = theta[3];
    const T2& ka    = theta[4];
    const T2& mtt   = theta[5];
    const T2& circ0 = theta[6];
    const T2& gamma = theta[7];
    const T2& alpha = theta[8];

    T2 k10 = CL / V1;
    T2 k12 = Q / V1;
    T2 k21 = Q / V2;
    T2 ktr = 4 / mtt;
    
    Eigen::Matrix<scalar_type, -1, 1> dxdt(8);
    dxdt[0] = -ka * x[0];
    dxdt[1] = ka * x[0] - (k10 + k12) * x[1] + k21 * x[2];
    dxdt[2] = k12 * x[1] - k21 * x[2];    
    scalar_type conc = x[1]/V1;
    scalar_type EDrug = alpha * conc;

    // x[4], x[5], x[6], x[7] and x[8] are differences from circ0.
    scalar_type prol = x[3] + circ0;
    scalar_type transit1 = x[4] + circ0;
    scalar_type transit2 = x[5] + circ0;
    scalar_type transit3 = x[6] + circ0;
    scalar_type circ = stan::math::fmax(stan::math::machine_precision(), x[7] + circ0);

    dxdt[3] = ktr * prol * ((1 - EDrug) * (pow(circ0 / circ, gamma) - 1));
    dxdt[4] = ktr * (prol - transit1);
    dxdt[5] = ktr * (transit1 - transit2);
    dxdt[6] = ktr * (transit2 - transit3);
    dxdt[7] = ktr * (transit3 - circ);

    return dxdt;
  }
};

struct TorstenOdeTest_neutropenia : public TorstenOdeTest {
  const TwoCptNeutModelODE f;
  const TwoCptNeutModelODE_eigen f_eigen;
  std::vector<double> ts;
  std::vector<double> theta;
  std::vector<double> y0;
  Eigen::VectorXd y0_vec;

  using F = TwoCptNeutModelODE;
  using F_eigen = TwoCptNeutModelODE_eigen;

  TorstenOdeTest_neutropenia() :
    f(),
    f_eigen(),
    ts        {0.083, 0.167, 0.25, 0.5, 0.75, 1, 1.5, 2,3,4,6,8,12,24,36,48,60,72,200},
    theta     {10, 15, 35, 105, 2, 125, 5, 0.17, 2.0e-4},
    y0        {100.0, 10.0, 10.0, 10.0, 10.0, 10.0, 10.0, 10.0},
    y0_vec(8)
  {
    y0_vec << 100.0, 10.0, 10.0, 10.0, 10.0, 10.0, 10.0, 10.0;
  }
};

#endif
