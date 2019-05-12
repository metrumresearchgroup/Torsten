#ifndef TEST_UNIT_TORSTEN_PMX_ODE_TEST_FIXTURE
#define TEST_UNIT_TORSTEN_PMX_ODE_TEST_FIXTURE

#include <stan/math/rev/mat.hpp>
#include <gtest/gtest.h>
#include <boost/numeric/odeint.hpp>
#include <test/unit/math/torsten/test_util.hpp>
#include <test/unit/math/prim/arr/functor/harmonic_oscillator.hpp>
#include <test/unit/math/prim/arr/functor/lorenz.hpp>
#include <iostream>
#include <sstream>
#include <vector>
#include <string>

using torsten::dsolve::PMXCvodesSystem;
using torsten::dsolve::PMXCvodesFwdSystem;

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
                          const std::vector<double>& y0,
                          const std::vector<double>& theta,
                          const std::vector<double>& ts) {
    auto user_data = ode.to_user_data();
    Ode* odep = static_cast<Ode* >(user_data);

    EXPECT_EQ(odep -> n(), y0.size());
    EXPECT_EQ(odep -> n_par(), theta.size());
    auto y = odep -> nv_y();
    for (size_t i = 0; i < ode.n(); ++i) {
      EXPECT_EQ(NV_Ith_S(y, i), y0[i]);
    }

    N_Vector res = N_VNew_Serial(y0.size());
    EXPECT_EQ(torsten::dsolve::cvodes_rhs<Ode>()(ts.back(), y, res, user_data), 0);
    auto fval = ode.f()(ts.back(), y0, theta, x_r, x_i, msgs);
    for (size_t i = 0; i < y0.size(); ++i)
      EXPECT_EQ(NV_Ith_S(res, i), fval[i]);
  }

  void test_cvodes_fwd_sens(std::vector<stan::math::var>& theta,
                            std::vector<std::vector<stan::math::var>>& pk_y,
                            std::vector<std::vector<stan::math::var>>& stan_y,
                            double fval_eps,
                            double sens_eps) {
    torsten::test::test_grad(theta, pk_y, stan_y, fval_eps, sens_eps);
  }

  void test_cvodes_fwd_sens(std::vector<stan::math::var>& theta,
                            Eigen::MatrixXd& pk_y,
                            std::vector<std::vector<stan::math::var>>& stan_y,
                            double fval_eps,
                            double sens_eps) {
    EXPECT_EQ(pk_y.rows(), stan_y.size());
    EXPECT_EQ(pk_y.cols() % stan_y[0].size(), 0);
    const int n = pk_y.cols() / stan_y[0].size();
    EXPECT_EQ(n, 1 + theta.size());

    std::vector<double> g;
    for (int i = 0; i < pk_y.rows(); ++i) {
      for (size_t j = 0; j < stan_y[0].size(); ++j) {
        EXPECT_NEAR(pk_y(i, n * j), stan_y[i][j].val(), fval_eps);
        stan::math::set_zero_all_adjoints();
        stan_y[i][j].grad(theta, g);
        for (size_t k = 0; k < g.size(); ++k) {
          EXPECT_NEAR(pk_y(i, n * j + k + 1), g[k], sens_eps);
        }
      }
    }
  }
};

struct TorstenOdeTest_sho : public TorstenOdeTest {
  const harm_osc_ode_fun f;
  std::vector<double> ts;
  std::vector<double> theta;
  std::vector<double> y0;

  using F = harm_osc_ode_fun;

  TorstenOdeTest_sho() :
    f(),
    ts           {0.1, 0.2, 0.3, 10},
    theta        {0.15},
    y0           {1.0, 0.0}
  {}
};

struct TorstenOdeTest_chem : public TorstenOdeTest {
  const chemical_kinetics f;
  std::vector<double> ts;
  std::vector<double> theta;
  std::vector<double> y0;

  using F = chemical_kinetics;

  TorstenOdeTest_chem() :
    f(),
    ts          {0.4, 4.0, 40.0},
    theta       {0.04, 1.E4, 3.E7},
    y0          {1.0, 0.0, 0.0}
  {}
};

struct TorstenOdeTest_lorenz : public TorstenOdeTest {
  const lorenz_ode_fun f;
  std::vector<double> ts;
  std::vector<double> theta;
  std::vector<double> y0;

  using F = lorenz_ode_fun;

  TorstenOdeTest_lorenz() :
    f(),
    ts        {0.1, 1.0, 10.0},
    theta     {10.0, 28.0, 8.0/3.0},
    y0        {10.0, 1.0, 1.0}
  {}
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

struct TorstenOdeTest_neutropenia : public TorstenOdeTest {
  const TwoCptNeutModelODE f;
  std::vector<double> ts;
  std::vector<double> theta;
  std::vector<double> y0;

  using F = TwoCptNeutModelODE;

  TorstenOdeTest_neutropenia() :
    f(),
    ts        {0.083, 0.167, 0.25, 0.5, 0.75, 1, 1.5, 2,3,4,6,8,12,24,36,48,60,72,200},
    theta     {10, 15, 35, 105, 2, 125, 5, 0.17, 2.0e-4},
    y0        {100.0, 10.0, 10.0, 10.0, 10.0, 10.0, 10.0, 10.0}
  {}
};

#endif
