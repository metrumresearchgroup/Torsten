#ifndef TEST_UNIT_TORSTEN_PK_ODE_TEST_FIXTURE
#define TEST_UNIT_TORSTEN_PK_ODE_TEST_FIXTURE

#include <stan/math/rev/mat.hpp>
#include <gtest/gtest.h>
#include <boost/numeric/odeint.hpp>
#include <stan/math/torsten/dsolve/pk_cvodes_fwd_system.hpp>
#include <stan/math/torsten/dsolve/pk_cvodes_integrator.hpp>
#include <test/unit/math/prim/arr/functor/harmonic_oscillator.hpp>
#include <test/unit/math/prim/arr/functor/lorenz.hpp>
#include <iostream>
#include <sstream>
#include <vector>
#include <string>

using torsten::dsolve::PKCvodesSystem;
using torsten::dsolve::PKCvodesFwdSystem;

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
  const double t0;
  const std::vector<double> x_r;
  const std::vector<int> x_i;
  double rtol;
  double atol;
  int max_num_steps;
  std::ostream* msgs;

  void SetUp() {
    // make sure memory's clean before starting each test
    stan::math::recover_memory();
  }
  TorstenOdeTest() :
    t0(0.0),
    rtol             {1.E-10},
    atol             {1.E-10},
    max_num_steps    {100000},
    msgs             {nullptr} {
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
    EXPECT_EQ(ode.rhs()(ts.back(), y, res, user_data), 0);
    auto fval = ode.f()(ts.back(), y0, theta, x_r, x_i, msgs);
    for (size_t i = 0; i < y0.size(); ++i)
      EXPECT_EQ(NV_Ith_S(res, i), fval[i]);
  }

  void test_cvodes_fwd_sens(std::vector<stan::math::var>& theta,
                            std::vector<
                            std::vector<stan::math::var>>& pk_y,
                            std::vector<
                            std::vector<stan::math::var>>& stan_y,
                            double fval_eps,
                            double sens_eps) {
    EXPECT_EQ(pk_y.size(), stan_y.size());
    EXPECT_EQ(pk_y[0].size(), stan_y[0].size());

    for (size_t i = 0; i < pk_y.size(); ++i) {
      for (size_t j = 0; j < pk_y[0].size(); ++j) {
        EXPECT_NEAR(pk_y[i][j].val(), stan_y[i][j].val(), fval_eps);
      }
    }

    std::vector<double> g, g1;
    for (size_t i = 0; i < pk_y.size(); ++i) {
      for (size_t j = 0; j < pk_y[0].size(); ++j) {
        stan::math::set_zero_all_adjoints();
        pk_y[i][j].grad(theta, g);
        stan::math::set_zero_all_adjoints();
        stan_y[i][j].grad(theta, g1);
        for (size_t m = 0; m < g.size(); ++m) {
          EXPECT_NEAR(g[m], g1[m], sens_eps);
        }
      }
    }
  }
};

struct TorstenOdeTest_sho : public TorstenOdeTest {
  const harm_osc_ode_fun f;
  const std::vector<double> ts;
  const std::vector<double> theta;
  const std::vector<double> y0;

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
  const std::vector<double> ts;
  const std::vector<double> theta;
  const std::vector<double> y0;

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
  const std::vector<double> ts;
  const std::vector<double> theta;
  const std::vector<double> y0;

  using F = lorenz_ode_fun;

  TorstenOdeTest_lorenz() :
    f(),
    ts        {0.1, 1.0, 10.0},
    theta     {10.0, 28.0, 8.0/3.0},
    y0        {10.0, 1.0, 1.0}
  {}
};



#endif
