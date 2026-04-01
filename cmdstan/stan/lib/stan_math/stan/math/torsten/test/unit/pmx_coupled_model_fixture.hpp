#ifndef TEST_UNIT_TORSTEN_PK_CPT_MODEL_TEST_FIXTURE
#define TEST_UNIT_TORSTEN_PK_CPT_MODEL_TEST_FIXTURE

//
// #include <gtest/gtest.h>
#include <boost/numeric/odeint.hpp>
#include <stan/math/torsten/pmx_onecpt_model.hpp>
#include <stan/math/torsten/pmx_twocpt_model.hpp>
#include <stan/math/torsten/pmx_linode_model.hpp>
#include <stan/math/torsten/pmx_ode_model.hpp>
#include <stan/math/torsten/test/unit/pmx_ode_test_fixture.hpp>
#include <iostream>
#include <sstream>
#include <vector>
#include <string>

struct CoupledOneCptODE {
  /*
   * Coupled model functor
   */
  template <typename T0, typename T1, typename T2, typename T3>
  inline
  Eigen::Matrix<typename stan::return_type_t<T0, T1, T2, T3>, -1, 1>
  operator()(const T0& t,
             const Eigen::Matrix<T1, -1, 1>& x,
             const Eigen::Matrix<T2, -1, 1>& x_pk,
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

    Eigen::Matrix<scalar, -1, 1> dxdt(3);
    dxdt << ktr * prol * ((1 - Edrug) * pow(circ0 / circ, gamma) - 1),
      ktr * (prol - transit),
      ktr * (transit - circ);

    return dxdt;
  }

  /*
   * full ODE model functor
   */
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

struct CoupledTwoCptODE {
  /*
   * Coupled model functor
   */
  template <typename T0, typename T1, typename T2, typename T3>
  Eigen::Matrix<typename stan::return_type_t<T0, T1, T2, T3>, -1, 1>
  operator()(const T0& t,
             const Eigen::Matrix<T1, -1, 1>& y,
             const Eigen::Matrix<T2, -1, 1>& y_pk,
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

    Eigen::Matrix<scalar, -1, 1> dxdt(3);
    dxdt << ktr * prol * ((1 - Edrug) * pow(circ0 / circ, gamma) - 1),
      ktr * (prol - transit),
      ktr * (transit - circ);

    return dxdt;
  }

  /*
   * full ODE model functor
   */
  template <typename T0, typename T1, typename T2>
  Eigen::Matrix<typename stan::return_type_t<T0, T1, T2>, -1, 1>
  operator()(const T0& t,
             const Eigen::Matrix<T1, -1, 1> & x,
             std::ostream* pstream_,
             const std::vector<T2>& theta,
             const std::vector<double>& x_r,
             const std::vector<int>& x_i) const {
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

    Eigen::Matrix<scalar, -1, 1> dxdt(6);
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

struct TorstenCoupledOneCptTest : public testing::Test {
  // for events generation
  const int nt;
  const int nOde;
  const int nPD;
  std::vector<double> time;
  std::vector<double> amt;
  std::vector<double> rate;
  std::vector<int> cmt;
  std::vector<int> evid;
  std::vector<double> ii;
  std::vector<int> addl;
  std::vector<int> ss;
  std::vector<std::vector<double> > parameters;
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
  TorstenCoupledOneCptTest() :
    nt(10),
    nOde(5),
    nPD(3),
    time(nt),
    amt(nt, 0),
    rate(nt, 0),
    cmt(nt, 2),
    evid(nt, 0),
    ii(nt, 0),
    addl(nt, 0),
    ss(nt, 0),
    // params: // CL // VC // ka // Mtt // Circ0 // alpha // gamma
    parameters{ {10, 35, 2.0, 125, 5, 3e-4, 0.17} },
    biovar{ { 1, 1, 1, 1, 1 } },
    tlag{ { 0, 0, 0, 0, 0 } },
    t0(0.0),
    rtol             {1.E-10},
    atol             {1.E-10},
    max_num_steps    {100000},
    msgs             {nullptr} {
      for (int i = 0; i < nt; ++i) {
        time[i] = i * 0.25;
      }
      time.back() = 4.0;
      amt[0]  = 10000;
      cmt[0]  = 1;
      evid[0] = 1;
      SetUp();
  }
};

struct TorstenCoupledTwoCptTest : public testing::Test {
  // for events generation
  const int nt;
  const int nOde;
  const int nPD;
  std::vector<double> time;
  std::vector<double> amt;
  std::vector<double> rate;
  std::vector<int> cmt;
  std::vector<int> evid;
  std::vector<double> ii;
  std::vector<int> addl;
  std::vector<int> ss;
  std::vector<std::vector<double> > parameters;
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
  TorstenCoupledTwoCptTest() :
    nt(10),
    nOde(6),
    nPD(3),
    time(nt),
    amt(nt, 0),
    rate(nt, 0),
    cmt(nt, 2),
    evid(nt, 0),
    ii(nt, 0),
    addl(nt, 0),
    ss(nt, 0),
    // CL // Q // VC // VP // ka // Mtt // Circ0 // alpha // gamma
    parameters{ {10, 15, 35, 105, 2.0, 125, 5, 3e-4, 0.17} },
    biovar{ { 1, 1, 1, 1, 1, 1 } },
    tlag{ { 0, 0, 0, 0, 0, 0 } },
    t0(0.0),
    rtol             {1.E-10},
    atol             {1.E-10},
    max_num_steps    {100000},
    msgs             {nullptr} {
      for (int i = 0; i < nt; ++i) {
        time[i] = i * 0.25;
      }
      time.back() = 4.0;
      amt[0]  = 10000;
      cmt[0]  = 1;
      evid[0] = 1;
      SetUp();
  }
};

#endif
