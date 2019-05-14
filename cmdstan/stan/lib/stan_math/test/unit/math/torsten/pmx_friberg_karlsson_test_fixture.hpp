#ifndef TEST_UNIT_TORSTEN_PK_FRIBERG_KARLSSON_MODEL_TEST_FIXTURE
#define TEST_UNIT_TORSTEN_PK_FRIBERG_KARLSSON_MODEL_TEST_FIXTURE

#include <stan/math/rev/mat.hpp>
#include <gtest/gtest.h>
#include <boost/numeric/odeint.hpp>
#include <stan/math/torsten/pmx_onecpt_model.hpp>
#include <stan/math/torsten/pmx_ode_model.hpp>
#include <iostream>
#include <sstream>
#include <vector>
#include <string>

struct FribergKarlsson {
   // parms contains both the PK and the PD parameters.
   // x contains both the PK and the PD states.
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

struct FribergKarlssonTest : public testing::Test {
  const FribergKarlsson f;
  int nCmt;
  int nt;
  std::vector<std::vector<double> > theta;
  std::vector<std::vector<double> > biovar;
  std::vector<std::vector<double> > tlag;
  std::vector<double> time;
  std::vector<double> amt;
  std::vector<double> rate;
  std::vector<int> cmt;
  std::vector<int> evid;
  std::vector<double> ii;
  std::vector<int> addl;
  std::vector<int> ss;

  void SetUp() {
    // make sure memory's clean before starting each test
    stan::math::recover_memory();
  }

  FribergKarlssonTest() :
    f(),
    nCmt(8),
    nt(10),
    // CL , Q , Vc , Vp , ka , MTT , Circ0 , alpha , gamma
    theta(1, {10, 15, 35, 105, 2.0, 125, 5, 3e-4, 0.17}),
    biovar(1, {1, 1, 1, 1, 1, 1, 1, 1}),
    tlag(1, {0, 0, 0, 0, 0, 0, 0, 0}),
    time(nt),
    amt(nt, 0),
    rate(nt, 0),
    cmt(nt, 2),
    evid(nt, 0),
    ii(nt, 0),
    addl(nt, 0),
    ss(nt, 0)
  {
    time[0] = 0.0;
    for(int i = 1; i < 10; i++) time[i] = time[i - 1] + 1.25;

    amt[0] = 80 * 1000.0;
    cmt[0] = 1;
    evid[0] = 1;
    ii[0] = 12;

    SetUp();
  }

  void resize(int n) {
    nt = n;
    time.resize(nt);
    amt .resize(nt);
    rate.resize(nt);
    cmt .resize(nt);
    evid.resize(nt);
    ii  .resize(nt);
    addl.resize(nt);
    ss  .resize(nt);
  }
};

#endif
