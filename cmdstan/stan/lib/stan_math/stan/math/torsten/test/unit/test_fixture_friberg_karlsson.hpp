#ifndef TEST_UNIT_TORSTEN_FRIEBERG_KARLSSON_TEST_FIXTURE
#define TEST_UNIT_TORSTEN_FRIEBERG_KARLSSON_TEST_FIXTURE

#include <stan/math/rev/core/recover_memory.hpp>
#include <gtest/gtest.h>
#include <boost/numeric/odeint.hpp>
#include <stan/math/torsten/test/unit/test_fixture_model.hpp>
#include <stan/math/torsten/pmx_ode_model.hpp>
#include <iostream>
#include <sstream>
#include <vector>
#include <string>

struct FribergKarlssonFunc {
   // parms contains both the PK and the PD parameters.
   // x contains both the PK and the PD states.
  template <typename T0, typename T1, typename T2>
  Eigen::Matrix<typename stan::return_type_t<T0, T1, T2>, -1, 1>
  operator()(const T0& t,
             const Eigen::Matrix<T1, -1, 1>& x,
             std::ostream* pstream__,
             const std::vector<T2>& parms,
             const std::vector<double>& x_r,
             const std::vector<int>& x_i) const {
    using scalar = typename stan::return_type_t<T0, T1, T2>;

    // PK variables
    T2
      CL = parms[0],
      Q = parms[1],
      VC = parms[2],
      VP = parms[3],
      ka = parms[4],
      k10 = CL / VC,
      k12 = Q / VC,
      k21 = Q / VP;

    // PD variables
    T2
      MTT = parms[5],
      circ0 = parms[6],
      gamma = parms[7],
      alpha = parms[8],
      ktr = 4 / MTT;
    typename stan::return_type_t<T1, T2>
      prol = x[3] + circ0,
      transit1 = x[4] + circ0,
      transit2 = x[5] + circ0,
      transit3 = x[6] + circ0,
      circ = stan::math::fmax(stan::math::machine_precision(), x[7] + circ0);

    Eigen::Matrix<scalar, -1, 1> dxdt(8);
    dxdt[0] = -ka * x[0];
    dxdt[1] = ka * x[0] - (k10 + k12) * x[1] + k21 * x[2];
    dxdt[2] = k12 * x[1] - k21 * x[2];

    scalar conc = x[1] / VC;
    scalar Edrug = alpha * conc;

    dxdt[3] = ktr * prol * (((1 - Edrug) * pow((circ0 / circ), gamma)) - 1);
    dxdt[4] = ktr * (prol - transit1);
    dxdt[5] = ktr * (transit1 - transit2);
    dxdt[6] = ktr * (transit2 - transit3);
    dxdt[7] = ktr * (transit3 - circ);

    return dxdt;
  }
};

template<typename T>
struct test_fk : public TorstenPMXTest<test_fk<T> > {
  test_fk() {
    this -> ncmt = 8;

    this -> reset_events(10);

    this -> theta.resize(1);

    // CL , Q , Vc , Vp , ka , MTT , Circ0 , alpha , gamma
    this -> theta[0] = {10, 15, 35, 105, 2.0, 125, 5, 0.17, 3e-4};

    this -> biovar.resize(1);
    this -> biovar[0] = {1, 1, 1, 1, 1, 1, 1, 1};
    this -> tlag.resize(1);
    this -> tlag[0] = {0, 0, 0, 0, 0, 0, 0, 0};

    for(int i = 0; i < this -> nt; i++) {
      this -> time[i] = i * 1.25; 
    }

    this -> amt[0] = 80000.0;
    this -> cmt[0] = 1;
    this -> evid[0] = 1;
    this -> ii[0] = 12;
    this -> addl[0] = 14;
  }
};

#endif
