#ifndef TEST_UNIT_TORSTEN_TWOCPT_TEST_FIXTURE
#define TEST_UNIT_TORSTEN_TWOCPT_TEST_FIXTURE

#include <stan/math/rev/core/recover_memory.hpp>
#include <gtest/gtest.h>
#include <boost/numeric/odeint.hpp>
#include <stan/math/torsten/test/unit/test_fixture_model.hpp>
#include <stan/math/torsten/pmx_twocpt_model.hpp>
#include <stan/math/torsten/pmx_ode_model.hpp>
#include <iostream>
#include <sstream>
#include <vector>
#include <string>

template<typename T>
struct test_twocpt : public TorstenPMXTest<test_twocpt<T> > {
  test_twocpt() {
    this -> ncmt = 3;

    this -> reset_events(10);

    this -> theta.resize(1);
    this -> theta[0] = {5, 8, 20, 70, 1.2}; // CL, Q, V2, V3, ka
    this -> biovar.resize(1);
    this -> biovar[0] = {1, 1, 1};
    this -> tlag.resize(1);
    this -> tlag[0] = {0, 0, 0};
    this -> pMatrix.resize(1);
    this -> pMatrix[0].resize(3, 3);

    auto CL = this -> theta[0][0];
    auto Q = this -> theta[0][1];
    auto V2 = this -> theta[0][2];
    auto V3 = this -> theta[0][3];
    auto ka = this -> theta[0][4];
    auto k10 = CL / V2;
    auto k12 = Q / V2;
    auto k21 = Q / V3;
    this -> pMatrix[0] << -ka, 0.0, 0.0, ka, -(k10 + k12), k21, 0.0, k12, -k21;

    for(int i = 0; i < this -> nt; i++) {
      this -> time[i] = i * 0.25; 
    }
    this -> time[this -> nt - 1] = 4.0;

    this -> amt[0] = 1000.0;
    this -> cmt[0] = 1;
    this -> evid[0] = 1;
    this -> ii[0] = 12;
    this -> addl[0] = 14;
  }
};

#endif
