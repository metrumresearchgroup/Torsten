#ifndef TEST_UNIT_TORSTEN_PK_ONECPT_EFFCPT_MODEL_TEST_FIXTURE
#define TEST_UNIT_TORSTEN_PK_ONECPT_EFFCPT_MODEL_TEST_FIXTURE

#include <stan/math/rev/core/recover_memory.hpp>
#include <gtest/gtest.h>
#include <boost/numeric/odeint.hpp>
#include <stan/math/torsten/test/unit/test_fixture_model.hpp>
#include <stan/math/torsten/pmx_solve_onecpt_effcpt.hpp>
#include <stan/math/torsten/pmx_ode_model.hpp>
#include <iostream>
#include <sstream>
#include <vector>
#include <string>

template<typename T>
struct test_onecpt_effcpt : public TorstenPMXTest<test_onecpt_effcpt<T>> {
  test_onecpt_effcpt() {
    this -> ncmt = 3;

    this -> reset_events(10);

    this -> theta.resize(1);
    this -> theta[0] = {10, 80, 1.2, 0.5};
    this -> biovar.resize(1);
    this -> biovar[0] = {1, 1, 1};
    this -> tlag.resize(1);
    this -> tlag[0] = {0, 0, 0};
    this -> pMatrix.resize(1);
    this -> pMatrix[0].resize(3, 3);
    this -> pMatrix[0] << - this -> theta[0][2], 0.0, 0.0,
      this -> theta[0][2], - this -> theta[0][0]/ this -> theta[0][1], 0.0,
      0.0, this -> theta[0][3], - this -> theta[0][3];

    this -> time[0] = 0.0;
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
