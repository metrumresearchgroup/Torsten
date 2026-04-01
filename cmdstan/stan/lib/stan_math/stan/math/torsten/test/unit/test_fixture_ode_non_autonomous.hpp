#ifndef TEST_UNIT_TORSTEN_ONECPT_NON_AUTONOMOUS_TEST_FIXTURE
#define TEST_UNIT_TORSTEN_ONECPT_NON_AUTONOMOUS_TEST_FIXTURE

#include <stan/math/rev/core/recover_memory.hpp>
#include <gtest/gtest.h>
#include <boost/numeric/odeint.hpp>
#include <stan/math/torsten/test/unit/test_fixture_model.hpp>
#include <stan/math/torsten/pmx_ode_model.hpp>
#include <iostream>
#include <sstream>
#include <vector>
#include <string>

template <typename T0, typename T1, typename T2, typename T3>
Eigen::Matrix<typename stan::return_type_t<T0, T1, T2, T3>, -1, 1>
onecpt_abstime(const T0& t,
                       const Eigen::Matrix<T1, -1, 1>& x,
                       std::ostream* pstream__,
                       const std::vector<T2>& parms,
                       const std::vector<T3>& x_r,
                       const std::vector<int>& x_i) {
  typedef typename stan::return_type_t<T0, T1, T2, T3> scalar;

  scalar CL0 = parms[0], V1 = parms[1], ka = parms[2], CLSS = parms[3],
    K = parms[4];

  scalar CL = CL0 + (CLSS - CL0) * (1 - stan::math::exp(-K * t));
  scalar k10 = CL / V1;

  Eigen::Matrix<scalar, -1, 1> y(2);

  y[0] = -ka * x[0];
  y[1] = ka * x[0] - k10 * x[1];

  return y;
}

struct onecpt_abstime_functor {
  template <typename T0, typename T1, typename T2, typename T3>
  Eigen::Matrix<typename stan::return_type_t<T0, T1, T2, T3>, -1, 1>
  operator()(const T0& t,
             const Eigen::Matrix<T1, -1, 1>& x,
             std::ostream* pstream__,
             const std::vector<T2>& parms,
             const std::vector<T3>& x_r,
             const std::vector<int>& x_i) const {
    return onecpt_abstime(t, x, pstream__, parms, x_r, x_i);
  }
};

template<typename T>
struct test_onecpt_non_autonomous : public TorstenPMXTest<test_onecpt_non_autonomous<T> > {
  test_onecpt_non_autonomous() {
    this -> ncmt = 2;

    this -> reset_events(10);

    this -> theta.resize(1);
    this -> theta[0] = {10, 80, 1.2, 2, 1}; // CL0, VC, ka, CLSS, K
    this -> biovar.resize(1);
    this -> biovar[0] = {1, 1};
    this -> tlag.resize(1);
    this -> tlag[0] = {0, 0};

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
