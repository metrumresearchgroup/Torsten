#ifndef TEST_UNIT_TORSTEN_COUPLED_ONECPT_MODEL_TEST_FIXTURE
#define TEST_UNIT_TORSTEN_COUPLED_ONECPT_MODEL_TEST_FIXTURE

#include <stan/math/rev/core/recover_memory.hpp>
#include <gtest/gtest.h>
#include <boost/numeric/odeint.hpp>
#include <stan/math/torsten/test/unit/test_fixture_model.hpp>
#include <stan/math/torsten/pmx_onecpt_model.hpp>
#include <stan/math/torsten/pmx_solve_onecpt_rk45.hpp>
#include <stan/math/torsten/pmx_solve_onecpt_bdf.hpp>
#include <stan/math/torsten/pmx_ode_model.hpp>
#include <iostream>
#include <sstream>
#include <vector>
#include <string>

struct CoupledOneCptODE {
  /**
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

  /**
   * full ODE model functor
   */
  template <typename T0, typename T1, typename T2>
  inline
  Eigen::Matrix<typename stan::return_type_t<T0, T1, T2>, -1, 1>
  operator()(const T0& t,
             const Eigen::Matrix<T1, -1, 1>& x,
             std::ostream* pstream_,
             const std::vector<T2>& theta,
             const std::vector<double>& x_r,
             const std::vector<int>& x_i) const {
    using scalar = typename stan::return_type_t<T0, T1, T2>;
    
    T2 CL = theta[0];
    T2 VC = theta[1];
    T2 ka = theta[2];
    T2 Mtt = theta[3];
    T2 circ0 = theta[4];
    T2 alpha = theta[5];
    T2 gamma = theta[6];
    T2 ktr = 4.0 / Mtt;
    scalar prol = x[2] + circ0;
    scalar transit = x[3] + circ0;
    scalar circ = x[4] + circ0;
    scalar Edrug = alpha * x[1] / VC;

    Eigen::Matrix<scalar, -1, 1>  dxdt(5);
    dxdt << - ka * x[0],
      ka * x[0] - CL / VC * x[1],
      ktr * prol * ((1 - Edrug) * pow(circ0 / circ, gamma) - 1),
      ktr * (prol - transit),
      ktr * (transit - circ);
    
    return dxdt;
  }
};

template<typename T>
struct test_coupled_onecpt : public TorstenPMXTest<test_coupled_onecpt<T> > {
  const size_t nOde;
  test_coupled_onecpt() : nOde(5) {
    this -> ncmt = nOde;           // # of PD compartments
    this -> reset_events(10);
    this -> theta.resize(1);
    this -> theta[0] = {10, 35, 2.0, 125, 5, 3e-4, 0.17};
    this -> biovar.resize(1);
    this -> biovar[0] = {1, 1, 1, 1, 1};
    this -> tlag.resize(1);
    this -> tlag[0] = {0, 0, 0, 0, 0};

    for(int i = 0; i < this -> nt; i++) {
      this -> time[i] = i * 0.25; 
    }
    this -> time[this -> nt - 1] = 4.0;

    this -> amt[0] = 10000.0;
    this -> cmt[0] = 1;
    this -> evid[0] = 1;
  }
};

#endif
