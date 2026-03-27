#ifdef TORSTEN_MPI_DYN

#include <stan/math.hpp>
#include <stan/math/rev/core.hpp>
#include <test/unit/math/rev/fun/util.hpp>
#include <stan/math/torsten/dsolve/pmx_integrate_ode_adams.hpp>
#include <stan/math/torsten/dsolve/pmx_integrate_ode_bdf.hpp>
#include <stan/math/torsten/dsolve/pmx_integrate_ode_group_adams.hpp>
#include <stan/math/torsten/dsolve/pmx_integrate_ode_group_bdf.hpp>
#include <stan/math/torsten/dsolve/pmx_integrate_ode_group_rk45.hpp>
#include <stan/math/torsten/mpi.hpp>
#include <stan/math/torsten/test/unit/pmx_ode_test_fixture.hpp>
#include <stan/math/torsten/test/unit/test_macros.hpp>
#include <stan/math/rev/functor/integrate_ode_bdf.hpp>
#include <nvector/nvector_serial.h>
#include <stan/math/torsten/mpi/environment.hpp>
#include <test/unit/util.hpp>
#include <gtest/gtest.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <memory>
#include <iostream>
#include <fstream>
#include <limits>
#include <vector>
#include <chrono>
#include <ctime>


using stan::math::integrate_ode_bdf;
using stan::math::integrate_ode_adams;
using stan::math::integrate_ode_rk45;
using torsten::pmx_integrate_ode_group_bdf;
using torsten::pmx_integrate_ode_group_adams;
using torsten::pmx_integrate_ode_group_rk45;
using torsten::pmx_integrate_ode_rk45;
using torsten::pmx_integrate_ode_bdf;
using torsten::pmx_integrate_ode_adams;
using torsten::dsolve::pmx_ode_group_mpi_functor_id;
using torsten::mpi::PMXDynamicLoad;
using stan::math::matrix_v;
using stan::math::var;
using std::vector;

TORSTEN_MPI_SESSION_INIT;

namespace torsten {
  namespace dsolve {
    template<typename... Args>
    inline auto torsten::dsolve::pmx_ode_group_mpi_functor::operator()(Args&&... args) const {
      if (id == 0) { const TwoCptNeutModelODE f; return f(std::forward<Args>(args)...); }
      if (id == 1) { const chemical_kinetics f; return f(std::forward<Args>(args)...); }
      if (id == 2) { const lorenz_ode_fun f; return f(std::forward<Args>(args)...); }

      // return default
      TwoCptNeutModelODE f;
      return f(std::forward<Args>(args)...);
    }

    template<>
    struct pmx_ode_group_mpi_functor_id<TwoCptNeutModelODE> { static constexpr int value = 0; };
    template<>
    struct pmx_ode_group_mpi_functor_id<chemical_kinetics> { static constexpr int value = 1; };
    template<>
    struct pmx_ode_group_mpi_functor_id<lorenz_ode_fun> { static constexpr int value = 2; };
  }
}

TEST_F(TorstenOdeTest_chem, mpi_dynamic_load_cvodes_ivp_system_bdf_mpi) {
  stan::math::mpi::Envionment::init();

  const stan::math::mpi::Communicator& pmx_comm =
    torsten::mpi::Session::ode_parm_comm();

  torsten::mpi::PMXDynamicLoad<TORSTEN_MPI_DYN_CLEANER> master(pmx_comm);

  if (pmx_comm.rank() > 0) {
    torsten::mpi::PMXDynamicLoad<TORSTEN_MPI_DYN_SLAVE> load(pmx_comm);
    load.slave();
  } else {
    const int np = 10;

    vector<int> len(np, ts.size());
    vector<double> ts_m;
    ts_m.reserve(np * ts.size());
    for (int i = 0; i < np; ++i) ts_m.insert(ts_m.end(), ts.begin(), ts.end());

    vector<vector<double> > y0_m (np, y0);
    vector<vector<double> > theta_m (np, theta);
    vector<vector<double> > x_r_m (np, x_r);
    vector<vector<int> > x_i_m (np, x_i);

    {
      vector<vector<double> > y = integrate_ode_adams(f, y0, t0, ts, theta , x_r, x_i); // NOLINT
      Eigen::MatrixXd y_m = pmx_integrate_ode_group_adams(f, y0_m, t0, len, ts_m, theta_m , x_r_m, x_i_m);
      EXPECT_EQ(y_m.cols(), ts_m.size());
      int icol = 0;
      for (size_t i = 0; i < np; ++i) {
        Eigen::MatrixXd y_i = y_m.block(0, icol, y0.size(), len[i]);
        torsten::test::test_val(y, y_i);
        icol += len[i];
      }      
    }

    {
      vector<vector<double> > y = integrate_ode_bdf(f, y0, t0, ts, theta , x_r, x_i); // NOLINT
      Eigen::MatrixXd y_m = pmx_integrate_ode_group_bdf(f, y0_m, t0, len, ts_m, theta_m , x_r_m, x_i_m);
      EXPECT_EQ(y_m.cols(), ts_m.size());
      int icol = 0;
      for (size_t i = 0; i < np; ++i) {
        Eigen::MatrixXd y_i = y_m.block(0, icol, y0.size(), len[i]);
        torsten::test::test_val(y, y_i);
        icol += len[i];
      }      
    }

    {
      vector<vector<double> > y = integrate_ode_rk45(f, y0, t0, ts, theta , x_r, x_i); // NOLINT
      Eigen::MatrixXd y_m = pmx_integrate_ode_group_rk45(f, y0_m, t0, len, ts_m, theta_m , x_r_m, x_i_m);
      EXPECT_EQ(y_m.cols(), ts_m.size());
      int icol = 0;
      for (size_t i = 0; i < np; ++i) {
        Eigen::MatrixXd y_i = y_m.block(0, icol, y0.size(), len[i]);
        torsten::test::test_val(y, y_i);
        icol += len[i];
      }      
    }
  }
}

TEST_F(TorstenOdeTest_neutropenia, mpi_dynamic_load_fwd_sensitivity_uniform_theta) {
  stan::math::mpi::Envionment::init();

  const stan::math::mpi::Communicator& pmx_comm =
    torsten::mpi::Session::ode_parm_comm();

  torsten::mpi::PMXDynamicLoad<TORSTEN_MPI_DYN_CLEANER> master(pmx_comm);

  if (pmx_comm.rank() > 0) {
    torsten::mpi::PMXDynamicLoad<TORSTEN_MPI_DYN_SLAVE> load(pmx_comm);
    load.slave();
  } else {
    const int np = 10;

    vector<var> theta_var = stan::math::to_var(theta);

    vector<int> len(np, ts.size());
    vector<double> ts_m;
    ts_m.reserve(np * ts.size());
    for (int i = 0; i < np; ++i) ts_m.insert(ts_m.end(), ts.begin(), ts.end());

    vector<vector<double> > y0_m (np, y0);
    vector<vector<var> > theta_var_m (np, stan::math::to_var(theta));
    vector<vector<double> > x_r_m (np, x_r);
    vector<vector<int> > x_i_m (np, x_i);

    {
      vector<vector<var> > y = integrate_ode_adams(f, y0, t0, ts, theta_var, x_r, x_i);
      Eigen::Matrix<var, -1, -1> y_m = pmx_integrate_ode_group_adams(f, y0_m, t0, len, ts_m, theta_var_m , x_r_m, x_i_m);
      EXPECT_EQ(y_m.cols(), ts_m.size());
      int icol = 0;
      for (int i = 0; i < np; ++i) {
        matrix_v y_i = y_m.block(0, icol, y0.size(), len[i]);
        torsten::test::test_grad(theta_var, theta_var_m[i], y, y_i, 1e-7, 2e-7);
        icol += len[i];
      }      
    }

    {
      vector<vector<var> > y = integrate_ode_bdf(f, y0, t0, ts, theta_var, x_r, x_i);
      Eigen::Matrix<var, -1, -1> y_m = pmx_integrate_ode_group_bdf(f, y0_m, t0, len, ts_m, theta_var_m , x_r_m, x_i_m);
      EXPECT_EQ(y_m.cols(), ts_m.size());
      int icol = 0;
      for (int i = 0; i < np; ++i) {
        matrix_v y_i = y_m.block(0, icol, y0.size(), len[i]);
        torsten::test::test_grad(theta_var, theta_var_m[i], y, y_i, 1e-7, 2e-7);
        icol += len[i];
      }      
    }

    {
      vector<vector<var> > y = pmx_integrate_ode_rk45(f, y0, t0, ts, theta_var, x_r, x_i);
      Eigen::Matrix<var, -1, -1> y_m = pmx_integrate_ode_group_rk45(f, y0_m, t0, len, ts_m, theta_var_m , x_r_m, x_i_m);
      EXPECT_EQ(y_m.cols(), ts_m.size());
      int icol = 0;
      for (int i = 0; i < np; ++i) {
        matrix_v y_i = y_m.block(0, icol, y0.size(), len[i]);
        torsten::test::test_grad(theta_var, theta_var_m[i], y, y_i, 1e-7, 2e-7);
        icol += len[i];
      }      
    }
  }
}

TEST_F(TorstenOdeTest_neutropenia, mpi_dynamic_load_fwd_sensitivity_uniform_theta_exception) {
  stan::math::mpi::Envionment::init();

  const stan::math::mpi::Communicator& pmx_comm =
    torsten::mpi::Session::ode_parm_comm();

  torsten::mpi::PMXDynamicLoad<TORSTEN_MPI_DYN_CLEANER> master(pmx_comm);

  if (pmx_comm.rank() > 0) {
    torsten::mpi::PMXDynamicLoad<TORSTEN_MPI_DYN_SLAVE> load(pmx_comm);
    load.slave();
  } else {
    const int np = 10;

    vector<var> theta_var = stan::math::to_var(theta);

    vector<int> len(np, ts.size());
    vector<double> ts_m;
    ts_m.reserve(np * ts.size());
    for (int i = 0; i < np; ++i) ts_m.insert(ts_m.end(), ts.begin(), ts.end());

    vector<vector<double> > y0_m (np, y0);
    vector<vector<var> > theta_var_m (np, stan::math::to_var(theta));
    vector<vector<double> > x_r_m (np, x_r);
    vector<vector<int> > x_i_m (np, x_i);

    std::stringstream expected_message;
    EXPECT_THROW(pmx_integrate_ode_group_adams(f, y0_m, t0, len, ts_m, theta_var_m , x_r_m, x_i_m,
                                               1e-10, 1e-10, 1),
                 std::exception);
    EXPECT_THROW(pmx_integrate_ode_group_bdf(f, y0_m, t0, len, ts_m, theta_var_m , x_r_m, x_i_m,
                                             1e-10, 1e-10, 1),
                 std::exception);
    EXPECT_THROW(pmx_integrate_ode_group_rk45(f, y0_m, t0, len, ts_m, theta_var_m , x_r_m, x_i_m,
                                              1e-10, 1e-10, 1),
                 std::exception);
  }
}

#endif
