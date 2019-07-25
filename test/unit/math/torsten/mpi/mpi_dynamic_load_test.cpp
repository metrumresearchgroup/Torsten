#ifdef TORSTEN_MPI_DYN

#include <stan/math.hpp>
#include <stan/math/rev/core.hpp>
#include <stan/math/torsten/mpi/dynamic_load.hpp>
#include <test/unit/math/rev/mat/fun/util.hpp>
#include <test/unit/math/torsten/pmx_ode_test_fixture.hpp>
#include <test/unit/math/prim/arr/functor/harmonic_oscillator.hpp>
#include <nvector/nvector_serial.h>
#include <boost/mpi.hpp>
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
#include <random>

namespace torsten {
  namespace dsolve {

    template<typename... Args>
    inline auto torsten::dsolve::pmx_ode_group_mpi_functor::operator()(Args&&... args) const {
      if (id == 0) { const TwoCptNeutModelODE f; return f(std::forward<Args>(args)...); }

      // return default
      TwoCptNeutModelODE f;
      return f(std::forward<Args>(args)...);
    }

    template<>
    struct pmx_ode_group_mpi_functor_id<TwoCptNeutModelODE> { static constexpr int value = 0; };
  }
}


using torsten::dsolve::PMXCvodesFwdSystem;
using torsten::dsolve::pmx_ode_group_mpi_functor;
using torsten::dsolve::integrator_id;
using torsten::pmx_integrate_ode_group_adams;
using torsten::pmx_integrate_ode_adams;
using Eigen::MatrixXd;
using Eigen::Matrix;
using stan::math::var;
using std::vector;

TEST_F(TorstenOdeTest_neutropenia, mpi_dynamic_load_set_use_uniform_work) {
  torsten::mpi::Envionment::init();

  // size of population
  const int np = 8;
  std::vector<double> ts0 {ts};

  vector<var> theta_var = stan::math::to_var(theta);

  vector<int> len(np, ts0.size());
  vector<double> ts_m;
  ts_m.reserve(np * ts0.size());
  for (int i = 0; i < np; ++i) ts_m.insert(ts_m.end(), ts0.begin(), ts0.end());
  
  vector<vector<double> > y0_m (np, y0);
  vector<vector<var> > theta_var_m (np, theta_var);
  vector<vector<double> > theta_m (np, theta);
  vector<vector<double> > x_r_m (np, x_r);
  vector<vector<int> > x_i_m (np, x_i);

  torsten::mpi::Communicator pmx_comm(torsten::mpi::Session<NUM_TORSTEN_COMM>::env, MPI_COMM_WORLD);
  torsten::mpi::PMXDynamicLoad<TORSTEN_MPI_DYN_MASTER> load(pmx_comm);

  load.init_buf = std::vector<int>{0, 1, static_cast<int>(y0.size()), static_cast<int>(theta.size()), static_cast<int>(x_r.size()), static_cast<int>(x_i.size()), 0, 0, 0};
  load.set_work(1, y0_m, t0, len, ts_m, theta_m, x_r_m, x_i_m, 1.e-8, 1.e-10, 10000);
  std::vector<double> y0_slave;
  double t0_slave;
  std::vector<double> ts_slave;
  std::vector<double> theta_slave;
  std::vector<double> x_r_slave;
  std::vector<int> x_i_slave;
  double rtol_slave;
  double atol_slave;
  int max_num_step_slave;
  load.use_work(y0_slave, t0_slave, ts_slave, theta_slave, x_r_slave, x_i_slave, rtol_slave, atol_slave, max_num_step_slave);
  torsten::test::test_val(t0_slave, t0);
  torsten::test::test_val(y0_slave, y0);
  torsten::test::test_val(ts_slave, ts);
  torsten::test::test_val(theta_slave, theta);
  torsten::test::test_val(x_r_slave, x_r);
  torsten::test::test_val(rtol_slave, 1.e-8);
  torsten::test::test_val(atol_slave, 1.e-10);
  torsten::test::test_val(max_num_step_slave, 10000);
}

TEST_F(TorstenOdeTest_neutropenia, mpi_dynamic_load_set_use_non_uniform_work) {
  torsten::mpi::Envionment::init();

  // size of population
  const int np = 3;
  std::vector<std::vector<double>> tss(3, ts);
  tss[0].resize(ts.size() - 1);
  tss[2].resize(ts.size() - 3);

  vector<var> theta_var = stan::math::to_var(theta);

  vector<size_t> len{tss[0].size(), tss[1].size(), tss[2].size()};
  vector<double> ts_m;
  ts_m.reserve(tss[0].size() + tss[1].size() + tss[2].size());
  for (int i = 0; i < np; ++i) ts_m.insert(ts_m.end(), tss[i].begin(), tss[i].end());
  
  vector<vector<double> > y0_m (np, y0);
  vector<vector<var> > theta_var_m (np, theta_var);
  vector<vector<double> > theta_m (np, theta);
  x_r.push_back(1.0); vector<vector<double> > x_r_m (np, x_r);
  x_i.push_back(1.0); vector<vector<int> > x_i_m (np, x_i);

  torsten::mpi::Communicator pmx_comm(torsten::mpi::Session<NUM_TORSTEN_COMM>::env, MPI_COMM_WORLD);
  torsten::mpi::PMXDynamicLoad<TORSTEN_MPI_DYN_MASTER> load(pmx_comm);

  load.init_buf = std::vector<int>{0, 1, static_cast<int>(y0.size()), static_cast<int>(theta.size()), static_cast<int>(x_r.size()), static_cast<int>(x_i.size()), 0, 0, 0};
  load.set_work(1, y0_m, t0, len, ts_m, theta_m, x_r_m, x_i_m, 1.e-8, 1.e-10, 10000);
  std::vector<double> y0_slave;
  double t0_slave;
  std::vector<double> ts_slave;
  std::vector<double> theta_slave;
  std::vector<double> x_r_slave;
  std::vector<int> x_i_slave;
  double rtol_slave;
  double atol_slave;
  int max_num_step_slave;
  load.use_work(y0_slave, t0_slave, ts_slave, theta_slave, x_r_slave, x_i_slave, rtol_slave, atol_slave, max_num_step_slave);
  torsten::test::test_val(y0_slave, y0);
  torsten::test::test_val(theta_slave, theta);
  torsten::test::test_val(x_r_slave, x_r);
  torsten::test::test_val(rtol_slave, 1.e-8);
  torsten::test::test_val(atol_slave, 1.e-10);
  torsten::test::test_val(max_num_step_slave, 10000);

  for (size_t i = 0; i < tss.size(); ++i) {
    load.set_work(i, y0_m, t0, len, ts_m, theta_m, x_r_m, x_i_m, 1.e-8, 1.e-10, 10000);      
    load.use_work(y0_slave, t0_slave, ts_slave, theta_slave, x_r_slave, x_i_slave, rtol_slave, atol_slave, max_num_step_slave);    
    torsten::test::test_val(ts_slave, tss[i]);    
  }
}

TEST_F(TorstenOdeTest_neutropenia, mpi_dynamic_load_master_slave_single_work) {
  torsten::mpi::Envionment::init();

  // size of population
  const int np = 1;
  std::vector<double> ts0 {ts};

  vector<var> theta_var = stan::math::to_var(theta);

  vector<int> len(np, ts0.size());
  vector<double> ts_m;
  ts_m.reserve(np * ts0.size());
  for (int i = 0; i < np; ++i) ts_m.insert(ts_m.end(), ts0.begin(), ts0.end());
  
  vector<vector<double> > y0_m (np, y0);
  vector<vector<var> > theta_var_m (np, theta_var);
  vector<vector<double> > theta_m (np, theta);
  vector<vector<double> > x_r_m (np, x_r);
  vector<vector<int> > x_i_m (np, x_i);

  torsten::mpi::Communicator pmx_comm(torsten::mpi::Session<NUM_TORSTEN_COMM>::env, MPI_COMM_WORLD);

  torsten::dsolve::pmx_ode_group_mpi_functor fdyn(0);
  if (pmx_comm.rank == 0) {
    torsten::mpi::PMXDynamicLoad<TORSTEN_MPI_DYN_MASTER> load(pmx_comm);
    MatrixXd res = load.master(f, 1, y0_m, t0, len, ts_m, theta_m, x_r_m, x_i_m, 1.e-8, 1.e-8, 10000);
    MatrixXd sol = torsten::to_matrix(pmx_integrate_ode_adams(f, y0, t0, ts, theta, x_r, x_i, 0, 1.e-8, 1.e-8, 10000));
    torsten::test::test_val(res, sol);
    load.kill_slaves();
  } else {
    torsten::mpi::PMXDynamicLoad<TORSTEN_MPI_DYN_SLAVE> load(pmx_comm);
    load.slave();
  }
}

TEST_F(TorstenOdeTest_neutropenia, mpi_dynamic_load_master_slave_multiple_uniform_work) {
  torsten::mpi::Envionment::init();

  torsten::mpi::Communicator pmx_comm(torsten::mpi::Session<NUM_TORSTEN_COMM>::env, MPI_COMM_WORLD);
  std::vector<double> ts0 {ts};

  if (pmx_comm.rank > 0) {
    torsten::mpi::PMXDynamicLoad<TORSTEN_MPI_DYN_SLAVE> load(pmx_comm);
    load.slave();
  } else {
    torsten::mpi::PMXDynamicLoad<TORSTEN_MPI_DYN_MASTER> load(pmx_comm);
    for (int np = 1; np < 8; ++np) {
      vector<int> len(np, ts0.size());
      vector<double> ts_m;
      ts_m.reserve(np * ts0.size());
      for (int i = 0; i < np; ++i) ts_m.insert(ts_m.end(), ts0.begin(), ts0.end());
  
      vector<vector<double> > y0_m (np, y0);
      vector<vector<double> > theta_m (np, theta);
      vector<vector<double> > x_r_m (np, x_r);
      vector<vector<int> > x_i_m (np, x_i);

      MatrixXd res = load.master(f, 1, y0_m, t0, len, ts_m, theta_m, x_r_m, x_i_m, 1.e-8, 1.e-8, 10000);
      MatrixXd sol = torsten::to_matrix(pmx_integrate_ode_adams(f, y0, t0, ts, theta, x_r, x_i, 0, 1.e-8, 1.e-8, 10000));
      for (int i = 0; i < np; ++i) {
        MatrixXd res_i = res.block(0, ts.size() * i, y0.size(), ts.size());
        torsten::test::test_val(res_i, sol);
      }
    }

    load.kill_slaves();
  }
}

TEST_F(TorstenOdeTest_neutropenia, mpi_dynamic_load_bdf_master_slave_multiple_non_uniform_work) {
  torsten::mpi::Envionment::init();

  torsten::mpi::Communicator pmx_comm(torsten::mpi::Session<NUM_TORSTEN_COMM>::env, MPI_COMM_WORLD);
  std::vector<double> ts0 {ts};

  if (pmx_comm.rank > 0) {
    torsten::mpi::PMXDynamicLoad<TORSTEN_MPI_DYN_SLAVE> load(pmx_comm);
    load.slave();
  } else {
    torsten::mpi::PMXDynamicLoad<TORSTEN_MPI_DYN_MASTER> load(pmx_comm);
    const int np = 7;
    vector<int> len{19, 17, 10, 14, 18, 17, 16};
    vector<double> ts_m;
    for (int i = 0; i < np; ++i) {
      std::vector<double> ts_i(ts0.begin(), ts0.begin() + len[i]);
      ts_m.insert(ts_m.end(), ts_i.begin(), ts_i.end());
    }
  
    vector<vector<double> > y0_m (np, y0);
    vector<vector<double> > theta_m (np, theta);
    vector<vector<double> > x_r_m (np, x_r);
    vector<vector<int> > x_i_m (np, x_i);

    size_t integ_id = integrator_id<PMXCvodesFwdSystem<TwoCptNeutModelODE, double, double, double, CV_BDF, AD>>::value;
    pmx_ode_group_mpi_functor fdyn(0);
    MatrixXd res = load.master(f, integ_id, y0_m, t0, len, ts_m, theta_m, x_r_m, x_i_m, 1.e-8, 1.e-8, 10000);
    int ic = 0;
    for (int i = 0; i < np; ++i) {
      std::vector<double> ts_i(ts0.begin(), ts0.begin() + len[i]);
      MatrixXd sol = torsten::to_matrix(pmx_integrate_ode_bdf(f, y0, t0, ts_i, theta, x_r, x_i, 0, 1.e-8, 1.e-8, 10000));
      MatrixXd res_i = res.block(0, ic, y0.size(), len[i]);
      torsten::test::test_val(res_i, sol);
      ic += len[i];
    }
    load.kill_slaves();
  }
}

TEST_F(TorstenOdeTest_neutropenia, mpi_dynamic_load_master_slave_ts_par_multiple_non_uniform_work) {
  torsten::mpi::Envionment::init();

  torsten::mpi::Communicator pmx_comm(torsten::mpi::Session<NUM_TORSTEN_COMM>::env, MPI_COMM_WORLD);
  std::vector<double> ts0 {ts};

  if (pmx_comm.rank > 0) {
    torsten::mpi::PMXDynamicLoad<TORSTEN_MPI_DYN_SLAVE> load(pmx_comm);
    load.slave();
  } else {
    torsten::mpi::PMXDynamicLoad<TORSTEN_MPI_DYN_MASTER> load(pmx_comm);
    vector<int> len{19, 17, 10, 14, 18, 17, 16};
    const int np = len.size();
    vector<var> ts_m;
    for (int i = 0; i < np; ++i) {
      std::vector<var> ts_i(ts0.begin(), ts0.begin() + len[i]);
      ts_m.insert(ts_m.end(), ts_i.begin(), ts_i.end());
    }

    vector<vector<double> > y0_m (np, y0);
    vector<vector<double> > theta_m (np, theta);
    vector<vector<double> > x_r_m (np, x_r);
    vector<vector<int> > x_i_m (np, x_i);

    {
      using Ode = PMXCvodesFwdSystem<TwoCptNeutModelODE, var, double, double, CV_BDF, AD>;
      size_t integ_id = integrator_id<Ode>::value;
      Matrix<var, -1, -1> res = load.master(f, integ_id, y0_m, t0, len, ts_m, theta_m, x_r_m, x_i_m, 1.e-8, 1.e-8, 10000);
      std::vector<var>::const_iterator its = ts_m.begin();
      int ic = 0;
      for (int i = 0; i < np; ++i) {
        std::vector<var> ts1(its, its + len[i]);
        std::vector<var> ts2(ts0.begin(), ts0.begin() + len[i]);
        Matrix<var, -1, -1> sol = torsten::to_matrix(pmx_integrate_ode_bdf(f, y0, t0, ts2, theta, x_r, x_i, 0, 1.e-8, 1.e-8, 10000));
        Matrix<var, -1, -1> res_i = res.block(0, ic, y0.size(), len[i]);
        torsten::test::test_grad(ts1, ts2, res_i, sol);
        ic += len[i];
        its += len[i];
      }
    }

    {
      using Ode = PMXCvodesFwdSystem<TwoCptNeutModelODE, var, double, double, CV_ADAMS, AD>;
      size_t integ_id = integrator_id<Ode>::value;
      Matrix<var, -1, -1> res = load.master(f, integ_id, y0_m, t0, len, ts_m, theta_m, x_r_m, x_i_m, 1.e-8, 1.e-8, 10000);
      std::vector<var>::const_iterator its = ts_m.begin();
      int ic = 0;
      for (int i = 0; i < np; ++i) {
        std::vector<var> ts1(its, its + len[i]);
        std::vector<var> ts2(ts0.begin(), ts0.begin() + len[i]);
        Matrix<var, -1, -1> sol = torsten::to_matrix(pmx_integrate_ode_adams(f, y0, t0, ts2, theta, x_r, x_i, 0, 1.e-8, 1.e-8, 10000));
        Matrix<var, -1, -1> res_i = res.block(0, ic, y0.size(), len[i]);
        torsten::test::test_grad(ts1, ts2, res_i, sol);
        ic += len[i];
        its += len[i];
      }
    }

    {
      using scheme_t = boost::numeric::odeint::runge_kutta_dopri5<std::vector<double>, double, std::vector<double>, double>;
      using Ode = dsolve::PMXOdeintSystem<TwoCptNeutModelODE, var, double, double>;
      size_t integ_id = integrator_id<Ode>::value;
      Matrix<var, -1, -1> res = load.master(f, integ_id, y0_m, t0, len, ts_m, theta_m, x_r_m, x_i_m, 1.e-8, 1.e-8, 10000);
      std::vector<var>::const_iterator its = ts_m.begin();
      int ic = 0;
      for (int i = 0; i < np; ++i) {
        std::vector<var> ts1(its, its + len[i]);
        std::vector<var> ts2(ts0.begin(), ts0.begin() + len[i]);
        Matrix<var, -1, -1> sol = torsten::to_matrix(pmx_integrate_ode_rk45(f, y0, t0, ts2, theta, x_r, x_i, 0, 1.e-8, 1.e-8, 10000));
        Matrix<var, -1, -1> res_i = res.block(0, ic, y0.size(), len[i]);
        torsten::test::test_grad(ts1, ts2, res_i, sol);
        ic += len[i];
        its += len[i];
      }
    }
    load.kill_slaves();
  }
}

TEST_F(TorstenOdeTest_neutropenia, mpi_dynamic_load_master_slave_theta_par_multiple_non_uniform_work) {
  torsten::mpi::Envionment::init();

  torsten::mpi::Communicator pmx_comm(torsten::mpi::Session<NUM_TORSTEN_COMM>::env, MPI_COMM_WORLD);
  std::vector<double> ts0 {ts};

  if (pmx_comm.rank > 0) {
    torsten::mpi::PMXDynamicLoad<TORSTEN_MPI_DYN_SLAVE> load(pmx_comm);
    load.slave();
  } else {
    torsten::mpi::PMXDynamicLoad<TORSTEN_MPI_DYN_MASTER> load(pmx_comm);
    vector<int> len{19, 17, 10, 14, 18, 17, 16};
    const int np = len.size();
    vector<double> ts_m;
    for (int i = 0; i < np; ++i) {
      std::vector<double> ts_i(ts0.begin(), ts0.begin() + len[i]);
      ts_m.insert(ts_m.end(), ts_i.begin(), ts_i.end());
    }

    vector<vector<double> > y0_m (np, y0);
    vector<vector<var> > theta_m (np);
    vector<vector<double> > x_r_m (np, x_r);
    vector<vector<int> > x_i_m (np, x_i);

    // perturb theta_m
    theta_m[0] = std::vector<var>{10, 15, 35, 105, 2, 125, 5, 0.17, 2.0e-4};
    theta_m[1] = std::vector<var>{9, 14, 34, 104, 2, 125, 5, 0.17, 2.0e-4};
    theta_m[2] = std::vector<var>{10, 14, 34, 104, 2, 125, 5, 0.17, 2.0e-4};
    theta_m[3] = std::vector<var>{10, 14, 34, 101, 2, 124, 4, 0.17, 2.0e-4};
    theta_m[4] = std::vector<var>{11, 14, 34, 101, 2, 124, 4, 0.17, 2.0e-4};
    theta_m[5] = std::vector<var>{11, 14, 35, 101, 2, 124, 4, 0.17, 2.0e-4};
    theta_m[6] = std::vector<var>{10, 14, 35, 103, 2, 120, 4, 0.17, 2.0e-4};
    
    {
      using Ode = PMXCvodesFwdSystem<TwoCptNeutModelODE, double, double, var, CV_BDF, AD>;
      size_t integ_id = integrator_id<Ode>::value;
      Matrix<var, -1, -1> res = load.master(f, integ_id, y0_m, t0, len, ts_m, theta_m, x_r_m, x_i_m, 1.e-8, 1.e-8, 10000);
      std::vector<double>::const_iterator its = ts_m.begin();
      int ic = 0;
      for (int i = 0; i < np; ++i) {
        std::vector<double> ts_i(its, its + len[i]);
        Matrix<var, -1, -1> sol = torsten::to_matrix(pmx_integrate_ode_bdf(f, y0, t0, ts_i, theta_m[i], x_r, x_i, 0, 1.e-8, 1.e-8, 10000));
        Matrix<var, -1, -1> res_i = res.block(0, ic, y0.size(), len[i]);
        torsten::test::test_grad(theta_m[i], theta_m[i], res_i, sol);
        ic += len[i];
        its += len[i];
      }
    }

    {
      using Ode = PMXCvodesFwdSystem<TwoCptNeutModelODE, double, double, var, CV_ADAMS, AD>;
      size_t integ_id = integrator_id<Ode>::value;
      Matrix<var, -1, -1> res = load.master(f, integ_id, y0_m, t0, len, ts_m, theta_m, x_r_m, x_i_m, 1.e-8, 1.e-8, 10000);
      std::vector<double>::const_iterator its = ts_m.begin();
      int ic = 0;
      for (int i = 0; i < np; ++i) {
        std::vector<double> ts_i(its, its + len[i]);
        Matrix<var, -1, -1> sol = torsten::to_matrix(pmx_integrate_ode_adams(f, y0, t0, ts_i, theta_m[i], x_r, x_i, 0, 1.e-8, 1.e-8, 10000));
        Matrix<var, -1, -1> res_i = res.block(0, ic, y0.size(), len[i]);
        torsten::test::test_grad(theta_m[i], theta_m[i], res_i, sol);
        ic += len[i];
        its += len[i];
      }
    }

    {
      using scheme_t = boost::numeric::odeint::runge_kutta_dopri5<std::vector<double>, double, std::vector<double>, double>;
      using Ode = dsolve::PMXOdeintSystem<TwoCptNeutModelODE, double, double, var>;
      size_t integ_id = integrator_id<Ode>::value;
      Matrix<var, -1, -1> res = load.master(f, integ_id, y0_m, t0, len, ts_m, theta_m, x_r_m, x_i_m, 1.e-8, 1.e-8, 10000);
      std::vector<double>::const_iterator its = ts_m.begin();
      int ic = 0;
      for (int i = 0; i < np; ++i) {
        std::vector<double> ts_i(its, its + len[i]);
        Matrix<var, -1, -1> sol = torsten::to_matrix(pmx_integrate_ode_rk45(f, y0, t0, ts_i, theta_m[i], x_r, x_i, 0, 1.e-8, 1.e-8, 10000));
        Matrix<var, -1, -1> res_i = res.block(0, ic, y0.size(), len[i]);
        torsten::test::test_grad(theta_m[i], theta_m[i], res_i, sol);
        ic += len[i];
        its += len[i];
      }
    }
    load.kill_slaves();
  }
}

TEST_F(TorstenOdeTest_neutropenia, mpi_dynamic_load_master_slave_y0_par_multiple_non_uniform_work) {
  torsten::mpi::Envionment::init();

  torsten::mpi::Communicator pmx_comm(torsten::mpi::Session<NUM_TORSTEN_COMM>::env, MPI_COMM_WORLD);
  std::vector<double> ts0 {ts};

  if (pmx_comm.rank > 0) {
    torsten::mpi::PMXDynamicLoad<TORSTEN_MPI_DYN_SLAVE> load(pmx_comm);
    load.slave();
  } else {
    torsten::mpi::PMXDynamicLoad<TORSTEN_MPI_DYN_MASTER> load(pmx_comm);
    vector<int> len{19, 17, 10, 14, 18, 17, 16};
    const int np = len.size();
    vector<double> ts_m;
    for (int i = 0; i < np; ++i) {
      std::vector<double> ts_i(ts0.begin(), ts0.begin() + len[i]);
      ts_m.insert(ts_m.end(), ts_i.begin(), ts_i.end());
    }

    vector<vector<var> > y0_m (np);
    vector<vector<double> > theta_m (np);
    vector<vector<double> > x_r_m (np, x_r);
    vector<vector<int> > x_i_m (np, x_i);

    // perturb theta_m
    theta_m[0] = std::vector<double>{10, 15, 35, 105, 2, 125, 5, 0.17, 2.0e-4};
    theta_m[1] = std::vector<double>{9, 14, 34, 104, 2, 125, 5, 0.17, 2.0e-4};
    theta_m[2] = std::vector<double>{10, 14, 34, 104, 2, 125, 5, 0.17, 2.0e-4};
    theta_m[3] = std::vector<double>{10, 14, 34, 101, 2, 124, 4, 0.17, 2.0e-4};
    theta_m[4] = std::vector<double>{11, 14, 34, 101, 2, 124, 4, 0.17, 2.0e-4};
    theta_m[5] = std::vector<double>{11, 14, 35, 101, 2, 124, 4, 0.17, 2.0e-4};
    theta_m[6] = std::vector<double>{10, 14, 35, 103, 2, 120, 4, 0.17, 2.0e-4};
    
    y0_m[0] = std::vector<var>{100.0, 10.0, 10.0, 10.0, 10.0, 10.0, 10.0, 10.0};
    y0_m[1] = std::vector<var>{88.0, 10.0, 10.0, 10.0, 10.0, 10.0, 10.0, 10.0};
    y0_m[2] = std::vector<var>{88.0, 12.0, 12.0, 13.0, 10.0, 10.0, 10.0, 10.0};
    y0_m[3] = std::vector<var>{88.0, 12.0, 12.0, 13.0, 9.0, 11.0, 9.0, 10.0};
    y0_m[4] = std::vector<var>{95.0, 10.0, 10.0, 10.0, 10.0, 10.0, 10.0, 10.0};
    y0_m[5] = std::vector<var>{95.0, 12.0, 12.0, 13.0, 10.0, 10.0, 10.0, 10.0};
    y0_m[6] = std::vector<var>{106.0, 12.0, 12.0, 13.0, 9.0, 11.0, 9.0, 10.0};

    {
      using Ode = PMXCvodesFwdSystem<TwoCptNeutModelODE, double, var, double, CV_BDF, AD>;
      size_t integ_id = integrator_id<Ode>::value;
      Matrix<var, -1, -1> res = load.master(f, integ_id, y0_m, t0, len, ts_m, theta_m, x_r_m, x_i_m, 1.e-8, 1.e-8, 10000);
      std::vector<double>::const_iterator its = ts_m.begin();
      int ic = 0;
      for (int i = 0; i < np; ++i) {
        std::vector<double> ts_i(its, its + len[i]);
        Matrix<var, -1, -1> sol = torsten::to_matrix(pmx_integrate_ode_bdf(f, y0_m[i], t0, ts_i, theta_m[i], x_r, x_i, 0, 1.e-8, 1.e-8, 10000));
        Matrix<var, -1, -1> res_i = res.block(0, ic, y0.size(), len[i]);
        torsten::test::test_grad(y0_m[i], y0_m[i], res_i, sol);
        ic += len[i];
        its += len[i];
      }
    }

    {
      using Ode = PMXCvodesFwdSystem<TwoCptNeutModelODE, double, var, double, CV_ADAMS, AD>;
      size_t integ_id = integrator_id<Ode>::value;
      Matrix<var, -1, -1> res = load.master(f, integ_id, y0_m, t0, len, ts_m, theta_m, x_r_m, x_i_m, 1.e-8, 1.e-8, 10000);
      std::vector<double>::const_iterator its = ts_m.begin();
      int ic = 0;
      for (int i = 0; i < np; ++i) {
        std::vector<double> ts_i(its, its + len[i]);
        Matrix<var, -1, -1> sol = torsten::to_matrix(pmx_integrate_ode_adams(f, y0_m[i], t0, ts_i, theta_m[i], x_r, x_i, 0, 1.e-8, 1.e-8, 10000));
        Matrix<var, -1, -1> res_i = res.block(0, ic, y0.size(), len[i]);
        torsten::test::test_grad(y0_m[i], y0_m[i], res_i, sol);
        ic += len[i];
        its += len[i];
      }
    }

    {
      using scheme_t = boost::numeric::odeint::runge_kutta_dopri5<std::vector<double>, double, std::vector<double>, double>;
      using Ode = dsolve::PMXOdeintSystem<TwoCptNeutModelODE, double, var, double>;
      size_t integ_id = integrator_id<Ode>::value;
      Matrix<var, -1, -1> res = load.master(f, integ_id, y0_m, t0, len, ts_m, theta_m, x_r_m, x_i_m, 1.e-8, 1.e-8, 10000);
      std::vector<double>::const_iterator its = ts_m.begin();
      int ic = 0;
      for (int i = 0; i < np; ++i) {
        std::vector<double> ts_i(its, its + len[i]);
        Matrix<var, -1, -1> sol = torsten::to_matrix(pmx_integrate_ode_rk45(f, y0_m[i], t0, ts_i, theta_m[i], x_r, x_i, 0, 1.e-8, 1.e-8, 10000));
        Matrix<var, -1, -1> res_i = res.block(0, ic, y0.size(), len[i]);
        torsten::test::test_grad(y0_m[i], y0_m[i], res_i, sol);
        ic += len[i];
        its += len[i];
      }
    }
    load.kill_slaves();
  }
}

TEST_F(TorstenOdeTest_neutropenia, mpi_dynamic_load_master_slave_y0_theta_par_multiple_non_uniform_work) {
  torsten::mpi::Envionment::init();

  torsten::mpi::Communicator pmx_comm(torsten::mpi::Session<NUM_TORSTEN_COMM>::env, MPI_COMM_WORLD);
  std::vector<double> ts0 {ts};

  if (pmx_comm.rank > 0) {
    torsten::mpi::PMXDynamicLoad<TORSTEN_MPI_DYN_SLAVE> load(pmx_comm);
    load.slave();
  } else {
    torsten::mpi::PMXDynamicLoad<TORSTEN_MPI_DYN_MASTER> load(pmx_comm);
    vector<int> len{19, 17, 10, 14, 18, 17, 16};
    const int np = len.size();
    vector<double> ts_m;
    for (int i = 0; i < np; ++i) {
      std::vector<double> ts_i(ts0.begin(), ts0.begin() + len[i]);
      ts_m.insert(ts_m.end(), ts_i.begin(), ts_i.end());
    }

    vector<vector<var> > y0_m (np);
    vector<vector<var> > theta_m (np);
    vector<vector<double> > x_r_m (np, x_r);
    vector<vector<int> > x_i_m (np, x_i);

    // perturb theta_m
    theta_m[0] = std::vector<var>{10, 15, 35, 105, 2, 125, 5, 0.17, 2.0e-4};
    theta_m[1] = std::vector<var>{9, 14, 34, 104, 2, 125, 5, 0.17, 2.0e-4};
    theta_m[2] = std::vector<var>{10, 14, 34, 104, 2, 125, 5, 0.17, 2.0e-4};
    theta_m[3] = std::vector<var>{10, 14, 34, 101, 2, 124, 4, 0.17, 2.0e-4};
    theta_m[4] = std::vector<var>{11, 14, 34, 101, 2, 124, 4, 0.17, 2.0e-4};
    theta_m[5] = std::vector<var>{11, 14, 35, 101, 2, 124, 4, 0.17, 2.0e-4};
    theta_m[6] = std::vector<var>{10, 14, 35, 103, 2, 120, 4, 0.17, 2.0e-4};
    
    y0_m[0] = std::vector<var>{100.0, 10.0, 10.0, 10.0, 10.0, 10.0, 10.0, 10.0};
    y0_m[1] = std::vector<var>{88.0, 10.0, 10.0, 10.0, 10.0, 10.0, 10.0, 10.0};
    y0_m[2] = std::vector<var>{88.0, 12.0, 12.0, 13.0, 10.0, 10.0, 10.0, 10.0};
    y0_m[3] = std::vector<var>{88.0, 12.0, 12.0, 13.0, 9.0, 11.0, 9.0, 10.0};
    y0_m[4] = std::vector<var>{95.0, 10.0, 10.0, 10.0, 10.0, 10.0, 10.0, 10.0};
    y0_m[5] = std::vector<var>{95.0, 12.0, 12.0, 13.0, 10.0, 10.0, 10.0, 10.0};
    y0_m[6] = std::vector<var>{106.0, 12.0, 12.0, 13.0, 9.0, 11.0, 9.0, 10.0};

    {
      using Ode = PMXCvodesFwdSystem<TwoCptNeutModelODE, double, var, var, CV_ADAMS, AD>;
      size_t integ_id = integrator_id<Ode>::value;
      Matrix<var, -1, -1> res = load.master(f, integ_id, y0_m, t0, len, ts_m, theta_m, x_r_m, x_i_m, 1.e-8, 1.e-8, 10000);
      std::vector<double>::const_iterator its = ts_m.begin();
      int ic = 0;
      for (int i = 0; i < np; ++i) {
        std::vector<double> ts_i(its, its + len[i]);
        Matrix<var, -1, -1> sol = torsten::to_matrix(pmx_integrate_ode_adams(f, y0_m[i], t0, ts_i, theta_m[i], x_r, x_i, 0, 1.e-8, 1.e-8, 10000));
        Matrix<var, -1, -1> res_i = res.block(0, ic, y0.size(), len[i]);
        std::vector<var> pars;
        pars.insert(pars.end(), y0_m[i].begin(), y0_m[i].end());
        pars.insert(pars.end(), theta_m[i].begin(), theta_m[i].end());
        torsten::test::test_grad(pars, pars, res_i, sol);
        ic += len[i];
        its += len[i];
      }
    }

    {
      using Ode = PMXCvodesFwdSystem<TwoCptNeutModelODE, double, var, var, CV_BDF, AD>;
      size_t integ_id = integrator_id<Ode>::value;
      Matrix<var, -1, -1> res = load.master(f, integ_id, y0_m, t0, len, ts_m, theta_m, x_r_m, x_i_m, 1.e-8, 1.e-8, 10000);
      std::vector<double>::const_iterator its = ts_m.begin();
      int ic = 0;
      for (int i = 0; i < np; ++i) {
        std::vector<double> ts_i(its, its + len[i]);
        Matrix<var, -1, -1> sol = torsten::to_matrix(pmx_integrate_ode_bdf(f, y0_m[i], t0, ts_i, theta_m[i], x_r, x_i, 0, 1.e-8, 1.e-8, 10000));
        Matrix<var, -1, -1> res_i = res.block(0, ic, y0.size(), len[i]);
        std::vector<var> pars;
        pars.insert(pars.end(), y0_m[i].begin(), y0_m[i].end());
        pars.insert(pars.end(), theta_m[i].begin(), theta_m[i].end());
        torsten::test::test_grad(pars, pars, res_i, sol);
        ic += len[i];
        its += len[i];
      }
    }  

    {
      using scheme_t = boost::numeric::odeint::runge_kutta_dopri5<std::vector<double>, double, std::vector<double>, double>;
      using Ode = dsolve::PMXOdeintSystem<TwoCptNeutModelODE, double, var, var>;
      size_t integ_id = integrator_id<Ode>::value;
      Matrix<var, -1, -1> res = load.master(f, integ_id, y0_m, t0, len, ts_m, theta_m, x_r_m, x_i_m, 1.e-8, 1.e-8, 10000);
      std::vector<double>::const_iterator its = ts_m.begin();
      int ic = 0;
      for (int i = 0; i < np; ++i) {
        std::vector<double> ts_i(its, its + len[i]);
        Matrix<var, -1, -1> sol = torsten::to_matrix(pmx_integrate_ode_rk45(f, y0_m[i], t0, ts_i, theta_m[i], x_r, x_i, 0, 1.e-8, 1.e-8, 10000));
        Matrix<var, -1, -1> res_i = res.block(0, ic, y0.size(), len[i]);
        std::vector<var> pars;
        pars.insert(pars.end(), y0_m[i].begin(), y0_m[i].end());
        pars.insert(pars.end(), theta_m[i].begin(), theta_m[i].end());
        torsten::test::test_grad(pars, pars, res_i, sol);
        ic += len[i];
        its += len[i];
      }
    }
    load.kill_slaves();
  }
}

TEST_F(TorstenOdeTest_neutropenia, mpi_dynamic_load_master_slave_ts_y0_theta_par_multiple_non_uniform_work) {
  torsten::mpi::Envionment::init();

  torsten::mpi::Communicator pmx_comm(torsten::mpi::Session<NUM_TORSTEN_COMM>::env, MPI_COMM_WORLD);
  std::vector<double> ts0 {ts};

  if (pmx_comm.rank > 0) {
    torsten::mpi::PMXDynamicLoad<TORSTEN_MPI_DYN_SLAVE> load(pmx_comm);
    load.slave();
  } else {
    torsten::mpi::PMXDynamicLoad<TORSTEN_MPI_DYN_MASTER> load(pmx_comm);
    vector<int> len{19, 17, 10, 14, 18, 17, 16};
    const int np = len.size();
    vector<var> ts_m;
    for (int i = 0; i < np; ++i) {
      std::vector<var> ts_i(ts0.begin(), ts0.begin() + len[i]);
      ts_m.insert(ts_m.end(), ts_i.begin(), ts_i.end());
    }

    vector<vector<var> > y0_m (np);
    vector<vector<var> > theta_m (np);
    vector<vector<double> > x_r_m (np, x_r);
    vector<vector<int> > x_i_m (np, x_i);

    // perturb theta_m
    theta_m[0] = std::vector<var>{10, 15, 35, 105, 2, 125, 5, 0.17, 2.0e-4};
    theta_m[1] = std::vector<var>{9, 14, 34, 104, 2, 125, 5, 0.17, 2.0e-4};
    theta_m[2] = std::vector<var>{10, 14, 34, 104, 2, 125, 5, 0.17, 2.0e-4};
    theta_m[3] = std::vector<var>{10, 14, 34, 101, 2, 124, 4, 0.17, 2.0e-4};
    theta_m[4] = std::vector<var>{11, 14, 34, 101, 2, 124, 4, 0.17, 2.0e-4};
    theta_m[5] = std::vector<var>{11, 14, 35, 101, 2, 124, 4, 0.17, 2.0e-4};
    theta_m[6] = std::vector<var>{10, 14, 35, 103, 2, 120, 4, 0.17, 2.0e-4};
    
    y0_m[0] = std::vector<var>{100.0, 10.0, 10.0, 10.0, 10.0, 10.0, 10.0, 10.0};
    y0_m[1] = std::vector<var>{88.0, 10.0, 10.0, 10.0, 10.0, 10.0, 10.0, 10.0};
    y0_m[2] = std::vector<var>{88.0, 12.0, 12.0, 13.0, 10.0, 10.0, 10.0, 10.0};
    y0_m[3] = std::vector<var>{88.0, 12.0, 12.0, 13.0, 9.0, 11.0, 9.0, 10.0};
    y0_m[4] = std::vector<var>{95.0, 10.0, 10.0, 10.0, 10.0, 10.0, 10.0, 10.0};
    y0_m[5] = std::vector<var>{95.0, 12.0, 12.0, 13.0, 10.0, 10.0, 10.0, 10.0};
    y0_m[6] = std::vector<var>{106.0, 12.0, 12.0, 13.0, 9.0, 11.0, 9.0, 10.0};

    {
      using Ode = PMXCvodesFwdSystem<TwoCptNeutModelODE, var, var, var, CV_BDF, AD>;
      size_t integ_id = integrator_id<Ode>::value;
      Matrix<var, -1, -1> res = load.master(f, integ_id, y0_m, t0, len, ts_m, theta_m, x_r_m, x_i_m, 1.e-8, 1.e-8, 10000);
      std::vector<var>::const_iterator its = ts_m.begin();
      int ic = 0;
      for (int i = 0; i < np; ++i) {
        std::vector<var> ts_i(its, its + len[i]);
        Matrix<var, -1, -1> sol = torsten::to_matrix(pmx_integrate_ode_bdf(f, y0_m[i], t0, ts_i, theta_m[i], x_r, x_i, 0, 1.e-8, 1.e-8, 10000));
        Matrix<var, -1, -1> res_i = res.block(0, ic, y0.size(), len[i]);
        std::vector<var> pars;
        pars.insert(pars.end(), y0_m[i].begin(), y0_m[i].end());
        pars.insert(pars.end(), theta_m[i].begin(), theta_m[i].end());
        pars.insert(pars.end(), ts_i.begin(), ts_i.end());
        torsten::test::test_grad(pars, pars, res_i, sol);
        ic += len[i];
        its += len[i];
      }
    }

    {
      using Ode = PMXCvodesFwdSystem<TwoCptNeutModelODE, var, var, var, CV_ADAMS, AD>;
      size_t integ_id = integrator_id<Ode>::value;
      Matrix<var, -1, -1> res = load.master(f, integ_id, y0_m, t0, len, ts_m, theta_m, x_r_m, x_i_m, 1.e-8, 1.e-8, 10000);
      std::vector<var>::const_iterator its = ts_m.begin();
      int ic = 0;
      for (int i = 0; i < np; ++i) {
        std::vector<var> ts_i(its, its + len[i]);
        Matrix<var, -1, -1> sol = torsten::to_matrix(pmx_integrate_ode_adams(f, y0_m[i], t0, ts_i, theta_m[i], x_r, x_i, 0, 1.e-8, 1.e-8, 10000));
        Matrix<var, -1, -1> res_i = res.block(0, ic, y0.size(), len[i]);
        std::vector<var> pars;
        pars.insert(pars.end(), y0_m[i].begin(), y0_m[i].end());
        pars.insert(pars.end(), theta_m[i].begin(), theta_m[i].end());
        pars.insert(pars.end(), ts_i.begin(), ts_i.end());
        torsten::test::test_grad(pars, pars, res_i, sol);
        ic += len[i];
        its += len[i];
      }
    }

    {
      using scheme_t = boost::numeric::odeint::runge_kutta_dopri5<std::vector<double>, double, std::vector<double>, double>;
      using Ode = dsolve::PMXOdeintSystem<TwoCptNeutModelODE, var, var, var>;
      size_t integ_id = integrator_id<Ode>::value;
      Matrix<var, -1, -1> res = load.master(f, integ_id, y0_m, t0, len, ts_m, theta_m, x_r_m, x_i_m, 1.e-8, 1.e-8, 10000);
      std::vector<var>::const_iterator its = ts_m.begin();
      int ic = 0;
      for (int i = 0; i < np; ++i) {
        std::vector<var> ts_i(its, its + len[i]);
        Matrix<var, -1, -1> sol = torsten::to_matrix(pmx_integrate_ode_rk45(f, y0_m[i], t0, ts_i, theta_m[i], x_r, x_i, 0, 1.e-8, 1.e-8, 10000));
        Matrix<var, -1, -1> res_i = res.block(0, ic, y0.size(), len[i]);
        std::vector<var> pars;
        pars.insert(pars.end(), y0_m[i].begin(), y0_m[i].end());
        pars.insert(pars.end(), theta_m[i].begin(), theta_m[i].end());
        pars.insert(pars.end(), ts_i.begin(), ts_i.end());
        torsten::test::test_grad(pars, pars, res_i, sol);
        ic += len[i];
        its += len[i];
      }
    }
    load.kill_slaves();
  }
}

TEST_F(TorstenOdeTest_neutropenia, mpi_dynamic_load_master_slave_ts_y0_par_multiple_non_uniform_work) {
  torsten::mpi::Envionment::init();

  torsten::mpi::Communicator pmx_comm(torsten::mpi::Session<NUM_TORSTEN_COMM>::env, MPI_COMM_WORLD);
  std::vector<double> ts0 {ts};

  if (pmx_comm.rank > 0) {
    torsten::mpi::PMXDynamicLoad<TORSTEN_MPI_DYN_SLAVE> load(pmx_comm);
    load.slave();
  } else {
    torsten::mpi::PMXDynamicLoad<TORSTEN_MPI_DYN_MASTER> load(pmx_comm);
    vector<int> len{19, 17, 10, 14, 18, 17, 16};
    const int np = len.size();
    vector<var> ts_m;
    for (int i = 0; i < np; ++i) {
      std::vector<var> ts_i(ts0.begin(), ts0.begin() + len[i]);
      ts_m.insert(ts_m.end(), ts_i.begin(), ts_i.end());
    }

    vector<vector<var> > y0_m (np, stan::math::to_var(y0));
    vector<vector<double> > theta_m (np);
    vector<vector<double> > x_r_m (np, x_r);
    vector<vector<int> > x_i_m (np, x_i);

    // perturb theta_m
    theta_m[0] = std::vector<double>{10, 15, 35, 105, 2, 125, 5, 0.17, 2.0e-4};
    theta_m[1] = std::vector<double>{9, 14, 34, 104, 2, 125, 5, 0.17, 2.0e-4};
    theta_m[2] = std::vector<double>{10, 14, 34, 104, 2, 125, 5, 0.17, 2.0e-4};
    theta_m[3] = std::vector<double>{10, 14, 34, 101, 2, 124, 4, 0.17, 2.0e-4};
    theta_m[4] = std::vector<double>{11, 14, 34, 101, 2, 124, 4, 0.17, 2.0e-4};
    theta_m[5] = std::vector<double>{11, 14, 35, 101, 2, 124, 4, 0.17, 2.0e-4};
    theta_m[6] = std::vector<double>{10, 14, 35, 103, 2, 120, 4, 0.17, 2.0e-4};
    
    {
      using Ode = PMXCvodesFwdSystem<TwoCptNeutModelODE, var, var, double, CV_ADAMS, AD>;
      size_t integ_id = integrator_id<Ode>::value;
      Matrix<var, -1, -1> res = load.master(f, integ_id, y0_m, t0, len, ts_m, theta_m, x_r_m, x_i_m, 1.e-8, 1.e-8, 10000);
      std::vector<var>::const_iterator its = ts_m.begin();
      int ic = 0;
      for (int i = 0; i < np; ++i) {
        std::vector<var> ts_i(its, its + len[i]);
        Matrix<var, -1, -1> sol = torsten::to_matrix(pmx_integrate_ode_adams(f, y0_m[i], t0, ts_i, theta_m[i], x_r, x_i, 0, 1.e-8, 1.e-8, 10000));
        Matrix<var, -1, -1> res_i = res.block(0, ic, y0.size(), len[i]);
        std::vector<var> pars;
        pars.insert(pars.end(), y0_m[i].begin(), y0_m[i].end());
        pars.insert(pars.end(), ts_i.begin(), ts_i.end());
        torsten::test::test_grad(pars, pars, res_i, sol);
        ic += len[i];
        its += len[i];
      }
    }

    {
      using Ode = PMXCvodesFwdSystem<TwoCptNeutModelODE, var, var, double, CV_BDF, AD>;
      size_t integ_id = integrator_id<Ode>::value;
      Matrix<var, -1, -1> res = load.master(f, integ_id, y0_m, t0, len, ts_m, theta_m, x_r_m, x_i_m, 1.e-8, 1.e-8, 10000);
      std::vector<var>::const_iterator its = ts_m.begin();
      int ic = 0;
      for (int i = 0; i < np; ++i) {
        std::vector<var> ts_i(its, its + len[i]);
        Matrix<var, -1, -1> sol = torsten::to_matrix(pmx_integrate_ode_bdf(f, y0_m[i], t0, ts_i, theta_m[i], x_r, x_i, 0, 1.e-8, 1.e-8, 10000));
        Matrix<var, -1, -1> res_i = res.block(0, ic, y0.size(), len[i]);
        std::vector<var> pars;
        pars.insert(pars.end(), y0_m[i].begin(), y0_m[i].end());
        pars.insert(pars.end(), ts_i.begin(), ts_i.end());
        torsten::test::test_grad(pars, pars, res_i, sol);
        ic += len[i];
        its += len[i];
      }
    }

    {
      using scheme_t = boost::numeric::odeint::runge_kutta_dopri5<std::vector<double>, double, std::vector<double>, double>;
      using Ode = dsolve::PMXOdeintSystem<TwoCptNeutModelODE, var, var, double>;
      size_t integ_id = integrator_id<Ode>::value;
      Matrix<var, -1, -1> res = load.master(f, integ_id, y0_m, t0, len, ts_m, theta_m, x_r_m, x_i_m, 1.e-8, 1.e-8, 10000);
      std::vector<var>::const_iterator its = ts_m.begin();
      int ic = 0;
      for (int i = 0; i < np; ++i) {
        std::vector<var> ts_i(its, its + len[i]);
        Matrix<var, -1, -1> sol = torsten::to_matrix(pmx_integrate_ode_rk45(f, y0_m[i], t0, ts_i, theta_m[i], x_r, x_i, 0, 1.e-8, 1.e-8, 10000));
        Matrix<var, -1, -1> res_i = res.block(0, ic, y0.size(), len[i]);
        std::vector<var> pars;
        pars.insert(pars.end(), y0_m[i].begin(), y0_m[i].end());
        pars.insert(pars.end(), ts_i.begin(), ts_i.end());
        torsten::test::test_grad(pars, pars, res_i, sol);
        ic += len[i];
        its += len[i];
      }
    }
    load.kill_slaves();
  }
}


#endif
