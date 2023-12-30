#include <gtest/gtest.h>
#include <stan/math/rev/core.hpp>
#include <stan/math/rev.hpp>
#include <stan/math/torsten/dsolve/pmx_ode_system.hpp>
#include <stan/math/torsten/dsolve/pmx_odeint_integrator.hpp>
#include <stan/math/torsten/dsolve/pmx_integrate_ode_rk45.hpp>
#include <stan/math/torsten/dsolve/pmx_integrate_ode_bdf.hpp>
#include <stan/math/torsten/dsolve/pmx_ode_rk45.hpp>
#include <stan/math/torsten/dsolve/pmx_ode_bdf.hpp>
#include <stan/math/torsten/dsolve/pmx_ode_adams.hpp>
#include <stan/math/torsten/dsolve/pmx_ode_ckrk.hpp>
#include <stan/math/torsten/test/unit/pmx_ode_test_fixture.hpp>
#include <sundials/sundials_context.h>
#include <boost/numeric/odeint/external/eigen/eigen_algebra.hpp>
#include <boost/fusion/adapted/std_tuple.hpp>
#include <boost/fusion/include/algorithm.hpp>

TEST_F(TorstenOdeTest_sho, variadic_ode_system_odeint) {
  using torsten::dsolve::PMXVariadicOdeSystem;
  using stan::math::value_of;
  using stan::math::to_var;
  
  y0[1] = 0.5;
  y0_vec[1] = 0.5;
  Eigen::Matrix<stan::math::var, -1, 1> y0_var(to_var(y0_vec));
  std::vector<stan::math::var> theta_var(to_var(theta));

  PMXVariadicOdeSystem<harm_osc_ode_fun_eigen, double, stan::math::var,
                       std::vector<stan::math::var>, std::vector<double>, std::vector<int>>
    ode(f_eigen, 0.0, ts, y0_var, nullptr, theta_var, x_r, x_i);

  Eigen::VectorXd dydt = f_eigen(ts[0], y0_vec, nullptr, theta, x_r, x_i);

  Eigen::VectorXd ydot = ode.dbl_rhs_impl(ts[0], y0_vec);
  EXPECT_FLOAT_EQ(ydot[0], dydt[0]);
  EXPECT_FLOAT_EQ(ydot[1], dydt[1]);

  torsten::dsolve::PMXOdeSystem<harm_osc_ode_fun, double, stan::math::var, stan::math::var>
    ode0(f, 0.0, ts, to_var(y0), theta_var, x_r, x_i, nullptr);
  std::vector<double> dydt_vec(ode.system_size), dydt_vec2(ode.system_size);
  ode0(ode0.y0_fwd_system, dydt_vec, ts[0]);

  ode(ode.y0_fwd_system, dydt_vec2, ts[0]);
  EXPECT_EQ(ode.system_size, ode0.system_size);
  for (size_t i = 0; i < ode0.system_size; ++i) {
    EXPECT_FLOAT_EQ(dydt_vec2[i], dydt_vec[i]);
  }
}

TEST_F(TorstenOdeTest_sho, variadic_ode_system_cvodes) {
  using torsten::dsolve::PMXVariadicOdeSystem;
  
  using stan::math::value_of;
  using stan::math::to_var;
  
  y0[1] = 0.5;
  y0_vec[1] = 0.5;

  Eigen::Matrix<stan::math::var, -1, 1> y0_var(to_var(y0_vec));
  std::vector<stan::math::var> theta_var(to_var(theta));

  PMXVariadicOdeSystem<harm_osc_ode_fun_eigen, double, stan::math::var,
                       std::vector<stan::math::var>, std::vector<double>, std::vector<int>>
    ode(f_eigen, 0.0, ts, y0_var, nullptr, theta_var, x_r, x_i);

  Eigen::VectorXd dydt = f_eigen(ts[0], y0_vec, nullptr, theta, x_r, x_i);
  sundials::Context sundials_context_;
  N_Vector nv_y(N_VNew_Serial(2, sundials_context_));
  N_Vector ydot(N_VNew_Serial(2, sundials_context_));
  NV_Ith_S(nv_y, 0) = y0[0];
  NV_Ith_S(nv_y, 1) = y0[1];
  ode(ts[0], nv_y, ydot);
  EXPECT_FLOAT_EQ(NV_Ith_S(ydot, 0), dydt[0]);
  EXPECT_FLOAT_EQ(NV_Ith_S(ydot, 1), dydt[1]);
  N_VDestroy(nv_y);
  N_VDestroy(ydot);
}

TEST_F(TorstenOdeTest_sho, eigen_vector_rk45) {
  using torsten::dsolve::PMXVariadicOdeSystem;
  
  using stan::math::value_of;
  using torsten::dsolve::PMXOdeintIntegrator;

  {                             // data only
    auto y = torsten::pmx_ode_rk45_ctrl(f_eigen, y0_vec, t0, ts, msgs, rtol, atol, max_num_steps,
                                        theta, x_r, x_i);
    auto y_sol = torsten::pmx_integrate_ode_rk45(f, y0, t0, ts,
                                                 theta, x_r, x_i, rtol, atol, max_num_steps, msgs);
    for (size_t j = 0; j < ts.size(); ++j) {
      for (size_t i = 0; i < y0.size(); ++i) {
        EXPECT_FLOAT_EQ(y[j][i], y_sol[j][i]);
      }
    }
  }

  {                             // theat var
    std::vector<stan::math::var> theta_var(stan::math::to_var(theta));
    auto y = torsten::pmx_ode_rk45_ctrl(f_eigen, y0_vec, t0, ts, msgs, rtol, atol, max_num_steps,
                                        theta_var, x_r, x_i);
    auto y_sol = stan::math::ode_rk45_tol(f_eigen, y0_vec, t0, ts, rtol, atol, max_num_steps, msgs,
                                          theta_var, x_r, x_i);
    for (size_t j = 0; j < ts.size(); ++j) {
      torsten::test::test_grad(theta_var, y_sol[j], y[j], 1e-8, 1e-8);
    }
  }

  {                             // theat & y0 var
    std::vector<stan::math::var> theta_var(stan::math::to_var(theta));
    std::vector<stan::math::var> y0_var(stan::math::to_var(y0));
    Eigen::Matrix<stan::math::var, -1, 1> y0_vec_var(stan::math::to_vector(y0_var));

    auto y = torsten::pmx_ode_rk45_ctrl(f_eigen, y0_vec_var, t0, ts, msgs, rtol, atol, max_num_steps,
                                        theta_var, x_r, x_i);
    auto y_sol = stan::math::ode_rk45_tol(f_eigen, y0_vec_var, t0, ts, rtol, atol, max_num_steps, msgs,
                                          theta_var, x_r, x_i);

    for (size_t j = 0; j < ts.size(); ++j) {
      torsten::test::test_grad(theta_var, y_sol[j], y[j], 1e-8, 1e-8);
      torsten::test::test_grad(y0_var, y_sol[j], y[j], 1e-8, 1e-8);
    }
  }

  {                             // y0 & ts var
    std::vector<stan::math::var> y0_var(stan::math::to_var(y0));
    Eigen::Matrix<stan::math::var, -1, 1> y0_vec_var(stan::math::to_vector(y0_var));
    std::vector<stan::math::var> ts_var(stan::math::to_var(ts));

    auto y = torsten::pmx_ode_rk45_ctrl(f_eigen, y0_vec_var, t0, ts_var, msgs, rtol, atol, max_num_steps,
                                        theta, x_r, x_i);
    auto y_sol = stan::math::ode_rk45_tol(f_eigen, y0_vec_var, t0, ts_var, rtol, atol, max_num_steps, msgs,
                                          theta, x_r, x_i);

    for (size_t j = 0; j < ts.size(); ++j) {
      torsten::test::test_grad(ts_var, y_sol[j], y[j], 1e-8, 1e-8);
      torsten::test::test_grad(y0_var, y_sol[j], y[j], 1e-8, 1e-8);
    }
  }
}

TEST_F(TorstenOdeTest_chem, eigen_vector_bdf) {
  using torsten::dsolve::PMXVariadicOdeSystem;
  
  using stan::math::value_of;
  using torsten::dsolve::PMXOdeintIntegrator;

  {                             // data only
    auto y = torsten::pmx_ode_bdf_ctrl(f_eigen, y0_vec, t0, ts, msgs, rtol, atol, max_num_steps,
                                        theta, x_r, x_i);
    auto y_sol = torsten::pmx_integrate_ode_bdf(f, y0, t0, ts,
                                                 theta, x_r, x_i, rtol, atol, max_num_steps, msgs);
    for (size_t j = 0; j < ts.size(); ++j) {
      for (size_t i = 0; i < y0.size(); ++i) {
        EXPECT_FLOAT_EQ(y[j][i], y_sol[j][i]);
      }
    }
  }

  {                             // theat var
    std::vector<stan::math::var> theta_var(stan::math::to_var(theta));
    auto y = torsten::pmx_ode_bdf_ctrl(f_eigen, y0_vec, t0, ts, msgs, rtol, atol, max_num_steps,
                                        theta_var, x_r, x_i);
    auto y_sol = stan::math::ode_bdf_tol(f_eigen, y0_vec, t0, ts, rtol, atol, max_num_steps, msgs,
                                          theta_var, x_r, x_i);
    for (size_t j = 0; j < ts.size(); ++j) {
      torsten::test::test_grad(theta_var, y_sol[j], y[j], 1e-8, 1e-6);
    }
  }

  {                             // theat & y0 var
    std::vector<stan::math::var> theta_var(stan::math::to_var(theta));
    std::vector<stan::math::var> y0_var(stan::math::to_var(y0));
    Eigen::Matrix<stan::math::var, -1, 1> y0_vec_var(stan::math::to_vector(y0_var));

    auto y = torsten::pmx_ode_bdf_ctrl(f_eigen, y0_vec_var, t0, ts, msgs, rtol, atol, max_num_steps,
                                        theta_var, x_r, x_i);
    auto y_sol = stan::math::ode_bdf_tol(f_eigen, y0_vec_var, t0, ts, rtol, atol, max_num_steps, msgs,
                                          theta_var, x_r, x_i);

    for (size_t j = 0; j < ts.size(); ++j) {
      torsten::test::test_grad(theta_var, y_sol[j], y[j], 1e-8, 1e-6);
      torsten::test::test_grad(y0_var, y_sol[j], y[j], 1e-8, 1e-6);
    }
  }

  {                             // y0 & ts var
    std::vector<stan::math::var> y0_var(stan::math::to_var(y0));
    Eigen::Matrix<stan::math::var, -1, 1> y0_vec_var(stan::math::to_vector(y0_var));
    std::vector<stan::math::var> ts_var(stan::math::to_var(ts));

    auto y = torsten::pmx_ode_bdf_ctrl(f_eigen, y0_vec_var, t0, ts_var, msgs, rtol, atol, max_num_steps,
                                        theta, x_r, x_i);
    auto y_sol = stan::math::ode_bdf_tol(f_eigen, y0_vec_var, t0, ts_var, rtol, atol, max_num_steps, msgs,
                                          theta, x_r, x_i);

    for (size_t j = 0; j < ts.size(); ++j) {
      torsten::test::test_grad(ts_var, y_sol[j], y[j], 1e-8, 1e-6);
      torsten::test::test_grad(y0_var, y_sol[j], y[j], 1e-8, 1e-6);
    }
  }
}

TEST_F(TorstenOdeTest_lorenz, eigen_vector_adams) {
  using torsten::dsolve::PMXVariadicOdeSystem;
  
  using stan::math::value_of;
  using torsten::dsolve::PMXOdeintIntegrator;

  {                             // theat & y0 var
    std::vector<stan::math::var> theta_var(stan::math::to_var(theta));
    std::vector<stan::math::var> y0_var(stan::math::to_var(y0));
    Eigen::Matrix<stan::math::var, -1, 1> y0_vec_var(stan::math::to_vector(y0_var));

    auto y = torsten::pmx_ode_adams_ctrl(f_eigen, y0_vec_var, t0, ts, msgs, rtol, atol, max_num_steps,
                                        theta_var, x_r, x_i);
    auto y_sol = stan::math::ode_adams_tol(f_eigen, y0_vec_var, t0, ts, rtol, atol, max_num_steps, msgs,
                                          theta_var, x_r, x_i);

    for (size_t j = 0; j < ts.size(); ++j) {
      torsten::test::test_grad(theta_var, y_sol[j], y[j], 1e-6, 8e-6);
      torsten::test::test_grad(y0_var, y_sol[j], y[j], 1e-6, 8e-6);
    }
  }

  {                             // y0 & ts var
    std::vector<stan::math::var> y0_var(stan::math::to_var(y0));
    Eigen::Matrix<stan::math::var, -1, 1> y0_vec_var(stan::math::to_vector(y0_var));
    std::vector<stan::math::var> ts_var(stan::math::to_var(ts));

    auto y = torsten::pmx_ode_adams_ctrl(f_eigen, y0_vec_var, t0, ts_var, msgs, rtol, atol, max_num_steps,
                                        theta, x_r, x_i);
    auto y_sol = stan::math::ode_adams_tol(f_eigen, y0_vec_var, t0, ts_var, rtol, atol, max_num_steps, msgs,
                                          theta, x_r, x_i);

    for (size_t j = 0; j < ts.size(); ++j) {
      torsten::test::test_grad(ts_var, y_sol[j], y[j], 1e-6, 8e-6);
      torsten::test::test_grad(y0_var, y_sol[j], y[j], 1e-6, 8e-6);
    }
  }

  {                             // ts var
    std::vector<stan::math::var> ts_var(stan::math::to_var(ts));
    auto y = torsten::pmx_ode_adams_ctrl(f_eigen, y0_vec, t0, ts_var, msgs, rtol, atol, max_num_steps,
                                        theta, x_r, x_i);
    auto y_sol = stan::math::ode_adams_tol(f_eigen, y0_vec, t0, ts_var, rtol, atol, max_num_steps, msgs,
                                          theta, x_r, x_i);

    for (size_t j = 0; j < ts.size(); ++j) {
      torsten::test::test_grad(ts_var, y_sol[j], y[j], 1e-6, 8e-6);
    }
  }
}

TEST_F(TorstenOdeTest_neutropenia, eigen_vector_ckrk) {
  using torsten::dsolve::PMXVariadicOdeSystem;
  
  using stan::math::value_of;
  using torsten::dsolve::PMXOdeintIntegrator;

  {                             // theat & y0 var
    std::vector<stan::math::var> theta_var(stan::math::to_var(theta));
    std::vector<stan::math::var> y0_var(stan::math::to_var(y0));
    Eigen::Matrix<stan::math::var, -1, 1> y0_vec_var(stan::math::to_vector(y0_var));

    auto y = torsten::pmx_ode_ckrk_ctrl(f_eigen, y0_vec_var, t0, ts, msgs, rtol, atol, max_num_steps,
                                        theta_var, x_r, x_i);
    auto y_sol = stan::math::ode_rk45_tol(f_eigen, y0_vec_var, t0, ts, rtol, atol, max_num_steps, msgs,
                                          theta_var, x_r, x_i);

    for (size_t j = 0; j < ts.size(); ++j) {
      torsten::test::test_grad(theta_var, y_sol[j], y[j], 1e-8, 1e-8);
      torsten::test::test_grad(y0_var, y_sol[j], y[j], 1e-8, 1e-8);
    }
  }

  {                             // y0 & ts var
    std::vector<stan::math::var> y0_var(stan::math::to_var(y0));
    Eigen::Matrix<stan::math::var, -1, 1> y0_vec_var(stan::math::to_vector(y0_var));
    std::vector<stan::math::var> ts_var(stan::math::to_var(ts));

    auto y = torsten::pmx_ode_ckrk_ctrl(f_eigen, y0_vec_var, t0, ts_var, msgs, rtol, atol, max_num_steps,
                                        theta, x_r, x_i);
    auto y_sol = stan::math::ode_rk45_tol(f_eigen, y0_vec_var, t0, ts_var, rtol, atol, max_num_steps, msgs,
                                          theta, x_r, x_i);

    for (size_t j = 0; j < ts.size(); ++j) {
      torsten::test::test_grad(ts_var, y_sol[j], y[j], 1e-8, 1e-7);
      torsten::test::test_grad(y0_var, y_sol[j], y[j], 1e-8, 1e-7);
    }
  }

  {                             // ts var
    std::vector<stan::math::var> ts_var(stan::math::to_var(ts));
    auto y = torsten::pmx_ode_ckrk_ctrl(f_eigen, y0_vec, t0, ts_var, msgs, rtol, atol, max_num_steps,
                                        theta, x_r, x_i);
    auto y_sol = stan::math::ode_rk45_tol(f_eigen, y0_vec, t0, ts_var, rtol, atol, max_num_steps, msgs,
                                          theta, x_r, x_i);

    for (size_t j = 0; j < ts.size(); ++j) {
      torsten::test::test_grad(ts_var, y_sol[j], y[j], 1e-8, 1e-7);
    }
  }
}
