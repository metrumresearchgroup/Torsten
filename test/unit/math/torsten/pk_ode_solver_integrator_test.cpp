#include <stan/math.hpp>
#include <stan/math/rev/core.hpp>
#include <test/unit/math/rev/mat/fun/util.hpp>
#include <test/unit/math/torsten/pk_cpt_model_test_fixture.hpp>
#include <stan/math/torsten/pk_ode_model.hpp>
#include <test/unit/util.hpp>
#include <gtest/gtest.h>


TEST_F(TorstenCptOdeModelTest, general_ode_solver) {
  using stan::math::var;
  using stan::math::to_var;
  using refactor::PKTwoCptModel;
  using refactor::PKTwoCptODE;
  using refactor::PKODEModel;
  using refactor::PKOdeFunctorRateAdaptor;

  rate[0] = 1200;
  rate[1] = 2000;
  rate[2] = 3000;
  PKTwoCptModel<double, double, double, double> model0(t0, y0, rate, CL, Q, V2, V3, ka); // NOLINT
  std::vector<double> yvec(y0.data(), y0.data() + y0.size());
  PKOdeFunctorRateAdaptor<PKTwoCptODE, double> f1(model0.f());
  using model_t = PKODEModel<double, double, double, double, PKTwoCptODE, int>;
  model_t model(t0, y0, rate, model0.par(), model0.f(), model0.ncmt());

  Eigen::Matrix<torsten::scalar_t<model_t>, Eigen::Dynamic, 1> y;
  std::vector<std::vector<double> > y1;
  ts[0] = 20.0;
  ts.resize(1);

  PkOdeIntegrator<StanRk45> integ1(rtol, atol, max_num_steps, msgs);
  y = model.solve(ts[0], integ1);
  y1 = stan::math::integrate_ode_rk45(f1, yvec, t0, ts, model.par(), model.rate(), x_i, msgs); // NOLINT
  EXPECT_FLOAT_EQ(y(0), y1[0][0]);
  EXPECT_FLOAT_EQ(y(1), y1[0][1]);
  EXPECT_FLOAT_EQ(y(2), y1[0][2]);

  PkOdeIntegrator<StanAdams> integ2(rtol, atol, max_num_steps, msgs);
  y = model.solve(ts[0], integ2);
  y1 = stan::math::integrate_ode_adams(f1, yvec, t0, ts, model.par(), model.rate(), x_i, msgs); // NOLINT
  EXPECT_FLOAT_EQ(y(0), y1[0][0]);
  EXPECT_FLOAT_EQ(y(1), y1[0][1]);
  EXPECT_FLOAT_EQ(y(2), y1[0][2]);

  PkOdeIntegrator<StanBdf> integ3(rtol, atol, max_num_steps, msgs);
  y = model.solve(ts[0], integ3);
  y1 = stan::math::integrate_ode_bdf(f1, yvec, t0, ts, model.par(), model.rate(), x_i, msgs); // NOLINT
  EXPECT_FLOAT_EQ(y(0), y1[0][0]);
  EXPECT_FLOAT_EQ(y(1), y1[0][1]);
  EXPECT_FLOAT_EQ(y(2), y1[0][2]);

  PkOdeIntegrator<PkAdams> integ4(rtol, atol, max_num_steps, msgs);
  y = model.solve(ts[0], integ4);
  y1 = torsten::dsolve::pk_integrate_ode_adams(f1, yvec, t0, ts, model.par(), model.rate(), x_i, msgs); // NOLINT
  EXPECT_FLOAT_EQ(y(0), y1[0][0]);
  EXPECT_FLOAT_EQ(y(1), y1[0][1]);
  EXPECT_FLOAT_EQ(y(2), y1[0][2]);

  PkOdeIntegrator<PkBdf> integ5(rtol, atol, max_num_steps, msgs);
  y = model.solve(ts[0], integ5);
  y1 = torsten::dsolve::pk_integrate_ode_bdf(f1, yvec, t0, ts, model.par(), model.rate(), x_i, msgs); // NOLINT
  EXPECT_FLOAT_EQ(y(0), y1[0][0]);
  EXPECT_FLOAT_EQ(y(1), y1[0][1]);
  EXPECT_FLOAT_EQ(y(2), y1[0][2]);
}

TEST_F(TorstenCptOdeModelTest, general_ode_solver_y0) {
  using stan::math::var;
  using stan::math::to_var;
  using refactor::PKTwoCptModel;
  using refactor::PKTwoCptODE;
  using refactor::PKODEModel;
  using refactor::PKOdeFunctorRateAdaptor;

  rate[0] = 1200;
  rate[1] = 2000;
  rate[2] = 3000;
  y0[0] = 800;
  y0[1] = 0;
  y0[2] = 8000;
  PKTwoCptModel<double, double, double, double> model0(t0, y0, rate, CL, Q, V2, V3, ka); // NOLINT
  std::vector<double> yvec(y0.data(), y0.data() + y0.size());
  PKOdeFunctorRateAdaptor<PKTwoCptODE, double> f1(model0.f());
  using model_t = PKODEModel<double, double, double, double, PKTwoCptODE, int>;
  model_t model(t0, y0, rate, model0.par(), model0.f(), model0.ncmt());

  Eigen::Matrix<torsten::scalar_t<model_t>, Eigen::Dynamic, 1> y;
  std::vector<std::vector<double> > y1;
  ts[0] = 20.0;
  ts.resize(1);

  PkOdeIntegrator<StanRk45> integ1(rtol, atol, max_num_steps, msgs);
  y = model.solve(ts[0], integ1);
  y1 = stan::math::integrate_ode_rk45(f1, yvec, t0, ts, model.par(), model.rate(), x_i, msgs); // NOLINT
  EXPECT_FLOAT_EQ(y(0), y1[0][0]);
  EXPECT_FLOAT_EQ(y(1), y1[0][1]);
  EXPECT_FLOAT_EQ(y(2), y1[0][2]);

  PkOdeIntegrator<StanAdams> integ2(rtol, atol, max_num_steps, msgs);
  y = model.solve(ts[0], integ2);
  y1 = stan::math::integrate_ode_adams(f1, yvec, t0, ts, model.par(), model.rate(), x_i, msgs); // NOLINT
  EXPECT_FLOAT_EQ(y(0), y1[0][0]);
  EXPECT_FLOAT_EQ(y(1), y1[0][1]);
  EXPECT_FLOAT_EQ(y(2), y1[0][2]);

  PkOdeIntegrator<StanBdf> integ3(rtol, atol, max_num_steps, msgs);
  y = model.solve(ts[0], integ3);
  y1 = stan::math::integrate_ode_bdf(f1, yvec, t0, ts, model.par(), model.rate(), x_i, msgs); // NOLINT
  EXPECT_FLOAT_EQ(y(0), y1[0][0]);
  EXPECT_FLOAT_EQ(y(1), y1[0][1]);
  EXPECT_FLOAT_EQ(y(2), y1[0][2]);

  PkOdeIntegrator<PkAdams> integ4(rtol, atol, max_num_steps, msgs);
  y = model.solve(ts[0], integ4);
  y1 = torsten::dsolve::pk_integrate_ode_adams(f1, yvec, t0, ts, model.par(), model.rate(), x_i, msgs); // NOLINT
  EXPECT_FLOAT_EQ(y(0), y1[0][0]);
  EXPECT_FLOAT_EQ(y(1), y1[0][1]);
  EXPECT_FLOAT_EQ(y(2), y1[0][2]);

  PkOdeIntegrator<PkBdf> integ5(rtol, atol, max_num_steps, msgs);
  y = model.solve(ts[0], integ5);
  y1 = torsten::dsolve::pk_integrate_ode_bdf(f1, yvec, t0, ts, model.par(), model.rate(), x_i, msgs); // NOLINT
  EXPECT_FLOAT_EQ(y(0), y1[0][0]);
  EXPECT_FLOAT_EQ(y(1), y1[0][1]);
  EXPECT_FLOAT_EQ(y(2), y1[0][2]);
}

TEST_F(TorstenCptOdeModelTest, general_ode_solver_par_sens) {
  using stan::math::var;
  using stan::math::to_var;
  using refactor::PKTwoCptModel;
  using refactor::PKTwoCptODE;
  using refactor::PKODEModel;
  using refactor::PKOdeFunctorRateAdaptor;

  rate[0] = 1200;
  rate[1] = 2000;
  rate[2] = 3000;
  y0[0] = 800;
  y0[1] = 0;
  y0[2] = 8000;
  std::vector<stan::math::var> theta = to_var(par);

  PKTwoCptModel<double, double, double, var> model0(t0, y0, rate, theta);
  std::vector<double> yvec(y0.data(), y0.data() + y0.size());
  PKOdeFunctorRateAdaptor<PKTwoCptODE, double> f1(model0.f());
  using model_t = PKODEModel<double, double, double, var, PKTwoCptODE, int>;
  model_t model(t0, y0, rate, model0.par(), model0.f(), model0.ncmt());

  Eigen::Matrix<torsten::scalar_t<model_t>, Eigen::Dynamic, 1> y;
  std::vector<std::vector<torsten::scalar_t<model_t>> > y1;
  std::vector<double> g, g1;
  ts[0] = 20.0;
  ts.resize(1);

  auto test_it = [&]() {
    EXPECT_FLOAT_EQ(y(0).val(), y1[0][0].val());
    EXPECT_FLOAT_EQ(y(1).val(), y1[0][1].val());
    EXPECT_FLOAT_EQ(y(2).val(), y1[0][2].val());

    for (int i = 0; i < y0.size(); ++i) {
      stan::math::set_zero_all_adjoints();    
      y(i).grad(theta, g);
      stan::math::set_zero_all_adjoints();    
      y1[0][i].grad(theta, g1);
      for (size_t j = 0; j < theta.size(); ++j) {
        EXPECT_FLOAT_EQ(g[j], g1[j]);
      }
    }
  };

  PkOdeIntegrator<StanRk45> integ1(rtol, atol, max_num_steps, msgs);
  PkOdeIntegrator<StanAdams> integ2(rtol, atol, max_num_steps, msgs);
  PkOdeIntegrator<StanBdf> integ3(rtol, atol, max_num_steps, msgs);
  PkOdeIntegrator<PkAdams> integ4(rtol, atol, max_num_steps, msgs);
  PkOdeIntegrator<PkBdf> integ5(rtol, atol, max_num_steps, msgs);

  y = model.solve(ts[0], integ1);
  y1 = stan::math::integrate_ode_rk45(f1, yvec, t0, ts, model.par(), model.rate(), x_i, msgs); // NOLINT
  test_it();

  y = model.solve(ts[0], integ2);
  y1 = stan::math::integrate_ode_adams(f1, yvec, t0, ts, model.par(), model.rate(), x_i, msgs); // NOLINT
  test_it();

  y = model.solve(ts[0], integ3);
  y1 = stan::math::integrate_ode_bdf(f1, yvec, t0, ts, model.par(), model.rate(), x_i, msgs); // NOLINT
  test_it();

  y = model.solve(ts[0], integ4);
  y1 = torsten::dsolve::pk_integrate_ode_adams(f1, yvec, t0, ts, model.par(), model.rate(), x_i, msgs); // NOLINT
  test_it();

  y = model.solve(ts[0], integ5);
  y1 = torsten::dsolve::pk_integrate_ode_bdf(f1, yvec, t0, ts, model.par(), model.rate(), x_i, msgs); // NOLINT
  test_it();
}

TEST_F(TorstenCptOdeModelTest, general_ode_solver_par_rate_sens) {
  using stan::math::var;
  using stan::math::to_var;
  using refactor::PKTwoCptModel;
  using refactor::PKTwoCptODE;
  using refactor::PKODEModel;
  using refactor::PKOdeFunctorRateAdaptor;

  rate[0] = 1200;
  rate[1] = 2000;
  rate[2] = 3000;
  y0[0] = 800;
  y0[1] = 0;
  y0[2] = 8000;
  std::vector<stan::math::var> theta = to_var(par);
  std::vector<stan::math::var> rate_var = to_var(rate);
  // using model_t = PKTwoCptModel<double, double, var, var>;
  // model_t model(t0, y0, rate_var, theta);
  // std::vector<double> yvec(y0.data(), y0.data() + y0.size());
  // PKODERateAdaptor<model_t> adaptor(model);

  PKTwoCptModel<double, double, var, var> model0(t0, y0, rate_var, theta);
  std::vector<double> yvec(y0.data(), y0.data() + y0.size());
  PKOdeFunctorRateAdaptor<PKTwoCptODE, var> f1(model0.f(), theta.size());
  using model_t = PKODEModel<double, double, var, var, PKTwoCptODE, int>;
  model_t model(t0, y0, model0.rate(), model0.par(), model0.f(), model0.ncmt()); // NOLINT
  theta.insert(theta.end(), rate_var.begin(), rate_var.end());

  Eigen::Matrix<torsten::scalar_t<model_t>, Eigen::Dynamic, 1> y;
  std::vector<std::vector<torsten::scalar_t<model_t>> > y1;
  std::vector<double> g, g1;
  ts[0] = 20.0;
  ts.resize(1);

  auto test_it = [&]() {
    EXPECT_FLOAT_EQ(y(0).val(), y1[0][0].val());
    EXPECT_FLOAT_EQ(y(1).val(), y1[0][1].val());
    EXPECT_FLOAT_EQ(y(2).val(), y1[0][2].val());

    for (int i = 0; i < y0.size(); ++i) {
      stan::math::set_zero_all_adjoints();    
      y(i).grad(theta, g);
      stan::math::set_zero_all_adjoints();    
      y1[0][i].grad(theta, g1);
      for (size_t j = 0; j < theta.size(); ++j) {
        EXPECT_FLOAT_EQ(g[j], g1[j]);
      }
    }

    for (int i = 0; i < y0.size(); ++i) {
      stan::math::set_zero_all_adjoints();    
      y(i).grad(rate_var, g);
      stan::math::set_zero_all_adjoints();    
      y1[0][i].grad(rate_var, g1);
      for (size_t j = 0; j < theta.size(); ++j) {
        EXPECT_FLOAT_EQ(g[j], g1[j]);
      }
    }
  };

  PkOdeIntegrator<StanRk45> integ1(rtol, atol, max_num_steps, msgs);
  PkOdeIntegrator<StanAdams> integ2(rtol, atol, max_num_steps, msgs);
  PkOdeIntegrator<StanBdf> integ3(rtol, atol, max_num_steps, msgs);
  PkOdeIntegrator<PkAdams> integ4(rtol, atol, max_num_steps, msgs);
  PkOdeIntegrator<PkBdf> integ5(rtol, atol, max_num_steps, msgs);

  y = model.solve(ts[0], integ1);
  y1 = stan::math::integrate_ode_rk45(f1, yvec, t0, ts, theta, x_r, x_i, msgs); // NOLINT
  test_it();

  y = model.solve(ts[0], integ2);
  y1 = stan::math::integrate_ode_adams(f1, yvec, t0, ts, theta, x_r, x_i, msgs); // NOLINT
  test_it();

  y = model.solve(ts[0], integ3);
  y1 = stan::math::integrate_ode_bdf(f1, yvec, t0, ts, theta, x_r, x_i, msgs); // NOLINT
  test_it();

  y = model.solve(ts[0], integ4);
  y1 = torsten::dsolve::pk_integrate_ode_adams(f1, yvec, t0, ts, theta, x_r, x_i, msgs); // NOLINT
  test_it();

  y = model.solve(ts[0], integ5);
  y1 = torsten::dsolve::pk_integrate_ode_bdf(f1, yvec, t0, ts, theta, x_r, x_i, msgs); // NOLINT
  test_it();
}
