#include <stan/math.hpp>
#include <stan/math/torsten/pmx_ode_model.hpp>
#include <stan/math/rev/core.hpp>
#include <test/unit/math/rev/mat/fun/util.hpp>
#include <test/unit/math/torsten/pmx_cpt_model_test_fixture.hpp>

#include <test/unit/util.hpp>
#include <gtest/gtest.h>

TEST_F(TorstenTwoCptModelTest, pk_integrator_t0_var) {
  using stan::math::var;
  using stan::math::to_var;
  using refactor::PMXTwoCptModel;
  using refactor::PMXTwoCptODE;
  using refactor::PKODEModel;

  y0(0) = 100.0;
  y0(1) = 1000.0;
  y0(2) = 150.0;

  PMXTwoCptModel<double, double, double, double> model0(t0, y0, rate, CL, Q, V2, V3, ka); // NOLINT

  using model1_t = PKODEModel<double, double, double, double, PMXTwoCptODE>;
  using model2_t = PKODEModel<var, double, double, double, PMXTwoCptODE>;

  double t1 = 1.5;
  var t1_v = t1;

  {
    PMXOdeIntegrator<PkBdf> integ(rtol, atol, max_num_steps, msgs);
    auto f1 = [&] (std::vector<double>& x) {
      model1_t model(x[0], y0, rate, model0.par(), model0.f());
      return model.solve(t1, integ);
    };
    auto f2 = [&] (std::vector<stan::math::var>& x) {
      model2_t model(x[0], y0, rate, model0.par(), model0.f());
      return model.solve(t1_v, integ);
    };
    std::vector<double> t0v{0.0};
    torsten::test::test_grad(f1, f2, t0v, 2e-5, 1e-6, 1e-3, 1e-5);
  }

  {
    PMXOdeIntegrator<PkAdams> integ(rtol, atol, max_num_steps, msgs);
    auto f1 = [&] (std::vector<double>& x) {
      model1_t model(x[0], y0, rate, model0.par(), model0.f());
      return model.solve(t1, integ);
    };
    auto f2 = [&] (std::vector<stan::math::var>& x) {
      model2_t model(x[0], y0, rate, model0.par(), model0.f());
      return model.solve(t1_v, integ);
    };
    std::vector<double> t0v{0.0};
    torsten::test::test_grad(f1, f2, t0v, 2e-5, 1e-6, 1e-3, 1e-5);
  }
}

TEST_F(TorstenTwoCptModelTest, pk_integrator_rate_t0_var) {
  using stan::math::var;
  using stan::math::to_var;
  using refactor::PMXTwoCptModel;
  using refactor::PMXTwoCptODE;
  using refactor::PKODEModel;

  y0(0) = 100.0;
  y0(1) = 1000.0;
  y0(2) = 150.0;

  PMXTwoCptModel<double, double, double, double> model0(t0, y0, rate, CL, Q, V2, V3, ka); // NOLINT

  using model1_t = PKODEModel<double, double, double, double, PMXTwoCptODE>;
  using model2_t = PKODEModel<var, double, var, double, PMXTwoCptODE>;

  double t1 = 1.5;
  var t1_v = t1;

  {
    PMXOdeIntegrator<PkBdf> integ(rtol, atol, max_num_steps, msgs);
    auto f1 = [&] (std::vector<double>& x) {
      double t0 = x[0];
      std::vector<double> rate1(x.begin() + 1, x.end());
      model1_t model(t0, y0, rate1, model0.par(), model0.f());
      return model.solve(t1, integ);
    };
    auto f2 = [&] (std::vector<stan::math::var>& x) {
      var t0 = x[0];
      std::vector<var> rate1(x.begin() + 1, x.end());
      model2_t model(t0, y0, rate1, model0.par(), model0.f());
      return model.solve(t1_v, integ);
    };
    std::vector<double> param{t0, 100.0, 200.0, 0.0};
    torsten::test::test_grad(f1, f2, param, 2e-5, 1e-6, 1e-3, 1e-5);
  }

  {
    PMXOdeIntegrator<PkAdams> integ(rtol, atol, max_num_steps, msgs);
    auto f1 = [&] (std::vector<double>& x) {
      double t0 = x[0];
      std::vector<double> rate1(x.begin() + 1, x.end());
      model1_t model(t0, y0, rate1, model0.par(), model0.f());
      return model.solve(t1, integ);
    };
    auto f2 = [&] (std::vector<stan::math::var>& x) {
      var t0 = x[0];
      std::vector<var> rate1(x.begin() + 1, x.end());
      model2_t model(t0, y0, rate1, model0.par(), model0.f());
      return model.solve(t1_v, integ);
    };
    std::vector<double> param{t0, 100.0, 200.0, 0.0};
    torsten::test::test_grad(f1, f2, param, 2e-5, 1e-6, 1e-3, 1e-5);
  }
}

TEST_F(TorstenTwoCptModelTest, pk_integrator_ts_var) {
  using stan::math::var;
  using stan::math::to_var;
  using refactor::PMXTwoCptModel;
  using refactor::PMXTwoCptODE;
  using refactor::PKODEModel;

  y0(0) = 100.0;
  y0(1) = 1000.0;
  y0(2) = 0.0;

  refactor::PKRec<var> y0_v(stan::math::to_var(y0));

  PMXTwoCptModel<double, double, double, double> model0(t0, y0, rate, CL, Q, V2, V3, ka); // NOLINT

  using model1_t = PKODEModel<double, double, double, double, PMXTwoCptODE>;
  using model2_t = PKODEModel<var, var, double, double, PMXTwoCptODE>;
  using model3_t = PKODEModel<var, double, double, double, PMXTwoCptODE>;

  double t0 = 0.0;
  stan::math::var t0_v = 0.0;
  std::vector<double> dtv{1.5};

  {
    PMXOdeIntegrator<PkBdf> integ(rtol, atol, max_num_steps, msgs);
    auto f1 = [&] (std::vector<double>& x) {
      model1_t model(t0, y0, rate, model0.par(), model0.f());
      return model.solve(x[0], integ);
    };
    auto f2 = [&] (std::vector<stan::math::var>& x) {
      model2_t model(t0_v, y0_v, rate, model0.par(), model0.f());
      return model.solve(x[0], integ);
    };
    torsten::test::test_grad(f1, f2, dtv, 2e-5, 1e-6, 1e-5, 1e-5);
  }

  {
    PMXOdeIntegrator<PkBdf> integ(rtol, atol, max_num_steps, msgs);
    auto f1 = [&] (std::vector<double>& x) {
      model1_t model(t0, y0, rate, model0.par(), model0.f());
      return model.solve(x[0], integ);
    };
    auto f2 = [&] (std::vector<stan::math::var>& x) {
      model3_t model(t0_v, y0, rate, model0.par(), model0.f());
      return model.solve(x[0], integ);
    };
    torsten::test::test_grad(f1, f2, dtv, 2e-5, 1e-6, 1e-5, 1e-5);
  }

  {
    PMXOdeIntegrator<PkAdams> integ(rtol, atol, max_num_steps, msgs);
    auto f1 = [&] (std::vector<double>& x) {
      model1_t model(t0, y0, rate, model0.par(), model0.f());
      return model.solve(x[0], integ);
    };
    auto f2 = [&] (std::vector<stan::math::var>& x) {
      model2_t model(t0_v, y0_v, rate, model0.par(), model0.f());
      return model.solve(x[0], integ);
    };
    torsten::test::test_grad(f1, f2, dtv, 2e-5, 1e-6, 1e-4, 1e-5);
  }

  {
    PMXOdeIntegrator<PkAdams> integ(rtol, atol, max_num_steps, msgs);
    auto f1 = [&] (std::vector<double>& x) {
      model1_t model(t0, y0, rate, model0.par(), model0.f());
      return model.solve(x[0], integ);
    };
    auto f2 = [&] (std::vector<stan::math::var>& x) {
      model3_t model(t0_v, y0, rate, model0.par(), model0.f());
      return model.solve(x[0], integ);
    };
    torsten::test::test_grad(f1, f2, dtv, 2e-5, 1e-6, 1e-4, 1e-5);
  }
}

TEST_F(TorstenTwoCptModelTest, pk_integrator_y0_var) {
  using stan::math::var;
  using stan::math::to_var;
  using refactor::PMXTwoCptModel;
  using refactor::PMXTwoCptODE;
  using refactor::PKODEModel;

  y0(0) = 100.0;
  y0(1) = 1000.0;
  y0(2) = 150.0;

  std::vector<double> y0_vec(stan::math::to_array_1d(y0));

  PMXTwoCptModel<double, double, double, double> model0(t0, y0, rate, CL, Q, V2, V3, ka); // NOLINT

  using model1_t = PKODEModel<double, double, double, double, PMXTwoCptODE>;
  using model2_t = PKODEModel<double, var, double, double, PMXTwoCptODE>;
  using model3_t = PKODEModel<var, var, double, double, PMXTwoCptODE>;

  double t0 = 0.2;
  stan::math::var t0_v = t0;
  double t1 = 1.0;

  {
    PMXOdeIntegrator<PkBdf> integ(rtol, atol, max_num_steps, msgs);
    auto f1 = [&] (std::vector<double>& x) {
      refactor::PKRec<double> xr = stan::math::to_vector(x);
      model1_t model(t0, xr, rate, model0.par(), model0.f());
      return model.solve(t1, integ);
    };
    auto f2 = [&] (std::vector<var>& x) {
      refactor::PKRec<var> xr = stan::math::to_vector(x);
      model2_t model(t0, xr, rate, model0.par(), model0.f());
      return model.solve(t1, integ);
    };
    torsten::test::test_grad(f1, f2, y0_vec, 2e-5, 1e-6, 1e-3, 1e-3);
  }

  {
    PMXOdeIntegrator<PkBdf> integ(rtol, atol, max_num_steps, msgs);
    auto f1 = [&] (std::vector<double>& x) {
      refactor::PKRec<double> xr = stan::math::to_vector(x);
      model1_t model(t0, xr, rate, model0.par(), model0.f());
      return model.solve(t1, integ);
    };
    auto f2 = [&] (std::vector<var>& x) {
      refactor::PKRec<var> xr = stan::math::to_vector(x);
      model3_t model(t0_v, xr, rate, model0.par(), model0.f());
      return model.solve(t1, integ);
    };
    torsten::test::test_grad(f1, f2, y0_vec, 2e-5, 1e-6, 1e-3, 1e-3);
  }

  {
    PMXOdeIntegrator<PkAdams> integ(rtol, atol, max_num_steps, msgs);
    auto f1 = [&] (std::vector<double>& x) {
      refactor::PKRec<double> xr = stan::math::to_vector(x);
      model1_t model(t0, xr, rate, model0.par(), model0.f());
      return model.solve(t1, integ);
    };
    auto f2 = [&] (std::vector<var>& x) {
      refactor::PKRec<var> xr = stan::math::to_vector(x);
      model2_t model(t0, xr, rate, model0.par(), model0.f());
      return model.solve(t1, integ);
    };
    torsten::test::test_grad(f1, f2, y0_vec, 2e-5, 1e-6, 1e-3, 1e-3);
  }

  {
    PMXOdeIntegrator<PkAdams> integ(rtol, atol, max_num_steps, msgs);
    auto f1 = [&] (std::vector<double>& x) {
      refactor::PKRec<double> xr = stan::math::to_vector(x);
      model1_t model(t0, xr, rate, model0.par(), model0.f());
      return model.solve(t1, integ);
    };
    auto f2 = [&] (std::vector<var>& x) {
      refactor::PKRec<var> xr = stan::math::to_vector(x);
      model3_t model(t0_v, xr, rate, model0.par(), model0.f());
      return model.solve(t1, integ);
    };
    torsten::test::test_grad(f1, f2, y0_vec, 2e-5, 1e-6, 1e-3, 1e-3);
  }
}

TEST_F(TorstenTwoCptModelTest, pk_bdf_integrator_dt_var) {
  using stan::math::var;
  using stan::math::to_var;
  using refactor::PMXTwoCptModel;
  using refactor::PMXTwoCptODE;
  using refactor::PKODEModel;

  y0(0) = 100.0;
  y0(1) = 1000.0;
  y0(2) = 0.0;

  PMXOdeIntegrator<PkBdf> integ(rtol, atol, max_num_steps, msgs);
  PMXTwoCptModel<double, double, double, double> model0(t0, y0, rate, CL, Q, V2, V3, ka); // NOLINT

  PKODEModel<double, double, double, double, PMXTwoCptODE> model1(t0, y0, rate, model0.par(), model0.f());
  using model_t = PKODEModel<var, double, double, double, PMXTwoCptODE>;
  var t0_v = t0;
  model_t model2(t0_v, y0, rate, model0.par(), model0.f());

  std::vector<double> dtv{ts[0]};
  {
    auto f1 = [&] (std::vector<double>& x) { return model1.solve(x[0], integ); };
    auto f2 = [&] (std::vector<var>& x) { return model2.solve(x[0], integ); };
    torsten::test::test_grad(f1, f2, dtv, 2e-5, 1e-6, 1e-3, 1e-3);
  }
}

TEST_F(TorstenTwoCptModelTest, pk_adams_integrator_dt_var) {
  using stan::math::var;
  using stan::math::to_var;
  using refactor::PMXTwoCptModel;
  using refactor::PMXTwoCptODE;
  using refactor::PKODEModel;

  y0(0) = 100.0;
  y0(1) = 1000.0;
  y0(2) = 0.0;

  PMXOdeIntegrator<PkAdams> integ(rtol, atol, max_num_steps, msgs);
  PMXTwoCptModel<double, double, double, double> model0(t0, y0, rate, CL, Q, V2, V3, ka); // NOLINT

  PKODEModel<double, double, double, double, PMXTwoCptODE> model1(t0, y0, rate, model0.par(), model0.f());
  using model_t = PKODEModel<var, double, double, double, PMXTwoCptODE>;
  var t0_v = t0;
  model_t model2(t0_v, y0, rate, model0.par(), model0.f());

  std::vector<double> dtv{ts[0]};
  {
    auto f1 = [&] (std::vector<double>& x) { return model1.solve(x[0], integ); };
    auto f2 = [&] (std::vector<var>& x) { return model2.solve(x[0], integ); };
    torsten::test::test_grad(f1, f2, dtv, 2e-5, 1e-6, 1e-3, 1e-4);
  }
}

TEST_F(TorstenTwoCptModelTest, general_ode_solver) {
  using stan::math::var;
  using stan::math::to_var;
  using refactor::PMXTwoCptModel;
  using refactor::PMXTwoCptODE;
  using refactor::PKODEModel;
  using refactor::PMXOdeFunctorRateAdaptor;

  rate[0] = 1200;
  rate[1] = 2000;
  rate[2] = 3000;
  PMXTwoCptModel<double, double, double, double> model0(t0, y0, rate, CL, Q, V2, V3, ka); // NOLINT
  std::vector<double> yvec(y0.data(), y0.data() + y0.size());
  PMXOdeFunctorRateAdaptor<PMXTwoCptODE, double> f1(model0.f());
  using model_t = PKODEModel<double, double, double, double, PMXTwoCptODE>;
  model_t model(t0, y0, rate, model0.par(), model0.f());

  Eigen::Matrix<torsten::scalar_t<model_t>, Eigen::Dynamic, 1> y;
  std::vector<std::vector<double> > y1;
  ts[0] = 20.0;
  ts.resize(1);

  PMXOdeIntegrator<StanRk45> integ1(rtol, atol, max_num_steps, msgs);
  y = model.solve(ts[0], integ1);
  y1 = stan::math::integrate_ode_rk45(f1, yvec, t0, ts, model.par(), model.rate(), x_i, msgs); // NOLINT
  EXPECT_FLOAT_EQ(y(0), y1[0][0]);
  EXPECT_FLOAT_EQ(y(1), y1[0][1]);
  EXPECT_FLOAT_EQ(y(2), y1[0][2]);

  PMXOdeIntegrator<StanAdams> integ2(rtol, atol, max_num_steps, msgs);
  y = model.solve(ts[0], integ2);
  y1 = stan::math::integrate_ode_adams(f1, yvec, t0, ts, model.par(), model.rate(), x_i, msgs); // NOLINT
  EXPECT_FLOAT_EQ(y(0), y1[0][0]);
  EXPECT_FLOAT_EQ(y(1), y1[0][1]);
  EXPECT_FLOAT_EQ(y(2), y1[0][2]);

  PMXOdeIntegrator<StanBdf> integ3(rtol, atol, max_num_steps, msgs);
  y = model.solve(ts[0], integ3);
  y1 = stan::math::integrate_ode_bdf(f1, yvec, t0, ts, model.par(), model.rate(), x_i, msgs); // NOLINT
  EXPECT_FLOAT_EQ(y(0), y1[0][0]);
  EXPECT_FLOAT_EQ(y(1), y1[0][1]);
  EXPECT_FLOAT_EQ(y(2), y1[0][2]);

  PMXOdeIntegrator<PkAdams> integ4(rtol, atol, max_num_steps, msgs);
  y = model.solve(ts[0], integ4);
  y1 = torsten::pmx_integrate_ode_adams(f1, yvec, t0, ts, model.par(), model.rate(), x_i, msgs); // NOLINT
  EXPECT_FLOAT_EQ(y(0), y1[0][0]);
  EXPECT_FLOAT_EQ(y(1), y1[0][1]);
  EXPECT_FLOAT_EQ(y(2), y1[0][2]);

  PMXOdeIntegrator<PkBdf> integ5(rtol, atol, max_num_steps, msgs);
  y = model.solve(ts[0], integ5);
  y1 = torsten::pmx_integrate_ode_bdf(f1, yvec, t0, ts, model.par(), model.rate(), x_i, msgs); // NOLINT
  EXPECT_FLOAT_EQ(y(0), y1[0][0]);
  EXPECT_FLOAT_EQ(y(1), y1[0][1]);
  EXPECT_FLOAT_EQ(y(2), y1[0][2]);
}

TEST_F(TorstenTwoCptModelTest, general_ode_solver_y0) {
  using stan::math::var;
  using stan::math::to_var;
  using refactor::PMXTwoCptModel;
  using refactor::PMXTwoCptODE;
  using refactor::PKODEModel;
  using refactor::PMXOdeFunctorRateAdaptor;

  rate[0] = 1200;
  rate[1] = 2000;
  rate[2] = 3000;
  y0[0] = 800;
  y0[1] = 0;
  y0[2] = 8000;
  PMXTwoCptModel<double, double, double, double> model0(t0, y0, rate, CL, Q, V2, V3, ka); // NOLINT
  std::vector<double> yvec(y0.data(), y0.data() + y0.size());
  PMXOdeFunctorRateAdaptor<PMXTwoCptODE, double> f1(model0.f());
  using model_t = PKODEModel<double, double, double, double, PMXTwoCptODE>;
  model_t model(t0, y0, rate, model0.par(), model0.f());

  Eigen::Matrix<torsten::scalar_t<model_t>, Eigen::Dynamic, 1> y;
  std::vector<std::vector<double> > y1;
  ts[0] = 20.0;
  ts.resize(1);

  PMXOdeIntegrator<StanRk45> integ1(rtol, atol, max_num_steps, msgs);
  y = model.solve(ts[0], integ1);
  y1 = stan::math::integrate_ode_rk45(f1, yvec, t0, ts, model.par(), model.rate(), x_i, msgs); // NOLINT
  EXPECT_FLOAT_EQ(y(0), y1[0][0]);
  EXPECT_FLOAT_EQ(y(1), y1[0][1]);
  EXPECT_FLOAT_EQ(y(2), y1[0][2]);

  PMXOdeIntegrator<StanAdams> integ2(rtol, atol, max_num_steps, msgs);
  y = model.solve(ts[0], integ2);
  y1 = stan::math::integrate_ode_adams(f1, yvec, t0, ts, model.par(), model.rate(), x_i, msgs); // NOLINT
  EXPECT_FLOAT_EQ(y(0), y1[0][0]);
  EXPECT_FLOAT_EQ(y(1), y1[0][1]);
  EXPECT_FLOAT_EQ(y(2), y1[0][2]);

  PMXOdeIntegrator<StanBdf> integ3(rtol, atol, max_num_steps, msgs);
  y = model.solve(ts[0], integ3);
  y1 = stan::math::integrate_ode_bdf(f1, yvec, t0, ts, model.par(), model.rate(), x_i, msgs); // NOLINT
  EXPECT_FLOAT_EQ(y(0), y1[0][0]);
  EXPECT_FLOAT_EQ(y(1), y1[0][1]);
  EXPECT_FLOAT_EQ(y(2), y1[0][2]);

  PMXOdeIntegrator<PkAdams> integ4(rtol, atol, max_num_steps, msgs);
  y = model.solve(ts[0], integ4);
  y1 = torsten::pmx_integrate_ode_adams(f1, yvec, t0, ts, model.par(), model.rate(), x_i, msgs); // NOLINT
  EXPECT_FLOAT_EQ(y(0), y1[0][0]);
  EXPECT_FLOAT_EQ(y(1), y1[0][1]);
  EXPECT_FLOAT_EQ(y(2), y1[0][2]);

  PMXOdeIntegrator<PkBdf> integ5(rtol, atol, max_num_steps, msgs);
  y = model.solve(ts[0], integ5);
  y1 = torsten::pmx_integrate_ode_bdf(f1, yvec, t0, ts, model.par(), model.rate(), x_i, msgs); // NOLINT
  EXPECT_FLOAT_EQ(y(0), y1[0][0]);
  EXPECT_FLOAT_EQ(y(1), y1[0][1]);
  EXPECT_FLOAT_EQ(y(2), y1[0][2]);
}

TEST_F(TorstenTwoCptModelTest, general_ode_solver_par_sens) {
  using stan::math::var;
  using stan::math::to_var;
  using refactor::PMXTwoCptModel;
  using refactor::PMXTwoCptODE;
  using refactor::PKODEModel;
  using refactor::PMXOdeFunctorRateAdaptor;

  rate[0] = 1200;
  rate[1] = 2000;
  rate[2] = 3000;
  y0[0] = 800;
  y0[1] = 0;
  y0[2] = 8000;
  std::vector<stan::math::var> theta = to_var(par);

  PMXTwoCptModel<double, double, double, var> model0(t0, y0, rate, theta);
  std::vector<double> yvec(y0.data(), y0.data() + y0.size());
  PMXOdeFunctorRateAdaptor<PMXTwoCptODE, double> f1(model0.f());
  using model_t = PKODEModel<double, double, double, var, PMXTwoCptODE>;
  model_t model(t0, y0, rate, model0.par(), model0.f());

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

  PMXOdeIntegrator<StanRk45> integ1(rtol, atol, max_num_steps, msgs);
  PMXOdeIntegrator<StanAdams> integ2(rtol, atol, max_num_steps, msgs);
  PMXOdeIntegrator<StanBdf> integ3(rtol, atol, max_num_steps, msgs);
  PMXOdeIntegrator<PkAdams> integ4(rtol, atol, max_num_steps, msgs);
  PMXOdeIntegrator<PkBdf> integ5(rtol, atol, max_num_steps, msgs);

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
  y1 = torsten::pmx_integrate_ode_adams(f1, yvec, t0, ts, model.par(), model.rate(), x_i, msgs); // NOLINT
  test_it();

  y = model.solve(ts[0], integ5);
  y1 = torsten::pmx_integrate_ode_bdf(f1, yvec, t0, ts, model.par(), model.rate(), x_i, msgs); // NOLINT
  test_it();
}

TEST_F(TorstenTwoCptModelTest, general_ode_solver_par_rate_sens) {
  using stan::math::var;
  using stan::math::to_var;
  using refactor::PMXTwoCptModel;
  using refactor::PMXTwoCptODE;
  using refactor::PKODEModel;
  using refactor::PMXOdeFunctorRateAdaptor;

  rate[0] = 1200;
  rate[1] = 2000;
  rate[2] = 3000;
  y0[0] = 800;
  y0[1] = 0;
  y0[2] = 8000;
  std::vector<stan::math::var> theta = to_var(par);
  std::vector<stan::math::var> rate_var = to_var(rate);
  // using model_t = PMXTwoCptModel<double, double, var, var>;
  // model_t model(t0, y0, rate_var, theta);
  // std::vector<double> yvec(y0.data(), y0.data() + y0.size());
  // PKODERateAdaptor<model_t> adaptor(model);

  PMXTwoCptModel<double, double, var, var> model0(t0, y0, rate_var, theta);
  std::vector<double> yvec(y0.data(), y0.data() + y0.size());
  PMXOdeFunctorRateAdaptor<PMXTwoCptODE, var> f1(model0.f(), theta.size());
  using model_t = PKODEModel<double, double, var, var, PMXTwoCptODE>;
  model_t model(t0, y0, model0.rate(), model0.par(), model0.f()); // NOLINT
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

  PMXOdeIntegrator<StanRk45> integ1(rtol, atol, max_num_steps, msgs);
  PMXOdeIntegrator<StanAdams> integ2(rtol, atol, max_num_steps, msgs);
  PMXOdeIntegrator<StanBdf> integ3(rtol, atol, max_num_steps, msgs);
  PMXOdeIntegrator<PkAdams> integ4(rtol, atol, max_num_steps, msgs);
  PMXOdeIntegrator<PkBdf> integ5(rtol, atol, max_num_steps, msgs);

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
  y1 = torsten::pmx_integrate_ode_adams(f1, yvec, t0, ts, theta, x_r, x_i, msgs); // NOLINT
  test_it();

  y = model.solve(ts[0], integ5);
  y1 = torsten::pmx_integrate_ode_bdf(f1, yvec, t0, ts, theta, x_r, x_i, msgs); // NOLINT
  test_it();
}
