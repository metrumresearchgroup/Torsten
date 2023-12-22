#include <stan/math.hpp>
#include <stan/math/torsten/pmx_ode_model.hpp>
#include <stan/math/torsten/dsolve.hpp>
#include <stan/math/rev/core.hpp>
#include <test/unit/math/rev/fun/util.hpp>
#include <stan/math/torsten/test/unit/pmx_cpt_model_test_fixture.hpp>

#include <test/unit/util.hpp>
#include <gtest/gtest.h>

using stan::math::var;
using stan::math::to_var;
using torsten::PKRec;
using torsten::PMXTwoCptModel;
using torsten::PMXTwoCptODE;
using torsten::PKODEModel;
using torsten::PMXOdeFunctorRateAdaptor;
using torsten::dsolve::PMXCvodesIntegrator;
using torsten::dsolve::PMXOdeintIntegrator;
using torsten::dsolve::PMXVariadicOdeSystem;

PMXOdeIntegrator<PMXVariadicOdeSystem,  PMXCvodesIntegrator<CV_ADAMS, CV_STAGGERED>> integrator_adams;
PMXOdeIntegrator<PMXVariadicOdeSystem, PMXCvodesIntegrator<CV_BDF, CV_STAGGERED>> integrator_bdf;

using scheme_t = boost::numeric::odeint::runge_kutta_dopri5<std::vector<double>, double, std::vector<double>, double>;
PMXOdeIntegrator<PMXVariadicOdeSystem, PMXOdeintIntegrator<scheme_t>> integrator_rk45;

PMXTwoCptODE f0;

TEST_F(TorstenTwoCptModelTest, pk_integrator_t0_var) {
  y0(0) = 100.0;
  y0(1) = 1000.0;
  y0(2) = 150.0;

  PMXTwoCptModel<double> model0(CL, Q, V2, V3, ka); // NOLINT

  using model1_t = PKODEModel<double, PMXTwoCptODE>;
  using model2_t = PKODEModel<double, PMXTwoCptODE>;

  double t1 = 1.5;
  var t1_v = t1;

  {
    PMXOdeIntegrator<PMXVariadicOdeSystem, PMXCvodesIntegrator<CV_BDF, CV_STAGGERED>> integ(rtol, atol, max_num_steps, msgs);
    auto f1 = [&] (const std::vector<double>& x) {
      model1_t model(model0.par(), y0.size(), model0.f());
      PKRec<double> y(y0);
      model.solve(y, x[0], t1, rate, integ); return y;
    };
    auto f2 = [&] (const std::vector<stan::math::var>& x) {
      model2_t model(model0.par(), y0.size(), model0.f());
      PKRec<var> y(to_var(y0));
      model.solve(y, x[0], t1_v, rate, integ); return y;
    };
    std::vector<double> t0v{0.0};
    torsten::test::test_grad(f1, f2, t0v, 2e-5, 1e-6, 1e-3, 1e-5);
  }

  {
    torsten::dsolve::PMXOdeIntegrator<PMXVariadicOdeSystem,  PMXCvodesIntegrator<CV_ADAMS, CV_STAGGERED>> integ(rtol, atol, max_num_steps, msgs);
    auto f1 = [&] (const std::vector<double>& x) {
      model1_t model(model0.par(), y0.size(), model0.f());
      PKRec<double> y(y0);
      model.solve(y, x[0], t1, rate, integ); return y;
    };
    auto f2 = [&] (const std::vector<stan::math::var>& x) {
      model2_t model(model0.par(), y0.size(), model0.f());
      PKRec<var> y(to_var(y0));
      model.solve(y, x[0], t1_v, rate, integ); return y;
    };
    std::vector<double> t0v{0.0};
    torsten::test::test_grad(f1, f2, t0v, 2e-5, 1e-6, 1e-3, 1e-5);
  }
}

TEST_F(TorstenTwoCptModelTest, pk_integrator_rate_t0_var) {
  y0(0) = 100.0;
  y0(1) = 1000.0;
  y0(2) = 150.0;

  PMXTwoCptModel<double> model0(CL, Q, V2, V3, ka); // NOLINT

  using model1_t = PKODEModel<double, PMXTwoCptODE>;
  using model2_t = PKODEModel<double, PMXTwoCptODE>;

  double t1 = 1.5;
  var t1_v = t1;

  {
    PMXOdeIntegrator<PMXVariadicOdeSystem, PMXCvodesIntegrator<CV_BDF, CV_STAGGERED>> integ(rtol, atol, max_num_steps, msgs);
    auto f1 = [&] (const std::vector<double>& x) {
      double t0 = x[0];
      std::vector<double> rate1(x.begin() + 1, x.end());
      model1_t model(model0.par(), y0.size(), model0.f());
      PKRec<double> y(y0);
      model.solve(y, t0, t1, rate1, integ); return y;
    };
    auto f2 = [&] (const std::vector<stan::math::var>& x) {
      var t0 = x[0];
      std::vector<var> rate1(x.begin() + 1, x.end());
      model2_t model(model0.par(), y0.size(), model0.f());
      PKRec<var> y(to_var(y0));
      model.solve(y, t0, t1_v, rate1, integ); return y;
    };
    std::vector<double> param{t0, 100.0, 200.0, 0.0};
    torsten::test::test_grad(f1, f2, param, 2e-5, 1e-6, 1e-3, 1e-5);
  }

  {
    torsten::dsolve::PMXOdeIntegrator<PMXVariadicOdeSystem,  PMXCvodesIntegrator<CV_ADAMS, CV_STAGGERED>> integ(rtol, atol, max_num_steps, msgs);
    auto f1 = [&] (const std::vector<double>& x) {
      double t0 = x[0];
      std::vector<double> rate1(x.begin() + 1, x.end());
      model1_t model(model0.par(), y0.size(), model0.f());
      PKRec<double> y(y0);
      model.solve(y, t0, t1, rate1, integ); return y;
    };
    auto f2 = [&] (const std::vector<stan::math::var>& x) {
      var t0 = x[0];
      std::vector<var> rate1(x.begin() + 1, x.end());
      model2_t model(model0.par(), y0.size(), model0.f());
      PKRec<var> y(to_var(y0));
      model.solve(y, t0, t1_v, rate1, integ); return y;
    };
    std::vector<double> param{t0, 100.0, 200.0, 0.0};
    torsten::test::test_grad(f1, f2, param, 2e-5, 1e-6, 1e-3, 1e-5);
  }
}

TEST_F(TorstenTwoCptModelTest, pk_integrator_ts_var) {
  y0(0) = 100.0;
  y0(1) = 1000.0;
  y0(2) = 0.0;

  torsten::PKRec<var> y0_v(stan::math::to_var(y0));

  PMXTwoCptModel<double> model0(CL, Q, V2, V3, ka); // NOLINT

  using model1_t = PKODEModel<double, PMXTwoCptODE>;
  using model2_t = PKODEModel<double, PMXTwoCptODE>;
  using model3_t = PKODEModel<double, PMXTwoCptODE>;

  double t0 = 0.0;
  stan::math::var t0_v = 0.0;
  std::vector<double> dtv{1.5};

  {
    PMXOdeIntegrator<PMXVariadicOdeSystem, PMXCvodesIntegrator<CV_BDF, CV_STAGGERED>> integ(rtol, atol, max_num_steps, msgs);
    auto f1 = [&] (const std::vector<double>& x) {
      model1_t model(model0.par(), y0.size(), model0.f());
      PKRec<double> y(y0);
      model.solve(y, t0, x[0], rate, integ); return y;
    };
    auto f2 = [&] (const std::vector<stan::math::var>& x) {
      model2_t model(model0.par(), y0.size(), model0.f());
      PKRec<var> y(y0_v);
      model.solve(y, t0, x[0], rate, integ); return y;
    };
    torsten::test::test_grad(f1, f2, dtv, 2e-5, 1e-6, 1e-5, 1e-5);
  }

  {
    PMXOdeIntegrator<PMXVariadicOdeSystem, PMXCvodesIntegrator<CV_BDF, CV_STAGGERED>> integ(rtol, atol, max_num_steps, msgs);
    auto f1 = [&] (const std::vector<double>& x) {
      model1_t model(model0.par(), y0.size(), model0.f());
      PKRec<double> y(y0);
      model.solve(y, t0, x[0], rate, integ); return y;
    };
    auto f2 = [&] (const std::vector<stan::math::var>& x) {
      model3_t model(model0.par(), y0.size(), model0.f());
      PKRec<var> y(to_var(y0));
      model.solve(y, t0, x[0], rate, integ); return y;
    };
    torsten::test::test_grad(f1, f2, dtv, 2e-5, 1e-6, 1e-5, 1e-5);
  }

  {
    torsten::dsolve::PMXOdeIntegrator<PMXVariadicOdeSystem,  PMXCvodesIntegrator<CV_ADAMS, CV_STAGGERED>> integ(rtol, atol, max_num_steps, msgs);
    auto f1 = [&] (const std::vector<double>& x) {
      model1_t model(model0.par(), y0.size(), model0.f());
      PKRec<double> y(y0);
      model.solve(y, t0, x[0], rate, integ); return y;
    };
    auto f2 = [&] (const std::vector<stan::math::var>& x) {
      model2_t model(model0.par(), y0.size(), model0.f());
      PKRec<var> y(y0_v);
      model.solve(y, t0, x[0], rate, integ); return y;
    };
    torsten::test::test_grad(f1, f2, dtv, 2e-5, 1e-6, 1e-4, 1e-5);
  }

  {
    torsten::dsolve::PMXOdeIntegrator<PMXVariadicOdeSystem,  PMXCvodesIntegrator<CV_ADAMS, CV_STAGGERED>> integ(rtol, atol, max_num_steps, msgs);
    auto f1 = [&] (const std::vector<double>& x) {
      model1_t model(model0.par(), y0.size(), model0.f());
      PKRec<double> y(y0);
      model.solve(y, t0, x[0], rate, integ); return y;
    };
    auto f2 = [&] (const std::vector<stan::math::var>& x) {
      model3_t model(model0.par(), y0.size(), model0.f());
      PKRec<var> y(to_var(y0));
      model.solve(y, t0, x[0], rate, integ); return y;
    };
    torsten::test::test_grad(f1, f2, dtv, 2e-5, 1e-6, 1e-4, 1e-5);
  }
}

TEST_F(TorstenTwoCptModelTest, pk_integrator_y0_var) {
  y0(0) = 100.0;
  y0(1) = 1000.0;
  y0(2) = 150.0;

  std::vector<double> y0_vec(stan::math::to_array_1d(y0));

  PMXTwoCptModel<double> model0(CL, Q, V2, V3, ka); // NOLINT

  using model1_t = PKODEModel<double, PMXTwoCptODE>;
  using model2_t = PKODEModel<double, PMXTwoCptODE>;
  using model3_t = PKODEModel<double, PMXTwoCptODE>;

  double t0 = 0.2;
  stan::math::var t0_v = t0;
  double t1 = 1.0;

  {
    PMXOdeIntegrator<PMXVariadicOdeSystem, PMXCvodesIntegrator<CV_BDF, CV_STAGGERED>> integ(rtol, atol, max_num_steps, msgs);
    auto f1 = [&] (const std::vector<double>& x) {
      torsten::PKRec<double> y = stan::math::to_vector(x);
      model1_t model(model0.par(), y0.size(), model0.f());
      model.solve(y, t0, t1, rate, integ); return y;
    };
    auto f2 = [&] (const std::vector<var>& x) {
      torsten::PKRec<var> y = stan::math::to_vector(x);
      model2_t model(model0.par(), y0.size(), model0.f());
      model.solve(y, t0, t1, rate, integ); return y;
    };
    torsten::test::test_grad(f1, f2, y0_vec, 2e-5, 1e-6, 1e-3, 1e-3);
  }

  {
    PMXOdeIntegrator<PMXVariadicOdeSystem, PMXCvodesIntegrator<CV_BDF, CV_STAGGERED>> integ(rtol, atol, max_num_steps, msgs);
    auto f1 = [&] (const std::vector<double>& x) {
      torsten::PKRec<double> y = stan::math::to_vector(x);
      model1_t model(model0.par(), y0.size(), model0.f());
      model.solve(y, t0, t1, rate, integ); return y;
    };
    auto f2 = [&] (const std::vector<var>& x) {
      torsten::PKRec<var> y = stan::math::to_vector(x);
      model3_t model(model0.par(), y0.size(), model0.f());
      model.solve(y, t0, t1, rate, integ); return y;
    };
    torsten::test::test_grad(f1, f2, y0_vec, 2e-5, 1e-6, 1e-3, 1e-3);
  }

  {
    torsten::dsolve::PMXOdeIntegrator<PMXVariadicOdeSystem,  PMXCvodesIntegrator<CV_ADAMS, CV_STAGGERED>> integ(rtol, atol, max_num_steps, msgs);
    auto f1 = [&] (const std::vector<double>& x) {
      torsten::PKRec<double> y = stan::math::to_vector(x);
      model1_t model(model0.par(), y0.size(), model0.f());
      model.solve(y, t0, t1, rate, integ); return y;
    };
    auto f2 = [&] (const std::vector<var>& x) {
      torsten::PKRec<var> y = stan::math::to_vector(x);
      model2_t model(model0.par(), y0.size(), model0.f());
      model.solve(y, t0, t1, rate, integ); return y;
    };
    torsten::test::test_grad(f1, f2, y0_vec, 2e-5, 1e-6, 1e-3, 1e-3);
  }

  {
    torsten::dsolve::PMXOdeIntegrator<PMXVariadicOdeSystem,  PMXCvodesIntegrator<CV_ADAMS, CV_STAGGERED>> integ(rtol, atol, max_num_steps, msgs);
    auto f1 = [&] (const std::vector<double>& x) {
      torsten::PKRec<double> y = stan::math::to_vector(x);
      model1_t model(model0.par(), y0.size(), model0.f());
      model.solve(y, t0, t1, rate, integ); return y;
    };
    auto f2 = [&] (const std::vector<var>& x) {
      torsten::PKRec<var> y = stan::math::to_vector(x);
      model3_t model(model0.par(), y0.size(), model0.f());
      model.solve(y, t0, t1, rate, integ); return y;
    };
    torsten::test::test_grad(f1, f2, y0_vec, 2e-5, 1e-6, 1e-3, 1e-3);
  }
}

TEST_F(TorstenTwoCptModelTest, pk_bdf_integrator_dt_var) {
  y0(0) = 100.0;
  y0(1) = 1000.0;
  y0(2) = 0.0;

  PMXTwoCptModel<double> model0(CL, Q, V2, V3, ka); // NOLINT

  PKODEModel<double, PMXTwoCptODE> model1(model0.par(), y0.size(), model0.f());
  using model_t = PKODEModel<double, PMXTwoCptODE>;
  var t0_v = t0;
  model_t model2(model0.par(), y0.size(), model0.f());

  std::vector<double> dtv{ts[0]};

  {
    PMXOdeIntegrator<PMXVariadicOdeSystem, PMXCvodesIntegrator<CV_BDF, CV_STAGGERED>> integ(rtol, atol, max_num_steps, msgs);
    auto f1 = [&] (const std::vector<double>& x) {PKRec<double> y(y0); model1.solve(y, t0, x[0], rate,integ); return y; };
    auto f2 = [&] (const std::vector<var>& x) {PKRec<var> y(to_var(y0)); model2.solve(y, t0_v, x[0], rate,integ); return y; };
    torsten::test::test_grad(f1, f2, dtv, 2e-5, 1e-6, 1e-3, 1e-3);
  }

  {
    torsten::dsolve::PMXOdeIntegrator<PMXVariadicOdeSystem,  PMXCvodesIntegrator<CV_ADAMS, CV_STAGGERED>> integ(rtol, atol, max_num_steps, msgs);
    auto f1 = [&] (const std::vector<double>& x) {PKRec<double> y(y0); model1.solve(y, t0, x[0], rate,integ); return y; };
    auto f2 = [&] (const std::vector<var>& x) {PKRec<var> y(to_var(y0)); model2.solve(y, t0_v, x[0], rate,integ); return y; };
    torsten::test::test_grad(f1, f2, dtv, 2e-5, 1e-6, 1e-3, 1e-4);
  }
}

TEST_F(TorstenTwoCptModelTest, general_ode_solver) {
  rate[0] = 1200;
  rate[1] = 2000;
  rate[2] = 3000;
  PMXTwoCptModel<double> model0(CL, Q, V2, V3, ka); // NOLINT
  PMXOdeFunctorRateAdaptor<PMXTwoCptODE> f1(f0);
  using model_t = PKODEModel<double, PMXTwoCptODE>;
  model_t model(model0.par(), y0.size(), model0.f());

  Eigen::Matrix<double, -1, 1> y;
  std::vector<Eigen::VectorXd> y1;
  ts[0] = 20.0;
  ts.resize(1);

  // torsten::dsolve::PMXOdeIntegrator<torsten::StanRk45> integ1(rtol, atol, max_num_steps, msgs);
  // y = y0; model.solve(y, t0, ts[0], rate, integ1);
  // y1 = stan::math::integrate_ode_rk45(f1, yvec, t0, ts, f1.adaptor.adapted_param(), {}, x_i, msgs); // NOLINT
  // EXPECT_FLOAT_EQ(y(0), y1[0][0]);
  // EXPECT_FLOAT_EQ(y(1), y1[0][1]);
  // EXPECT_FLOAT_EQ(y(2), y1[0][2]);

  // torsten::dsolve::PMXOdeIntegrator<torsten::StanAdams> integ2(rtol, atol, max_num_steps, msgs);
  // y = y0; model.solve(y, t0, ts[0], rate, integ2);
  // y1 = stan::math::integrate_ode_adams(f1, yvec, t0, ts, f1.adaptor.adapted_param(), {}, x_i, msgs); // NOLINT
  // EXPECT_FLOAT_EQ(y(0), y1[0][0]);
  // EXPECT_FLOAT_EQ(y(1), y1[0][1]);
  // EXPECT_FLOAT_EQ(y(2), y1[0][2]);

  // torsten::dsolve::PMXOdeIntegrator<torsten::StanBdf> integ3(rtol, atol, max_num_steps, msgs);
  // y = y0; model.solve(y, t0, ts[0], rate, integ3);
  // y1 = stan::math::integrate_ode_bdf(f1, yvec, t0, ts, f1.adaptor.adapted_param(), {}, x_i, msgs); // NOLINT
  // EXPECT_FLOAT_EQ(y(0), y1[0][0]);
  // EXPECT_FLOAT_EQ(y(1), y1[0][1]);
  // EXPECT_FLOAT_EQ(y(2), y1[0][2]);

  torsten::dsolve::PMXOdeIntegrator<PMXVariadicOdeSystem,  PMXCvodesIntegrator<CV_ADAMS, CV_STAGGERED>> integ4(rtol, atol, max_num_steps, msgs);
  y = y0; model.solve(y, t0, ts[0], rate, integ4);
  y1 = torsten::pmx_ode_adams(f1, y0, t0, ts, msgs, model0.par(), rate, x_r, x_i);
  EXPECT_FLOAT_EQ(y(0), y1[0][0]);
  EXPECT_FLOAT_EQ(y(1), y1[0][1]);
  EXPECT_FLOAT_EQ(y(2), y1[0][2]);

  PMXOdeIntegrator<PMXVariadicOdeSystem, PMXCvodesIntegrator<CV_BDF, CV_STAGGERED>> integ5(rtol, atol, max_num_steps, msgs);
  y = y0; model.solve(y, t0, ts[0], rate, integ5);
  y1 = torsten::pmx_ode_bdf(f1, y0, t0, ts, msgs, model0.par(), rate, x_r, x_i);
  EXPECT_FLOAT_EQ(y(0), y1[0][0]);
  EXPECT_FLOAT_EQ(y(1), y1[0][1]);
  EXPECT_FLOAT_EQ(y(2), y1[0][2]);
}

TEST_F(TorstenTwoCptModelTest, general_ode_solver_y0) {
  rate[0] = 1200;
  rate[1] = 2000;
  rate[2] = 3000;
  y0[0] = 800;
  y0[1] = 0;
  y0[2] = 8000;
  PMXTwoCptModel<double> model0(CL, Q, V2, V3, ka); // NOLINT
  PMXOdeFunctorRateAdaptor<PMXTwoCptODE> f1(f0);
  using model_t = PKODEModel<double, PMXTwoCptODE>;
  model_t model(model0.par(), y0.size(), model0.f());

  Eigen::Matrix<double, -1, 1> y;
  std::vector<Eigen::VectorXd> y1;
  ts[0] = 20.0;
  ts.resize(1);

  // torsten::dsolve::PMXOdeIntegrator<torsten::StanRk45> integ1(rtol, atol, max_num_steps, msgs);
  // y = y0; model.solve(y, t0, ts[0], rate, integ1);
  // y1 = stan::math::integrate_ode_rk45(f1, yvec, t0, ts, f1.adaptor.adapted_param(), {}, x_i, msgs); // NOLINT
  // EXPECT_FLOAT_EQ(y(0), y1[0][0]);
  // EXPECT_FLOAT_EQ(y(1), y1[0][1]);
  // EXPECT_FLOAT_EQ(y(2), y1[0][2]);

  // torsten::dsolve::PMXOdeIntegrator<torsten::StanAdams> integ2(rtol, atol, max_num_steps, msgs);
  // y = y0; model.solve(y, t0, ts[0], rate, integ2);
  // y1 = stan::math::integrate_ode_adams(f1, yvec, t0, ts, f1.adaptor.adapted_param(), {}, x_i, msgs); // NOLINT
  // EXPECT_FLOAT_EQ(y(0), y1[0][0]);
  // EXPECT_FLOAT_EQ(y(1), y1[0][1]);
  // EXPECT_FLOAT_EQ(y(2), y1[0][2]);

  // torsten::dsolve::PMXOdeIntegrator<torsten::StanBdf> integ3(rtol, atol, max_num_steps, msgs);
  // y = y0; model.solve(y, t0, ts[0], rate, integ3);
  // y1 = stan::math::integrate_ode_bdf(f1, yvec, t0, ts, f1.adaptor.adapted_param(), {}, x_i, msgs); // NOLINT
  // EXPECT_FLOAT_EQ(y(0), y1[0][0]);
  // EXPECT_FLOAT_EQ(y(1), y1[0][1]);
  // EXPECT_FLOAT_EQ(y(2), y1[0][2]);

  torsten::dsolve::PMXOdeIntegrator<PMXVariadicOdeSystem,  PMXCvodesIntegrator<CV_ADAMS, CV_STAGGERED>> integ4(rtol, atol, max_num_steps, msgs);
  y = y0; model.solve(y, t0, ts[0], rate, integ4);
  y1 = torsten::pmx_ode_adams(f1, y0, t0, ts, msgs, model0.par(), rate, x_r, x_i);
  EXPECT_FLOAT_EQ(y(0), y1[0][0]);
  EXPECT_FLOAT_EQ(y(1), y1[0][1]);
  EXPECT_FLOAT_EQ(y(2), y1[0][2]);

  PMXOdeIntegrator<PMXVariadicOdeSystem, PMXCvodesIntegrator<CV_BDF, CV_STAGGERED>> integ5(rtol, atol, max_num_steps, msgs);
  y = y0; model.solve(y, t0, ts[0], rate, integ5);
  y1 = torsten::pmx_ode_bdf(f1, y0, t0, ts, msgs, model0.par(), rate, x_r, x_i);
  EXPECT_FLOAT_EQ(y(0), y1[0][0]);
  EXPECT_FLOAT_EQ(y(1), y1[0][1]);
  EXPECT_FLOAT_EQ(y(2), y1[0][2]);
}

TEST_F(TorstenTwoCptModelTest, general_ode_solver_par_sens) {
  rate[0] = 1200;
  rate[1] = 2000;
  rate[2] = 3000;
  y0[0] = 800;
  y0[1] = 0;
  y0[2] = 8000;
  std::vector<stan::math::var> theta = to_var(par);

  PMXTwoCptModel<var> model0(theta);
  PMXOdeFunctorRateAdaptor<PMXTwoCptODE> f1(f0);
  using model_t = PKODEModel<var, PMXTwoCptODE>;
  model_t model(model0.par(), y0.size(), model0.f());

  Eigen::Matrix<var, Eigen::Dynamic, 1> y;
  std::vector<Eigen::Matrix<var, -1, 1> > y1;
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

  // torsten::dsolve::PMXOdeIntegrator<torsten::StanRk45> integ1(rtol, atol, max_num_steps, msgs);
  // torsten::dsolve::PMXOdeIntegrator<torsten::StanAdams> integ2(rtol, atol, max_num_steps, msgs);
  // torsten::dsolve::PMXOdeIntegrator<torsten::StanBdf> integ3(rtol, atol, max_num_steps, msgs);
  PMXOdeIntegrator<PMXVariadicOdeSystem, PMXCvodesIntegrator<CV_ADAMS, CV_STAGGERED>> integ4(rtol, atol, max_num_steps, msgs);
  PMXOdeIntegrator<PMXVariadicOdeSystem, PMXCvodesIntegrator<CV_BDF, CV_STAGGERED>> integ5(rtol, atol, max_num_steps, msgs);

  // y = to_var(y0); model.solve(y, t0, ts[0], rate, integ1);
  // y1 = stan::math::integrate_ode_rk45(f1, yvec, t0, ts, f1.adaptor.adapted_param(), {}, x_i, msgs); // NOLINT
  // test_it();

  // y = to_var(y0); model.solve(y, t0, ts[0], rate, integ2);
  // y1 = stan::math::integrate_ode_adams(f1, yvec, t0, ts, f1.adaptor.adapted_param(), {}, x_i, msgs); // NOLINT
  // test_it();

  // y = to_var(y0); model.solve(y, t0, ts[0], rate, integ3);
  // y1 = stan::math::integrate_ode_bdf(f1, yvec, t0, ts, f1.adaptor.adapted_param(), {}, x_i, msgs); // NOLINT
  // test_it();

  y = to_var(y0); model.solve(y, t0, ts[0], rate, integ4);
  y1 = torsten::pmx_ode_adams(f1, y0, t0, ts, msgs, model0.par(), rate, x_r, x_i);
  test_it();

  y = to_var(y0); model.solve(y, t0, ts[0], rate, integ5);
  y1 = torsten::pmx_ode_bdf(f1, y0, t0, ts, msgs, model0.par(), rate, x_r, x_i);
  test_it();
}

TEST_F(TorstenTwoCptModelTest, general_ode_solver_par_rate_sens) {
  rate[0] = 1200;
  rate[1] = 2000;
  rate[2] = 3000;
  y0[0] = 800;
  y0[1] = 0;
  y0[2] = 8000;
  std::vector<stan::math::var> theta = to_var(par);
  std::vector<stan::math::var> rate_var = to_var(rate);

  PMXTwoCptModel<var> model0(theta);
  PMXOdeFunctorRateAdaptor<PMXTwoCptODE> f1(f0);
  using model_t = PKODEModel<var, PMXTwoCptODE>;
  model_t model(theta, y0.size(), model0.f());

  Eigen::Matrix<var, Eigen::Dynamic, 1> y;
  std::vector<Eigen::Matrix<var, -1, 1> > y1;
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

  // torsten::dsolve::PMXOdeIntegrator<torsten::StanRk45> integ1(rtol, atol, max_num_steps, msgs);
  // torsten::dsolve::PMXOdeIntegrator<torsten::StanAdams> integ2(rtol, atol, max_num_steps, msgs);
  // torsten::dsolve::PMXOdeIntegrator<torsten::StanBdf> integ3(rtol, atol, max_num_steps, msgs);
  PMXOdeIntegrator<PMXVariadicOdeSystem, PMXCvodesIntegrator<CV_ADAMS, CV_STAGGERED>> integ4(rtol, atol, max_num_steps, msgs);
  PMXOdeIntegrator<PMXVariadicOdeSystem, PMXCvodesIntegrator<CV_BDF, CV_STAGGERED>> integ5(rtol, atol, max_num_steps, msgs);

  // y = to_var(y0); model.solve(y, t0, ts[0], rate_var, integ1);
  // y1 = stan::math::integrate_ode_rk45(f1, yvec, t0, ts, theta, x_r, x_i, msgs); // NOLINT
  // test_it();

  // y = to_var(y0); model.solve(y, t0, ts[0], rate_var, integ2);
  // y1 = stan::math::integrate_ode_adams(f1, yvec, t0, ts, theta, x_r, x_i, msgs); // NOLINT
  // test_it();

  // y = to_var(y0); model.solve(y, t0, ts[0], rate_var, integ3);
  // y1 = stan::math::integrate_ode_bdf(f1, yvec, t0, ts, theta, x_r, x_i, msgs); // NOLINT
  // test_it();

  y = to_var(y0); model.solve(y, t0, ts[0], rate_var, integ4);
  y1 = torsten::pmx_ode_adams(f1, y0, t0, ts, msgs, theta, rate_var, x_r, x_i);
  test_it();

  y = to_var(y0); model.solve(y, t0, ts[0], rate_var, integ5);
  y1 = torsten::pmx_ode_bdf(f1, y0, t0, ts, msgs, theta, rate_var, x_r, x_i);
  test_it();
}
