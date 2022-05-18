#include <stan/math.hpp>
#include <stan/math/rev/core.hpp>
#include <stan/math/torsten/dsolve/pmx_ode_integrator.hpp>
#include <test/unit/math/rev/fun/util.hpp>
#include <stan/math/torsten/test/unit/pmx_onecpt_model_test_fixture.hpp>
#include <stan/math/torsten/test/unit/test_macros.hpp>
#include <test/unit/util.hpp>
#include <gtest/gtest.h>

using stan::math::var;
using stan::math::to_var;
using torsten::PMXOneCptModel;
using stan::math::integrate_ode_bdf;
using torsten::pmx_integrate_ode_bdf;
using torsten::PMXOneCptODE;
using torsten::PKODEModel;
using torsten::PMXOdeFunctorRateAdaptor;
using torsten::dsolve::PMXOdeIntegrator;
using torsten::dsolve::PMXVariadicOdeSystem;
using torsten::dsolve::PMXCvodesIntegrator;
using torsten::dsolve::PMXOdeintIntegrator;

PMXOneCptODE f0;

TEST_F(TorstenOneCptModelTest, ka_zero) {
  y0(0) = 745;
  y0(1) = 100;
  rate[0] = 1200;
  rate[1] = 200;
  ka = 0.0;
  using model_t = PMXOneCptModel<double>;
  model_t model(CL, V2, ka);
  torsten::PKRec<double> y(y0);
  model.solve(y, t0, ts[0], rate);
  EXPECT_FLOAT_EQ(y(0), 865.0);
  EXPECT_FLOAT_EQ(y(1), 113.32912);
}

TEST_F(TorstenOneCptModelTest, rate_dbl) {
  rate[0] = 1200;
  rate[1] = 200;
  using model_t = PMXOneCptModel<double>;
  model_t model(CL, V2, ka);
  torsten::PMXOdeFunctorRateAdaptor<PMXOneCptODE> f1(f0);

  Eigen::VectorXd y = f1(t0, y0, msgs, model.par(), rate, x_r, x_i);
  EXPECT_DOUBLE_EQ(y[0], rate[0]);
  EXPECT_DOUBLE_EQ(y[1], rate[1]);
}

TEST_F(TorstenOneCptModelTest, rate_var) {
  rate[0] = 1200;
  rate[1] = 200;
  var CLv = to_var(CL);
  var V2v = to_var(V2);
  var kav = to_var(ka);
  std::vector<stan::math::var> rate_var{to_var(rate)};
  using model_t = PMXOneCptModel<var>;
  model_t model(CLv, V2v, kav);
  std::vector<stan::math::var> theta(model.par());
  torsten::PMXOdeFunctorRateAdaptor<PMXOneCptODE> f1(f0);

  Eigen::Matrix<var, -1, 1> y = f1(t0, y0, msgs, theta, rate_var, x_r, x_i);
  EXPECT_FLOAT_EQ(y[0].val(), rate[0]);
  EXPECT_FLOAT_EQ(y[1].val(), rate[1]);
}

TEST_F(TorstenOneCptModelTest, rate_var_y0) {
  rate[0] = 1200;
  rate[1] = 200;
  y0[0] = 150;
  y0[1] = 50;
  var CLv = to_var(CL);
  var V2v = to_var(V2);
  var kav = to_var(ka);
  std::vector<stan::math::var> rate_var{to_var(rate)};
  using model_t = PMXOneCptModel<var>;
  model_t model(CLv, V2v, kav);
  std::vector<stan::math::var> theta(model.par());
  torsten::PMXOdeFunctorRateAdaptor<PMXOneCptODE> f1(f0);

  Eigen::Matrix<var, -1, 1> y = f1(t0, y0, msgs, theta, rate_var, x_r, x_i);

  EXPECT_TRUE((std::is_same<torsten::f_t<model_t>, torsten::PMXOneCptODE>::value));
}

TEST_F(TorstenOneCptModelTest, sd_solver) {
  rate[0] = 1200;
  rate[1] = 200;
  y0[0] = 150;
  y0[1] = 50;
  ts[0] = 10.0;
  ts.resize(1);
  var CLv = to_var(CL);
  var V2v = to_var(V2);
  var kav = to_var(ka);
  std::vector<stan::math::var> rate_var{to_var(rate)};
  using model_t = PMXOneCptModel<var>;
  model_t model(CLv, V2v, kav);
  std::vector<stan::math::var> theta(model.par());
  torsten::PMXOdeFunctorRateAdaptor<PMXOneCptODE> f1(f0);

  auto y1 = pmx_ode_bdf(f1, y0, t0, ts, msgs, theta, rate_var, x_r, x_i);
  torsten::PKRec<var> y2(to_var(y0));
  model.solve(y2, t0, ts[0], rate_var);

  stan::math::vector_v y1_v = stan::math::to_vector(y1[0]);

  torsten::test::test_grad(theta, y1_v, y2, 1.e-6, 1.e-6);
  torsten::test::test_grad(rate_var, y1_v, y2, 1.e-6, 1.e-8);
}

TEST_F(TorstenOneCptModelTest, infusion_theta_grad) {
  rate[0] = 1100;
  rate[1] = 810;
  y0[0] = 150;
  y0[1] = 50;

  double dt = 2.5;
  
  auto f1 = [&](const std::vector<double>& pars) {
    using model_t = PMXOneCptModel<double>;
    model_t model(pars[0], pars[1], pars[2]);
    torsten::PKRec<double> y(y0);
    model.solve(y, t0, dt, rate);
    return y;
  };
  auto f2 = [&](const std::vector<double>& pars) {
    using model_t = PMXOneCptModel<double>;
    model_t model(pars[0], pars[1], pars[2]);
    torsten::PKRec<double> y(y0);
    model.solve_analytical(y, t0, dt, rate);
    return y;
  };
  auto f3 = [&](const std::vector<var>& pars) {
    using model_t = PMXOneCptModel<var>;
    model_t model(pars[0], pars[1], pars[2]);
    torsten::PKRec<var> y(to_var(y0));
    model.solve(y, t0, dt, rate);
    return y;
  };

  std::vector<double> pars{CL, V2, ka};
  torsten::test::test_grad(f1, f3, pars, 1.e-3, 1.e-16, 1.e-6, 2.e-10);
  torsten::test::test_grad(f2, f3, pars, 1.e-3, 1.e-12, 1.e-6, 2.e-10);
}

TEST_F(TorstenOneCptModelTest, ss_bolus_amt_grad) {
  rate[0] = 0;
  rate[1] = 0;
  y0[0] = 150;
  y0[1] = 50;
  ts[0] = 10.0;
  ts.resize(1);
  using model_t = PMXOneCptModel<double>;
  model_t model(CL, V2, ka);

  double ii = 12.0;
  
  int cmt = 0;
  auto f1 = [&](const std::vector<double>& amt_vec) {
    return model.solve(t0, amt_vec[0], rate[cmt-1], ii, cmt);
  };
  auto f2 = [&](const std::vector<var>& amt_vec) {
    return model.solve(t0, amt_vec[0], rate[cmt-1], ii, cmt);
  };
  std::vector<double> amt_vec{1000.0};
  cmt = 1; torsten::test::test_grad(f1, f2, amt_vec, 1.e-3, 1.e-16, 2.e-10, 1.e-12);
  cmt = 2; torsten::test::test_grad(f1, f2, amt_vec, 1.e-3, 1.e-16, 2.e-10, 1.e-12);

  // compare against textbook analytical sol.
  auto f3 = [&](const std::vector<var>& amt_vec) {
    return model.solve_analytical(t0, amt_vec[0], rate[cmt-1], ii, cmt);
  };
  torsten::PKRec<var> y2, y3;
  std::vector<var> amt_vec_var = stan::math::to_var(amt_vec);
  for (int i = 0; i < 2; ++i) {
    cmt = i + 1;
    y2 = f2(amt_vec_var);
    y3 = f3(amt_vec_var);
    torsten::test::test_grad(amt_vec_var, y2, y3, 1.e-18, 1.e-18);
  }
}

TEST_F(TorstenOneCptModelTest, ss_infusion_rate_grad) {
  y0[0] = 150;
  y0[1] = 50;
  ts[0] = 10.0;
  ts.resize(1);
  using model_t = PMXOneCptModel<double>;
  model_t model(CL, V2, ka);

  double amt = 1800;
  double ii = 12.0;
  
  int cmt = 0;
  auto f1 = [&](const std::vector<double>& rate_vec) {
    return model.solve(t0, amt, rate_vec[0], ii, cmt);
  };
  auto f2 = [&](const std::vector<var>& rate_vec) {
    return model.solve(t0, amt, rate_vec[0], ii, cmt);
  };

  std::vector<double> rate_vec{500.0};
  cmt = 1; torsten::test::test_grad(f1, f2, rate_vec, 1.e-3, 1.e-15, 2.e-10, 1.e-10);
  cmt = 2; torsten::test::test_grad(f1, f2, rate_vec, 1.e-3, 1.e-15, 2.e-10, 1.e-10);

  // compare against textbook analytical sol.
  auto f3 = [&](const std::vector<var>& rate_vec) {
    return model.solve_analytical(t0, amt, rate_vec[0], ii, cmt);
  };
  torsten::PKRec<var> y2, y3;
  std::vector<var> rate_var = stan::math::to_var(rate_vec);
  for (int i = 0; i < 2; ++i) {
    cmt = i + 1;
    y2 = f2(rate_var);
    y3 = f3(rate_var);
    torsten::test::test_grad(rate_var, y2, y3, 1.e-15, 1.e-15);
  }
}

TEST_F(TorstenOneCptModelTest, ss_bolus_grad_vs_long_run_sd) {
  rate[0] = 0;
  rate[1] = 0;
  y0[0] = 150;
  y0[1] = 50;
  ts[0] = 10.0;
  ts.resize(1);
  using model_t = PMXOneCptModel<double>;
  model_t model(CL, V2, ka);

  int cmt = 0;
  double ii = 12.0;
  
  auto f1 = [&](const std::vector<double>& amt_vec) {
    double t = t0;
    Eigen::Matrix<double, -1, 1> y = y0;
    for (int i = 0; i < 100; ++i) {
      model_t model_i(CL, V2, ka);
      double t_next = t + ii;
      model_i.solve(y, t, t_next, rate);
      y(cmt - 1) += amt_vec[0];
      t = t_next;
    }
    // steady state solution is the end of II dosing before
    // bolus is imposed, to check that we remove the
    // bolus(added in the last iteration) from the results
    y(cmt - 1) -= amt_vec[0];
    return y;
  };
  auto f2 = [&](const std::vector<var>& amt_vec) {
    return model.solve(t0, amt_vec[0], rate[cmt - 1], ii, cmt);
  };
  std::vector<double> amt_vec{1000.0};

  cmt = 1;
  torsten::test::test_grad(f1, f2, amt_vec, 1.e-3, 1.e-12, 1.e-8, 1.e-10);
  cmt = 2;
  torsten::test::test_grad(f1, f2, amt_vec, 1.e-3, 1.e-12, 1.e-7, 1.e-10);
}

TEST_F(TorstenOneCptModelTest, ss_infusion_grad_vs_long_run_sd) {
  y0[0] = 150;
  y0[1] = 50;
  using model_t = PMXOneCptModel<double>;

  int cmt = 0;
  double ii = 6.0;
  double amt = 1000;
  
  auto f1 = [&](std::vector<var>& rate_vec) {
    var t = t0;
    Eigen::Matrix<var, -1, 1> y = y0;
    var t_infus = amt/rate_vec[cmt - 1];
    const std::vector<var> rate_zero{0.0, 0.0};
    for (int i = 0; i < 200; ++i) {
      model_t model_i(CL, V2, ka);
      var t_next = t + t_infus;
      model_i.solve(y, t, t_next, rate_vec);
      t = t_next;
      model_t model_j(CL, V2, ka);
      t_next = t + ii - t_infus;
      model_j.solve(y, t, t_next, rate_zero);
    }
    return y;
  };

  auto f2 = [&](std::vector<var>& rate_vec) {
    PMXOneCptModel<double> model(CL, V2, ka);
    return model.solve(t0, amt, rate_vec[cmt - 1], ii, cmt);
  };

  std::vector<var> rate_vec(2, 0.0);

  {
    cmt = 1;
    rate_vec[cmt-1] = 300.0;
    // For SS we only compare gradient wrt the dosing compartment as the
    // non-dosing compartment doesn't enter the SS system
    Eigen::Matrix<var, -1, 1> y1 = f1(rate_vec);
    Eigen::Matrix<var, -1, 1> y2 = f2(rate_vec);
    std::vector<var> par{rate_vec[cmt - 1]};
    torsten::test::test_grad(par, y1, y2, 1.e-11, 1.e-10);
    rate_vec[cmt-1] = 0.0;
  }

  {
    cmt = 2;
    rate_vec[cmt-1] = 300.0;
    // For SS we only compare gradient wrt the dosing compartment as the
    // non-dosing compartment doesn't enter the SS system
    Eigen::Matrix<var, -1, 1> y1 = f1(rate_vec);
    Eigen::Matrix<var, -1, 1> y2 = f2(rate_vec);
    std::vector<var> par{rate_vec[cmt - 1]};
    torsten::test::test_grad(par, y1, y2, 1.e-11, 1.e-10);
    rate_vec[cmt-1] = 0.0;
  }
}

TEST_F(TorstenOneCptModelTest, ss_infusion_by_long_run_sd_vs_bdf_result) {
  y0[0] = 150;
  y0[1] = 50;
  using model_t = torsten::PMXOneCptModel<double>;

  int cmt = 0;
  const double ii = 11.9;
  const double amt = 1000;
  
  auto f1 = [&](std::vector<double>& rate_vec) {
    double t = t0;
    Eigen::Matrix<double, -1, 1> y = y0;
    double t_infus = amt/rate_vec[cmt - 1];
    const std::vector<double> rate_zero{0.0, 0.0};
    for (int i = 0; i < 50; ++i) {
      model_t model_i(CL, V2, ka);
      double t_next = t + t_infus;
      model_i.solve(y, t, t_next, rate_vec);
      t = t_next;
      model_t model_j(CL, V2, ka);
      t_next = t + ii - t_infus;
      model_j.solve(y, t, t_next, rate_zero);
    }
    return y;
  };

  PMXOneCptODE f2cpt;
  const std::vector<double> theta{CL, V2, ka};
  const PMXOdeIntegrator<PMXVariadicOdeSystem, PMXCvodesIntegrator<CV_BDF, CV_STAGGERED>> integrator;
  using ode_model_t = torsten::PKODEModel<double, PMXOneCptODE>;
  auto f2 = [&](std::vector<double>& rate_vec) {
    double t = t0;
    Eigen::Matrix<double, -1, 1> y = y0;
    double t_infus = amt/rate_vec[cmt - 1];
    const std::vector<double> rate_zero{0.0, 0.0};
    for (int i = 0; i < 50; ++i) {
      ode_model_t model_i(theta, y.size(), f2cpt);
      double t_next = t + t_infus;
      model_i.solve(y, t, t_next, rate_vec, integrator);
      t = t_next;
      ode_model_t model_j(theta, y.size(), f2cpt);
      t_next = t + ii - t_infus;
      model_j.solve(y, t, t_next, rate_zero, integrator);
    }
    return y;
  };
  
  std::vector<double> rate_vec(2, 0.0);

  {
    cmt = 1;
    rate_vec[cmt - 1] = 300.0;
    Eigen::Matrix<double, -1, 1> y1 = f1(rate_vec);
    Eigen::Matrix<double, -1, 1> y2 = f2(rate_vec);
    torsten::test::test_val(y1, y2);
    rate_vec[cmt - 1] = 0.0;
  }

  {
    cmt = 2;
    rate_vec[cmt - 1] = 300.0;
    Eigen::Matrix<double, -1, 1> y1 = f1(rate_vec);
    Eigen::Matrix<double, -1, 1> y2 = f2(rate_vec);
    torsten::test::test_val(y1, y2);
  }
}

TEST_F(TorstenOneCptModelTest, ss_infusion_grad_by_long_run_sd_vs_bdf_result) {
  y0[0] = 150;
  y0[1] = 50;
  using model_t = torsten::PMXOneCptModel<double>;

  int cmt = 0;
  const double ii = 11.9;
  const double amt = 1000;
  
  auto f1 = [&](std::vector<var>& rate_vec) {
    var t = t0;
    Eigen::Matrix<var, -1, 1> y = y0;
    var t_infus = amt/rate_vec[cmt - 1];
    const std::vector<var> rate_zero{0.0, 0.0};
    for (int i = 0; i < 50; ++i) {
      model_t model_i(CL, V2, ka);
      var t_next = t + t_infus;
      model_i.solve(y, t, t_next, rate_vec);
      t = t_next;
      model_t model_j(CL, V2, ka);
      t_next = t + ii - t_infus;
      model_j.solve(y, t, t_next, rate_zero);
    }
    return y;
  };

  PMXOneCptODE f2cpt;
  const std::vector<double> theta{CL, V2, ka};
  const PMXOdeIntegrator<PMXVariadicOdeSystem, PMXCvodesIntegrator<CV_BDF, CV_STAGGERED>> integrator;
  using ode_model_t = torsten::PKODEModel<double, PMXOneCptODE>;
  auto f2 = [&](std::vector<var>& rate_vec) {
    var t = t0;
    Eigen::Matrix<var, -1, 1> y = y0;
    var t_infus = amt/rate_vec[cmt - 1];
    const std::vector<var> rate_zero{0.0, 0.0};
    for (int i = 0; i < 50; ++i) {
      ode_model_t model_i(theta, y.size(), f2cpt);
      var t_next = t + t_infus;
      model_i.solve(y, t, t_next, rate_vec, integrator);
      t = t_next;
      ode_model_t model_j(theta, y.size(), f2cpt);
      t_next = t + ii - t_infus;
      model_j.solve(y, t, t_next, rate_zero, integrator);
    }
    return y;
  };
  
  std::vector<var> rate_vec(2, 0.0);

  {
    cmt = 1;
    rate_vec[cmt - 1] = 300.0;
    Eigen::Matrix<var, -1, 1> y1 = f1(rate_vec);
    Eigen::Matrix<var, -1, 1> y2 = f2(rate_vec);
    torsten::test::test_grad(rate_vec, y1, y2, 5.e-9, 5.e-11);
    rate_vec[cmt - 1] = 0.0;
  }

  {
    cmt = 2;
    rate_vec[cmt - 1] = 300.0;
    Eigen::Matrix<var, -1, 1> y1 = f1(rate_vec);
    Eigen::Matrix<var, -1, 1> y2 = f2(rate_vec);
    torsten::test::test_grad(rate_vec, y1, y2, 5.e-8, 5.e-10);
  }
}

TEST_F(TorstenOneCptModelTest, ss_bolus_grad_by_long_run_sd_vs_bdf_result) {
  y0[0] = 150;
  y0[1] = 50;
  rate[0] = 0;
  rate[1] = 0;
  using model_t = torsten::PMXOneCptModel<double>;

  int cmt = 0;
  const double ii = 11.9;
  
  auto f1 = [&](std::vector<var>& amt_vec) {
    double t = t0;
    Eigen::Matrix<var, -1, 1> y = y0;
    for (int i = 0; i < 100; ++i) {
      model_t model_i(CL, V2, ka);
      double t_next = t + ii;
      model_i.solve(y, t, t_next, rate);
      y(cmt - 1) += amt_vec[0];
      t = t_next;
    }
    // steady state solution is the end of II dosing before
    // bolus is imposed, to check that we remove the
    // bolus(added in the last iteration) from the results
    y(cmt - 1) -= amt_vec[0];
    return y;
  };

  PMXOneCptODE f1cpt;
  const std::vector<double> theta{CL, V2, ka};
  const PMXOdeIntegrator<PMXVariadicOdeSystem, PMXCvodesIntegrator<CV_BDF, CV_STAGGERED>> integrator;
  using ode_model_t = torsten::PKODEModel<double, PMXOneCptODE>;
  auto f2 = [&](std::vector<var>& amt_vec) {
    double t = t0;
    Eigen::Matrix<var, -1, 1> y = y0;
    for (int i = 0; i < 100; ++i) {
      ode_model_t model_i(theta, y0.size(), f1cpt);
      double t_next = t + ii;
      model_i.solve(y, t, t_next, rate, integrator);
      y(cmt - 1) += amt_vec[0];
      t = t_next;
    }
    // steady state solution is the end of II dosing before
    // bolus is imposed, to check that we remove the
    // bolus(added in the last iteration) from the results
    y(cmt - 1) -= amt_vec[0];
    return y;
  };
  
  std::vector<var> amt_vec{300.0};

  {
    cmt = 1;
    Eigen::Matrix<var, -1, 1> y1 = f1(amt_vec);
    Eigen::Matrix<var, -1, 1> y2 = f2(amt_vec);
    torsten::test::test_grad(amt_vec, y1, y2, 5.e-9, 5.e-11);
  }

  {
    cmt = 2;
    Eigen::Matrix<var, -1, 1> y1 = f1(amt_vec);
    Eigen::Matrix<var, -1, 1> y2 = f2(amt_vec);
    torsten::test::test_grad(amt_vec, y1, y2, 5.e-9, 5.e-11);
  }
}

TEST_F(TorstenOneCptModelTest, ss_bolus) {
  rate[0] = 0;
  rate[1] = 0;
  y0[0] = 150;
  y0[1] = 50;
  ts[0] = 10.0;
  ts.resize(1);
  var CLv = to_var(CL);
  var V2v = to_var(V2);
  var kav = to_var(ka);
  std::vector<stan::math::var> rate_var{to_var(rate)};
  using model_t = PMXOneCptModel<var>;
  model_t model(CLv, V2v, kav);
  std::vector<stan::math::var> theta(model.par());

  std::cout.precision(12);

  double amt = 1800;
  int cmt = 1;
  double ii = 12.0;
  
  auto y1 = model.solve(t0, amt, rate_var[cmt - 1], ii, cmt);
  EXPECT_FLOAT_EQ(y1(0).val(), 1.00330322392E-3);
  EXPECT_FLOAT_EQ(y1(1).val(), 2.07672937446E+0);

  std::vector<double> g1, g2;
  stan::math::set_zero_all_adjoints();
  y1(0).grad(theta, g1);
  EXPECT_NEAR(g1[0], 0.0, 1.e-18);
  EXPECT_NEAR(g1[1], 0.0, 1.e-18);
  EXPECT_FLOAT_EQ(g1[2], -0.0120396453978);
  stan::math::set_zero_all_adjoints();
  y1(1).grad(theta, g1);
  EXPECT_FLOAT_EQ(g1[0], -0.266849753086);
  EXPECT_FLOAT_EQ(g1[1], 0.166781095679);
  EXPECT_FLOAT_EQ(g1[2], -1.8559692314);

  cmt = 2;
  y1 = model.solve(t0, amt, rate_var[cmt - 1], ii, cmt);
  EXPECT_FLOAT_EQ(y1(0).val(), 0);
  EXPECT_FLOAT_EQ(y1(1).val(), 0.996102795153);

  stan::math::set_zero_all_adjoints();
  y1(0).grad(theta, g1);
  EXPECT_FLOAT_EQ(g1[0], 0.0);
  EXPECT_FLOAT_EQ(g1[1], 0.0);
  EXPECT_FLOAT_EQ(g1[2], 0.0);
  stan::math::set_zero_all_adjoints();
  y1(1).grad(theta, g1);
  EXPECT_FLOAT_EQ(g1[0], -0.149498104338);
  EXPECT_FLOAT_EQ(g1[1], 0.0934363152112);
  EXPECT_FLOAT_EQ(g1[2], 0.0);
}

TEST_F(TorstenOneCptModelTest, ss_multi_truncated_infusion) {
  rate[0] = 1100;
  rate[1] = 770;
  y0[0] = 150;
  y0[1] = 50;
  ts[0] = 10.0;
  ts.resize(1);
  var CLv = to_var(CL);
  var V2v = to_var(V2);
  var kav = to_var(ka);
  std::vector<stan::math::var> rate_var{to_var(rate)};
  using model_t = PMXOneCptModel<var>;
  model_t model(CLv, V2v, kav);
  std::vector<var> theta{model.par()};

  double amt = 1800;
  int cmt = 1;
  double ii = 12.0;
  
  auto y1 = model.solve(t0, amt, rate_var[cmt - 1], ii, cmt);
  EXPECT_FLOAT_EQ(y1(0).val(), 0.00312961339574);
  EXPECT_FLOAT_EQ(y1(1).val(), 3.61310672484);

  std::vector<double> g1, g2;
  stan::math::set_zero_all_adjoints();
  y1(0).grad(theta, g1);
  EXPECT_NEAR(g1[0], 0.0, 1.e-16);
  EXPECT_NEAR(g1[1], 0.0, 1.e-16);
  EXPECT_FLOAT_EQ(g1[2], -0.0342061212685);
  stan::math::set_zero_all_adjoints();
  y1(1).grad(theta, g1);
  EXPECT_FLOAT_EQ(g1[0], -0.421478621857);
  EXPECT_FLOAT_EQ(g1[1], 0.26342413866);
  EXPECT_FLOAT_EQ(g1[2], -3.20135491073);

  cmt = 2;
  y1 = model.solve(t0, amt, rate_var[cmt - 1], ii, cmt);
  EXPECT_FLOAT_EQ(y1(0).val(), 0.0);
  EXPECT_FLOAT_EQ(y1(1).val(), 2.25697891686);

  stan::math::set_zero_all_adjoints();
  y1(0).grad(theta, g1);
  EXPECT_FLOAT_EQ(g1[0], 0.0);
  EXPECT_FLOAT_EQ(g1[1], 0.0);
  EXPECT_FLOAT_EQ(g1[2], 0.0);
  stan::math::set_zero_all_adjoints();
  y1(1).grad(theta, g1);
  EXPECT_FLOAT_EQ(g1[0], -0.298001025911);
  EXPECT_FLOAT_EQ(g1[1], 0.186250641194);
  EXPECT_FLOAT_EQ(g1[2], 0.0);
}

TEST_F(TorstenOneCptModelTest, ss_const_infusion) {
  rate[0] = 1100;
  rate[1] = 770;
  y0[0] = 150;
  y0[1] = 50;
  ts[0] = 10.0;
  ts.resize(1);
  var CLv = to_var(CL);
  var V2v = to_var(V2);
  var kav = to_var(ka);
  std::vector<var> theta{CLv, V2v, kav};
  std::vector<stan::math::var> rate_var{to_var(rate)};
  using model_t = PMXOneCptModel<var>;
  model_t model(CLv, V2v, kav);

  std::cout.precision(12);

  double amt = 1800;
  int cmt = 1;
  double ii = 0.0;
  
  auto y1 = model.solve(t0, amt, rate_var[cmt - 1], ii, cmt);
  EXPECT_FLOAT_EQ(y1(0).val(), 916.666666667);
  EXPECT_FLOAT_EQ(y1(1).val(), 1760);

  std::vector<double> g1, g2;
  stan::math::set_zero_all_adjoints();
  y1(0).grad(theta, g1);
  EXPECT_FLOAT_EQ(g1[0], 0.0);
  EXPECT_FLOAT_EQ(g1[1], 0.0);
  EXPECT_FLOAT_EQ(g1[2], -763.888888889);
  stan::math::set_zero_all_adjoints();
  y1(1).grad(theta, g1);
  EXPECT_FLOAT_EQ(g1[0], -35.2);
  EXPECT_FLOAT_EQ(g1[1], 22);
  // EXPECT_FLOAT_EQ(g1[2], 0.0); // FIXME: fail for g++ but clang++
  EXPECT_NEAR(g1[2], 0.0, 1e-12);

  cmt = 2;
  y1 = model.solve(t0, amt, rate_var[cmt - 1], ii, cmt);
  EXPECT_FLOAT_EQ(y1(0).val(), 0.0);
  EXPECT_FLOAT_EQ(y1(1).val(), 1232);

  stan::math::set_zero_all_adjoints();
  y1(0).grad(theta, g1);
  EXPECT_FLOAT_EQ(g1[0], 0.0);
  EXPECT_FLOAT_EQ(g1[1], 0.0);
  EXPECT_FLOAT_EQ(g1[2], 0.0);
  stan::math::set_zero_all_adjoints();
  y1(1).grad(theta, g1);
  EXPECT_FLOAT_EQ(g1[0], -24.64);
  EXPECT_FLOAT_EQ(g1[1], 15.4);
  EXPECT_FLOAT_EQ(g1[2], 0.0);
}

TEST_F(TorstenOneCptModelTest, long_ss_infusion_vs_ode) {
  double amt = 1300;
  double r = 500;
  double ii = 1.2;
  std::vector<double> par{CL, V2, ka};
  std::vector<var> par_var(to_var(par));
  PMXOneCptModel<var> model1(par_var);
  PKODEModel<var, PMXOneCptODE> model2(par_var, model1.ncmt(), model1.f());

  const dsolve::PMXAnalyiticalIntegrator integ1;
  const PMXOdeIntegrator<PMXVariadicOdeSystem, PMXCvodesIntegrator<CV_BDF, CV_STAGGERED>> integ2;
  torsten::PKRec<var> y1 = model1.solve(ts[0], amt, r, ii, 1, integ1);
  torsten::PKRec<var> y2 = model2.solve(ts[0], amt, r, ii, 1, integ2);
  torsten::test::test_grad(par_var, y1, y2, 2e-6, 6e-6);
}
