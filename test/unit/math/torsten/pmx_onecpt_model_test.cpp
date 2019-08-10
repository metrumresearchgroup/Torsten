#include <stan/math.hpp>
#include <stan/math/rev/core.hpp>
#include <test/unit/math/rev/mat/fun/util.hpp>
#include <test/unit/math/torsten/pmx_onecpt_model_test_fixture.hpp>
#include <test/unit/math/torsten/test_util.hpp>
#include <test/unit/util.hpp>
#include <gtest/gtest.h>

using stan::math::var;
using stan::math::to_var;
using refactor::PMXOneCptModel;
using stan::math::integrate_ode_bdf;
using torsten::pmx_integrate_ode_bdf;
using refactor::PMXOneCptODE;
using refactor::PMXOdeFunctorRateAdaptor;

TEST_F(TorstenOneCptModelTest, rate_dbl) {
  rate[0] = 1200;
  rate[1] = 200;
  using model_t = PMXOneCptModel<double, double, double, double>;
  model_t model(t0, y0, rate, CL, V2, ka);
  std::vector<double> yvec(y0.data(), y0.data() + y0.size());
  PMXOdeFunctorRateAdaptor<PMXOneCptODE, double> f1(model.f());

  std::vector<double> y = f1(t0, yvec, model.par(), rate, x_i, msgs);
  EXPECT_DOUBLE_EQ(y[0], rate[0]);
  EXPECT_DOUBLE_EQ(y[1], rate[1]);
  EXPECT_FALSE(torsten::has_var_rate<model_t>::value);
}

TEST_F(TorstenOneCptModelTest, rate_var) {
  rate[0] = 1200;
  rate[1] = 200;
  var CLv = to_var(CL);
  var V2v = to_var(V2);
  var kav = to_var(ka);
  std::vector<stan::math::var> rate_var{to_var(rate)};
  using model_t = PMXOneCptModel<double, double, var, var>;
  model_t model(t0, y0, rate_var, CLv, V2v, kav);
  std::vector<double> yvec(y0.data(), y0.data() + y0.size());
  std::vector<stan::math::var> theta(model.par());
  PMXOdeFunctorRateAdaptor<PMXOneCptODE, var> f1(model.f(), theta.size());
  theta.insert(theta.end(), rate_var.begin(), rate_var.end());

  std::vector<var> y = f1(t0, yvec, theta, x_r, x_i, msgs);
  EXPECT_FLOAT_EQ(y[0].val(), rate[0]);
  EXPECT_FLOAT_EQ(y[1].val(), rate[1]);

  EXPECT_TRUE(torsten::has_var_rate<model_t>::value);
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
  using model_t = PMXOneCptModel<double, double, var, var>;
  model_t model(t0, y0, rate_var, CLv, V2v, kav);
  std::vector<stan::math::var> theta(model.par());
  std::vector<double> yvec(y0.data(), y0.data() + y0.size());
  PMXOdeFunctorRateAdaptor<PMXOneCptODE, var> f1(model.f(), theta.size());
  theta.insert(theta.end(), rate_var.begin(), rate_var.end());

  std::vector<var> y = f1(t0, yvec, theta, x_r, x_i, msgs);

  EXPECT_TRUE(torsten::has_var_rate<model_t>::value);
  EXPECT_FALSE(torsten::has_var_init<model_t>::value);
  EXPECT_TRUE(torsten::has_var_par<model_t>::value);
  EXPECT_TRUE((std::is_same<torsten::f_t<model_t>, refactor::PMXOneCptODE>::value));
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
  using model_t = PMXOneCptModel<double, double, var, var>;
  model_t model(t0, y0, rate_var, CLv, V2v, kav);
  std::vector<stan::math::var> theta(model.par());
  std::vector<double> yvec(y0.data(), y0.data() + y0.size());
  PMXOdeFunctorRateAdaptor<PMXOneCptODE, var> f1(model.f(), theta.size());
  theta.insert(theta.end(), rate_var.begin(), rate_var.end());

  auto y1 = pmx_integrate_ode_bdf(f1, yvec, t0, ts, theta, x_r, x_i, msgs);
  auto y2 = model.solve(ts[0]);

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
  
  auto f1 = [&](std::vector<double>& pars) {
    using model_t = PMXOneCptModel<double, double, double, double>;
    model_t model(t0, y0, rate, pars[0], pars[1], pars[2]);
    return model.solve(dt);
  };
  auto f2 = [&](std::vector<var>& pars) {
    using model_t = PMXOneCptModel<double, double, double, var>;
    model_t model(t0, y0, rate, pars[0], pars[1], pars[2]);
    return model.solve(dt);
  };

  std::vector<double> pars{CL, V2, ka};
  torsten::test::test_grad(f1, f2, pars, 1.e-3, 1.e-16, 1.e-6, 1.e-12);
}

TEST_F(TorstenOneCptModelTest, ss_bolus_amt_grad) {
  rate[0] = 0;
  rate[1] = 0;
  y0[0] = 150;
  y0[1] = 50;
  ts[0] = 10.0;
  ts.resize(1);
  using model_t = PMXOneCptModel<double, double, double, double>;
  model_t model(t0, y0, rate, CL, V2, ka);

  double ii = 12.0;
  
  int cmt = 0;
  auto f1 = [&](std::vector<double>& amt_vec) {
    return model.solve(amt_vec[0], rate[cmt-1], ii, cmt);
  };
  auto f2 = [&](std::vector<var>& amt_vec) {
    return model.solve(amt_vec[0], rate[cmt-1], ii, cmt);
  };
  std::vector<double> amt_vec{1000.0};
  cmt = 1; torsten::test::test_grad(f1, f2, amt_vec, 1.e-3, 1.e-16, 2.e-10, 1.e-12);
  cmt = 2; torsten::test::test_grad(f1, f2, amt_vec, 1.e-3, 1.e-16, 2.e-10, 1.e-12);
}

TEST_F(TorstenOneCptModelTest, ss_infusion_rate_grad) {
  y0[0] = 150;
  y0[1] = 50;
  ts[0] = 10.0;
  ts.resize(1);
  using model_t = PMXOneCptModel<double, double, double, double>;
  model_t model(t0, y0, rate, CL, V2, ka);

  double amt = 1800;
  double ii = 12.0;
  
  int cmt = 0;
  auto f1 = [&](std::vector<double>& rate_vec) {
    return model.solve(amt, rate_vec[0], ii, cmt);
  };
  auto f2 = [&](std::vector<var>& rate_vec) {
    return model.solve(amt, rate_vec[0], ii, cmt);
  };

  std::vector<double> rate_vec{500.0};
  cmt = 1; torsten::test::test_grad(f1, f2, rate_vec, 1.e-3, 1.e-15, 2.e-10, 1.e-10);
  cmt = 2; torsten::test::test_grad(f1, f2, rate_vec, 1.e-3, 1.e-15, 2.e-10, 1.e-10);
}

TEST_F(TorstenOneCptModelTest, ss_bolus_grad_vs_long_run_sd) {
  rate[0] = 0;
  rate[1] = 0;
  y0[0] = 150;
  y0[1] = 50;
  ts[0] = 10.0;
  ts.resize(1);
  using model_t = PMXOneCptModel<double, double, double, double>;
  model_t model(t0, y0, rate, CL, V2, ka);

  int cmt = 0;
  double ii = 12.0;
  
  auto f1 = [&](std::vector<double>& amt_vec) {
    double t = t0;
    Eigen::Matrix<double, -1, 1> y = y0;
    for (int i = 0; i < 100; ++i) {
      Eigen::Matrix<double, 1, -1> yt = y.transpose();
      model_t model_i(t, yt, rate, CL, V2, ka);
      double t_next = t + ii;
      Eigen::Matrix<double, -1, 1> ys = model_i.solve(t_next);
      ys(cmt - 1) += amt_vec[0];
      y = ys;
      t = t_next;
    }
    // steady state solution is the end of II dosing before
    // bolus is imposed, to check that we remove the
    // bolus(added in the last iteration) from the results
    y(cmt - 1) -= amt_vec[0];
    return y;
  };
  auto f2 = [&](std::vector<var>& amt_vec) {
    return model.solve(amt_vec[0], rate[cmt - 1], ii, cmt);
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
  using model_t = PMXOneCptModel<var, var, var, double>;

  int cmt = 0;
  double ii = 6.0;
  double amt = 1000;
  
  auto f1 = [&](std::vector<var>& rate_vec) {
    var t = t0;
    Eigen::Matrix<var, -1, 1> y = y0;
    var t_infus = amt/rate_vec[cmt - 1];
    const std::vector<var> rate_zero{0.0, 0.0};
    for (int i = 0; i < 200; ++i) {
      Eigen::Matrix<var, 1, -1> yt;
      Eigen::Matrix<var, -1, 1> ys;
      yt = y.transpose();
      model_t model_i(t, yt, rate_vec, CL, V2, ka);
      var t_next = t + t_infus;
      ys = model_i.solve(t_next);
      yt = ys.transpose();
      t = t_next;
      model_t model_j(t, yt, rate_zero, CL, V2, ka);
      t_next = t + ii - t_infus;
      ys = model_j.solve(t_next);
      y = ys;
    }
    return y;
  };

  auto f2 = [&](std::vector<var>& rate_vec) {
    PMXOneCptModel<double, double, double, double> model(t0, y0, rate, CL, V2, ka);
    return model.solve(amt, rate_vec[cmt - 1], ii, cmt);
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
  using model_t = refactor::PMXOneCptModel<double, double, double, double>;

  int cmt = 0;
  const double ii = 11.9;
  const double amt = 1000;
  
  auto f1 = [&](std::vector<double>& rate_vec) {
    double t = t0;
    Eigen::Matrix<double, -1, 1> y = y0;
    double t_infus = amt/rate_vec[cmt - 1];
    const std::vector<double> rate_zero{0.0, 0.0};
    for (int i = 0; i < 50; ++i) {
      Eigen::Matrix<double, 1, -1> yt;
      Eigen::Matrix<double, -1, 1> ys;
      yt = y.transpose();
      model_t model_i(t, yt, rate_vec, CL, V2, ka);
      double t_next = t + t_infus;
      ys = model_i.solve(t_next);
      yt = ys.transpose();
      t = t_next;
      model_t model_j(t, yt, rate_zero, CL, V2, ka);
      t_next = t + ii - t_infus;
      ys = model_j.solve(t_next);
      y = ys;
    }
    return y;
  };

  PMXOneCptODE f2cpt;
  const std::vector<double> theta{CL, V2, ka};
  const PMXOdeIntegrator<PkBdf> integrator;
  using ode_model_t = refactor::PKODEModel<double, double, double, double, PMXOneCptODE>;
  auto f2 = [&](std::vector<double>& rate_vec) {
    double t = t0;
    Eigen::Matrix<double, -1, 1> y = y0;
    double t_infus = amt/rate_vec[cmt - 1];
    const std::vector<double> rate_zero{0.0, 0.0};
    for (int i = 0; i < 50; ++i) {
      Eigen::Matrix<double, 1, -1> yt;
      Eigen::Matrix<double, -1, 1> ys;
      yt = y.transpose();
      ode_model_t model_i(t, yt, rate_vec, theta, f2cpt);
      double t_next = t + t_infus;
      ys = model_i.solve(t_next, integrator);
      yt = ys.transpose();
      t = t_next;
      ode_model_t model_j(t, yt, rate_zero, theta, f2cpt);
      t_next = t + ii - t_infus;
      ys = model_j.solve(t_next, integrator);
      y = ys;
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
  using model_t = refactor::PMXOneCptModel<var, var, var, double>;

  int cmt = 0;
  const double ii = 11.9;
  const double amt = 1000;
  
  auto f1 = [&](std::vector<var>& rate_vec) {
    var t = t0;
    Eigen::Matrix<var, -1, 1> y = y0;
    var t_infus = amt/rate_vec[cmt - 1];
    const std::vector<var> rate_zero{0.0, 0.0};
    for (int i = 0; i < 50; ++i) {
      Eigen::Matrix<var, 1, -1> yt;
      Eigen::Matrix<var, -1, 1> ys;
      yt = y.transpose();
      model_t model_i(t, yt, rate_vec, CL, V2, ka);
      var t_next = t + t_infus;
      ys = model_i.solve(t_next);
      yt = ys.transpose();
      t = t_next;
      model_t model_j(t, yt, rate_zero, CL, V2, ka);
      t_next = t + ii - t_infus;
      ys = model_j.solve(t_next);
      y = ys;
    }
    return y;
  };

  PMXOneCptODE f2cpt;
  const std::vector<double> theta{CL, V2, ka};
  const PMXOdeIntegrator<PkBdf> integrator;
  using ode_model_t = refactor::PKODEModel<var, var, var, double, PMXOneCptODE>;
  auto f2 = [&](std::vector<var>& rate_vec) {
    var t = t0;
    Eigen::Matrix<var, -1, 1> y = y0;
    var t_infus = amt/rate_vec[cmt - 1];
    const std::vector<var> rate_zero{0.0, 0.0};
    for (int i = 0; i < 50; ++i) {
      Eigen::Matrix<var, 1, -1> yt;
      Eigen::Matrix<var, -1, 1> ys;
      yt = y.transpose();
      ode_model_t model_i(t, yt, rate_vec, theta, f2cpt);
      var t_next = t + t_infus;
      ys = model_i.solve(t_next, integrator);
      yt = ys.transpose();
      t = t_next;
      ode_model_t model_j(t, yt, rate_zero, theta, f2cpt);
      t_next = t + ii - t_infus;
      ys = model_j.solve(t_next, integrator);
      y = ys;
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
  using model_t = refactor::PMXOneCptModel<double, var, double, double>;

  int cmt = 0;
  const double ii = 11.9;
  
  auto f1 = [&](std::vector<var>& amt_vec) {
    double t = t0;
    Eigen::Matrix<var, -1, 1> y = y0;
    for (int i = 0; i < 100; ++i) {
      Eigen::Matrix<var, 1, -1> yt = y.transpose();
      model_t model_i(t, yt, rate, CL, V2, ka);
      double t_next = t + ii;
      Eigen::Matrix<var, -1, 1> ys = model_i.solve(t_next);
      ys(cmt - 1) += amt_vec[0];
      y = ys;
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
  const PMXOdeIntegrator<PkBdf> integrator;
  using ode_model_t = refactor::PKODEModel<double, var, double, double, PMXOneCptODE>;
  auto f2 = [&](std::vector<var>& amt_vec) {
    double t = t0;
    Eigen::Matrix<var, -1, 1> y = y0;
    for (int i = 0; i < 100; ++i) {
      Eigen::Matrix<var, 1, -1> yt = y.transpose();
      ode_model_t model_i(t, yt, rate, theta, f1cpt);
      double t_next = t + ii;
      Eigen::Matrix<var, -1, 1> ys = model_i.solve(t_next, integrator);
      ys(cmt - 1) += amt_vec[0];
      y = ys;
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
  using model_t = PMXOneCptModel<double, double, var, var>;
  model_t model(t0, y0, rate_var, CLv, V2v, kav);
  std::vector<stan::math::var> theta(model.par());

  std::cout.precision(12);

  double amt = 1800;
  int cmt = 1;
  double ii = 12.0;
  
  auto y1 = model.solve(amt, rate_var[cmt - 1], ii, cmt);
  EXPECT_FLOAT_EQ(y1(0).val(), 1.00330322392E-3);
  EXPECT_FLOAT_EQ(y1(1).val(), 2.07672937446E+0);

  std::vector<double> g1, g2;
  stan::math::set_zero_all_adjoints();
  y1(0).grad(theta, g1);
  EXPECT_FLOAT_EQ(g1[0], 0.0);
  EXPECT_FLOAT_EQ(g1[1], 0.0);
  EXPECT_FLOAT_EQ(g1[2], -0.0120396453978);
  stan::math::set_zero_all_adjoints();
  y1(1).grad(theta, g1);
  EXPECT_FLOAT_EQ(g1[0], -0.266849753086);
  EXPECT_FLOAT_EQ(g1[1], 0.166781095679);
  EXPECT_FLOAT_EQ(g1[2], -1.8559692314);

  cmt = 2;
  y1 = model.solve(amt, rate_var[cmt - 1], ii, cmt);
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
  using model_t = PMXOneCptModel<double, double, var, var>;
  model_t model(t0, y0, rate_var, CLv, V2v, kav);
  std::vector<var> theta{model.par()};

  double amt = 1800;
  int cmt = 1;
  double ii = 12.0;
  
  auto y1 = model.solve(amt, rate_var[cmt - 1], ii, cmt);
  EXPECT_FLOAT_EQ(y1(0).val(), 0.00312961339574);
  EXPECT_FLOAT_EQ(y1(1).val(), 3.61310672484);

  std::vector<double> g1, g2;
  stan::math::set_zero_all_adjoints();
  y1(0).grad(theta, g1);
  EXPECT_FLOAT_EQ(g1[0], 0.0);
  EXPECT_FLOAT_EQ(g1[1], 0.0);
  EXPECT_FLOAT_EQ(g1[2], -0.0342061212685);
  stan::math::set_zero_all_adjoints();
  y1(1).grad(theta, g1);
  EXPECT_FLOAT_EQ(g1[0], -0.421478621857);
  EXPECT_FLOAT_EQ(g1[1], 0.26342413866);
  EXPECT_FLOAT_EQ(g1[2], -3.20135491073);

  cmt = 2;
  y1 = model.solve(amt, rate_var[cmt - 1], ii, cmt);
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
  using model_t = PMXOneCptModel<double, double, var, var>;
  model_t model(t0, y0, rate_var, CLv, V2v, kav);

  std::cout.precision(12);

  double amt = 1800;
  int cmt = 1;
  double ii = 0.0;
  
  auto y1 = model.solve(amt, rate_var[cmt - 1], ii, cmt);
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
  EXPECT_NEAR(g1[2], 0.0, 5e-13);

  cmt = 2;
  y1 = model.solve(amt, rate_var[cmt - 1], ii, cmt);
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
