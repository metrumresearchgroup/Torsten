#include <stan/math.hpp>
#include <stan/math/rev/core.hpp>
#include <test/unit/math/rev/mat/fun/util.hpp>
#include <test/unit/math/torsten/pmx_cpt_model_test_fixture.hpp>
#include <test/unit/util.hpp>
#include <gtest/gtest.h>

using stan::math::var;
using stan::math::to_var;
using stan::math::vector_v;
using stan::math::matrix_v;
using refactor::PMXTwoCptModel;
using torsten::pmx_integrate_ode_bdf;
using stan::math::integrate_ode_bdf;
using refactor::PMXTwoCptODE;
using refactor::PMXOdeFunctorRateAdaptor;

TEST_F(TorstenTwoCptModelTest, rate_dbl) {
  rate[0] = 1200;
  rate[1] = 200;
  rate[2] = 300;
  using model_t = PMXTwoCptModel<double, double, double, double>;
  model_t model(t0, y0, rate, CL, Q, V2, V3, ka);
  std::vector<double> yvec(y0.data(), y0.data() + y0.size());
  PMXOdeFunctorRateAdaptor<PMXTwoCptODE, double> f1(model.f());

  std::vector<double> y = f1(t0, yvec, model.par(), rate, x_i, msgs);
  EXPECT_FLOAT_EQ(y[0], rate[0]);
  EXPECT_FLOAT_EQ(y[1], rate[1]);
  EXPECT_FLOAT_EQ(y[2], rate[2]);
}

TEST_F(TorstenTwoCptModelTest, rate_var) {
  rate[0] = 1200;
  rate[1] = 200;
  rate[2] = 100;
  var CLv = to_var(CL);
  var Qv  = to_var(Q);
  var V2v = to_var(V2);
  var V3v = to_var(V3);
  var kav = to_var(ka);
  std::vector<stan::math::var> rate_var{to_var(rate)};
  using model_t = PMXTwoCptModel<double, double, var, var>;
  model_t model(t0, y0, rate_var, CLv, Qv, V2v, V3v, kav);
  std::vector<stan::math::var> theta(model.par());
  std::vector<double> yvec(y0.data(), y0.data() + y0.size());
  PMXOdeFunctorRateAdaptor<PMXTwoCptODE, var> f1(model.f(), theta.size());
  theta.insert(theta.end(), rate_var.begin(), rate_var.end());

  std::vector<var> y = f1(t0, yvec, theta, x_r, x_i, msgs);
  EXPECT_FLOAT_EQ(y[0].val(), rate[0]);
  EXPECT_FLOAT_EQ(y[1].val(), rate[1]);
  EXPECT_FLOAT_EQ(y[2].val(), rate[2]);
}

TEST_F(TorstenTwoCptModelTest, sd_solver) {
  rate[0] = 1200;
  rate[1] = 200;
  rate[2] = 100;
  y0[0] = 150;
  y0[1] = 50;
  y0[2] = 60;
  ts[0] = 10.0;
  ts.resize(1);
  var CLv = to_var(CL);
  var Qv  = to_var(Q);
  var V2v = to_var(V2);
  var V3v = to_var(V3);
  var kav = to_var(ka);
  std::vector<var> theta{CLv, Qv, V2v, V3v, kav};
  std::vector<stan::math::var> rate_var{to_var(rate)};
  using model_t = PMXTwoCptModel<double, double, var, var>;
  model_t model(t0, y0, rate_var, CLv, Qv, V2v, V3v, kav);
  std::vector<double> yvec(y0.data(), y0.data() + y0.size());
  PMXOdeFunctorRateAdaptor<PMXTwoCptODE, var> f1(model.f(), theta.size());
  theta.insert(theta.end(), rate_var.begin(), rate_var.end());

  auto y1 = pmx_integrate_ode_bdf(f1, yvec, t0, ts, theta, x_r, x_i, msgs);
  auto y2 = model.solve(ts[0]);
  stan::math::vector_v y1_v = stan::math::to_vector(y1[0]);

  torsten::test::test_grad(theta, y1_v, y2, 1.e-6, 1.e-6);
  torsten::test::test_grad(rate_var, y1_v, y2, 1.e-6, 1.e-7);
}

TEST_F(TorstenTwoCptModelTest, infusion_theta_grad) {
  rate[0] = 1100;
  rate[1] = 810;
  rate[2] = 200;
  y0[0] = 150;
  y0[1] = 50.0;
  y0[2] = 50.0;

  double dt = 2.5;
  
  auto f1 = [&](std::vector<double>& pars) {
    using model_t = PMXTwoCptModel<double, double, double, double>;
    model_t model(t0, y0, rate, pars[0], pars[1], pars[2], pars[3], pars[4]);
    return model.solve(dt);
  };
  auto f2 = [&](std::vector<var>& pars) {
    using model_t = PMXTwoCptModel<double, double, double, var>;
    model_t model(t0, y0, rate, pars[0], pars[1], pars[2], pars[3], pars[4]);
    return model.solve(dt);
  };

  std::vector<double> pars{CL, Q, V2, V3, ka};
  torsten::test::test_grad(f1, f2, pars, 1.e-3, 1.e-16, 1.e-6, 1.e-12);
}

TEST_F(TorstenTwoCptModelTest, ss_bolus_amt_grad) {
  rate[0] = 0;
  rate[1] = 0;
  rate[2] = 0;
  y0[0] = 150;
  y0[1] = 250;
  y0[2] = 350;
  ts[0] = 10.0;
  ts.resize(1);
  using model_t = PMXTwoCptModel<double, double, double, double>;
  model_t model(t0, y0, rate, CL, Q, V2, V3, ka);

  double ii = 12.5;
  
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
  cmt = 3; torsten::test::test_grad(f1, f2, amt_vec, 1.e-3, 1.e-16, 2.e-10, 1.e-12);
}

TEST_F(TorstenTwoCptModelTest, ss_infusion_rate_grad) {
  y0[0] = 150;
  y0[1] = 250;
  y0[2] = 350;
  ts[0] = 10.0;
  ts.resize(1);
  using model_t = PMXTwoCptModel<double, double, double, double>;
  model_t model(t0, y0, rate, CL, Q, V2, V3, ka);

  double amt = 1100;
  double ii = 12.5;
  
  int cmt = 0;
  auto f1 = [&](std::vector<double>& rate_vec) {
    return model.solve(amt, rate_vec[0], ii, cmt);
  };
  auto f2 = [&](std::vector<var>& rate_vec) {
    return model.solve(amt, rate_vec[0], ii, cmt);
  };
  std::vector<double> rate_vec{130.0};
  cmt = 1; torsten::test::test_grad(f1, f2, rate_vec, 1.e-3, 1.e-16, 1.e-8, 1.e-12);
  cmt = 2; torsten::test::test_grad(f1, f2, rate_vec, 1.e-3, 1.e-16, 1.e-9, 1.e-12);
  cmt = 3; torsten::test::test_grad(f1, f2, rate_vec, 1.e-3, 1.e-16, 1.e-9, 1.e-12);
}

TEST_F(TorstenTwoCptModelTest, ss_bolus_by_long_run_sd_vs_bdf_result) {
  using refactor::PKODEModel;

  rate[0] = 0;
  rate[1] = 0;
  rate[2] = 0;
  y0[0] = 150;
  y0[1] = 55;
  y0[2] = 120;
  using model_t = PMXTwoCptModel<double, double, double, double>;
  model_t model(t0, y0, rate, CL, Q, V2, V3, ka);

  int cmt = 0;
  const double ii = 8.5;
  
  auto f1 = [&](std::vector<double>& amt_vec) {
    double t = t0;
    Eigen::Matrix<double, -1, 1> y = y0;
    for (int i = 0; i < 100; ++i) {
      Eigen::Matrix<double, 1, -1> yt = y.transpose();
      model_t model_i(t, yt, rate, CL, Q, V2, V3, ka);
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

  PMXTwoCptODE f2cpt;
  const std::vector<double> theta{CL, Q, V2, V3, ka};
  const PMXOdeIntegrator<PkBdf> integrator;
  using ode_model_t = PKODEModel<double, double, double, double, PMXTwoCptODE>;
  auto f2 = [&](std::vector<double>& amt_vec) {
    double t = t0;
    Eigen::Matrix<double, -1, 1> y = y0;
    for (int i = 0; i < 100; ++i) {
      Eigen::Matrix<double, 1, -1> yt = y.transpose();
      ode_model_t model_i(t, yt, rate, theta, f2cpt);
      double t_next = t + ii;
      Eigen::Matrix<double, -1, 1> ys = model_i.solve(t_next,integrator);
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
  std::vector<double> amt_vec{1000.0};

  {
    cmt = 1;
    Eigen::Matrix<double, -1, 1> y1 = f1(amt_vec);
    Eigen::Matrix<double, -1, 1> y2 = f2(amt_vec);
    torsten::test::test_val(y1, y2);    
  }

  {
    cmt = 2;
    Eigen::Matrix<double, -1, 1> y1 = f1(amt_vec);
    Eigen::Matrix<double, -1, 1> y2 = f2(amt_vec);
    torsten::test::test_val(y1, y2);    
  }

  {
    // this is an unlikely dosing event but we still need to test the validity
    // of the solution
    cmt = 3;
    Eigen::Matrix<double, -1, 1> y1 = f1(amt_vec);
    Eigen::Matrix<double, -1, 1> y2 = f2(amt_vec);
    torsten::test::test_val(y1, y2);    
  }
}

TEST_F(TorstenTwoCptModelTest, ss_infusion_grad_vs_long_run_sd) {
  y0[0] = 150;
  y0[1] = 50;
  using model_t = PMXTwoCptModel<var, var, var, double>;

  int cmt = 0;
  double ii = 6.0;
  double amt = 1000;
  
  auto f1 = [&](std::vector<var>& rate_vec) {
    var t = t0;
    Eigen::Matrix<var, -1, 1> y = y0;
    var t_infus = amt/rate_vec[cmt - 1];
    const std::vector<var> rate_zero{0.0, 0.0, 0.0};
    for (int i = 0; i < 100; ++i) {
      Eigen::Matrix<var, 1, -1> yt;
      Eigen::Matrix<var, -1, 1> ys;
      yt = y.transpose();
      model_t model_i(t, yt, rate_vec, CL, Q, V2, V3, ka);
      var t_next = t + t_infus;
      ys = model_i.solve(t_next);
      yt = ys.transpose();
      t = t_next;
      model_t model_j(t, yt, rate_zero, CL, Q, V2, V3, ka);
      t_next = t + ii - t_infus;
      ys = model_j.solve(t_next);
      y = ys;
    }
    return y;
  };

  auto f2 = [&](std::vector<var>& rate_vec) {
    PMXTwoCptModel<double, double, double, double> model(t0, y0, rate, CL, Q, V2, V3, ka);
    return model.solve(amt, rate_vec[cmt - 1], ii, cmt);
  };

  std::vector<var> rate_vec(3, 0.0);

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

  {
    cmt = 3;
    rate_vec[cmt-1] = 300.0;
    // For SS we only compare gradient wrt the dosing compartment as the
    // non-dosing compartment doesn't enter the SS system
    Eigen::Matrix<var, -1, 1> y1 = f1(rate_vec);
    Eigen::Matrix<var, -1, 1> y2 = f2(rate_vec);
    std::vector<var> par{rate_vec[cmt - 1]};
    torsten::test::test_grad(par, y1, y2, 5.e-11, 1.e-10);
    rate_vec[cmt-1] = 0.0;
  }
}

TEST_F(TorstenTwoCptModelTest, ss_infusion_by_long_run_sd_vs_bdf_result) {
  y0[0] = 150;
  y0[1] = 50;
  using model_t = refactor::PMXTwoCptModel<double, double, double, double>;

  int cmt = 0;
  const double ii = 11.9;
  const double amt = 1000;
  
  auto f1 = [&](std::vector<double>& rate_vec) {
    double t = t0;
    Eigen::Matrix<double, -1, 1> y = y0;
    double t_infus = amt/rate_vec[cmt - 1];
    const std::vector<double> rate_zero{0.0, 0.0, 0.};
    for (int i = 0; i < 50; ++i) {
      Eigen::Matrix<double, 1, -1> yt;
      Eigen::Matrix<double, -1, 1> ys;
      yt = y.transpose();
      model_t model_i(t, yt, rate_vec, CL, Q, V2, V3, ka);
      double t_next = t + t_infus;
      ys = model_i.solve(t_next);
      yt = ys.transpose();
      t = t_next;
      model_t model_j(t, yt, rate_zero, CL, Q, V2, V3, ka);
      t_next = t + ii - t_infus;
      ys = model_j.solve(t_next);
      y = ys;
    }
    return y;
  };

  PMXTwoCptODE f2cpt;
  const std::vector<double> theta{CL, Q, V2, V3, ka};
  const PMXOdeIntegrator<PkBdf> integrator;
  using ode_model_t = refactor::PKODEModel<double, double, double, double, PMXTwoCptODE>;
  auto f2 = [&](std::vector<double>& rate_vec) {
    double t = t0;
    Eigen::Matrix<double, -1, 1> y = y0;
    double t_infus = amt/rate_vec[cmt - 1];
    const std::vector<double> rate_zero{0.0, 0.0, 0.0};
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
  
  std::vector<double> rate_vec(3, 0.0);

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

  {
    cmt = 3;
    rate_vec[cmt - 1] = 300.0;
    Eigen::Matrix<double, -1, 1> y1 = f1(rate_vec);
    Eigen::Matrix<double, -1, 1> y2 = f2(rate_vec);
    torsten::test::test_val(y1, y2);
  }
}

TEST_F(TorstenTwoCptModelTest, ss_infusion_grad_by_long_run_sd_vs_bdf_result) {
  y0[0] = 150;
  y0[1] = 50;
  using model_t = refactor::PMXTwoCptModel<var, var, var, double>;

  int cmt = 0;
  const double ii = 11.9;
  const double amt = 1000;
  
  auto f1 = [&](std::vector<var>& rate_vec) {
    var t = t0;
    Eigen::Matrix<var, -1, 1> y = y0;
    var t_infus = amt/rate_vec[cmt - 1];
    const std::vector<var> rate_zero(3, 0.0);
    for (int i = 0; i < 50; ++i) {
      Eigen::Matrix<var, 1, -1> yt;
      Eigen::Matrix<var, -1, 1> ys;
      yt = y.transpose();
      model_t model_i(t, yt, rate_vec, CL, Q, V2, V3, ka);
      var t_next = t + t_infus;
      ys = model_i.solve(t_next);
      yt = ys.transpose();
      t = t_next;
      model_t model_j(t, yt, rate_zero, CL, Q, V2, V3, ka);
      t_next = t + ii - t_infus;
      ys = model_j.solve(t_next);
      y = ys;
    }
    return y;
  };

  PMXTwoCptODE f2cpt;
  const std::vector<double> theta{CL, Q, V2, V3, ka};
  const PMXOdeIntegrator<PkBdf> integrator;
  using ode_model_t = refactor::PKODEModel<var, var, var, double, PMXTwoCptODE>;
  auto f2 = [&](std::vector<var>& rate_vec) {
    var t = t0;
    Eigen::Matrix<var, -1, 1> y = y0;
    var t_infus = amt/rate_vec[cmt - 1];
    const std::vector<var> rate_zero(3, 0.0);
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
  
  std::vector<var> rate_vec(3, 0.0);

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
    torsten::test::test_grad(rate_vec, y1, y2, 5.e-7, 1.e-9);
  }

  {
    cmt = 3;
    rate_vec[cmt - 1] = 300.0;
    Eigen::Matrix<var, -1, 1> y1 = f1(rate_vec);
    Eigen::Matrix<var, -1, 1> y2 = f2(rate_vec);
    torsten::test::test_grad(rate_vec, y1, y2, 5.e-7, 5.e-8);
  }
}

TEST_F(TorstenTwoCptModelTest, ss_bolus_grad_by_long_run_sd_vs_bdf_result) {
  y0[0] = 150;
  y0[1] = 50;
  y0[1] = 50;
  rate[0] = 0;
  rate[1] = 0;
  rate[2] = 0;
  using model_t = refactor::PMXTwoCptModel<double, var, double, double>;

  int cmt = 0;
  const double ii = 11.9;
  
  auto f1 = [&](std::vector<var>& amt_vec) {
    double t = t0;
    Eigen::Matrix<var, -1, 1> y = y0;
    for (int i = 0; i < 50; ++i) {
      Eigen::Matrix<var, 1, -1> yt = y.transpose();
      model_t model_i(t, yt, rate, CL, Q, V2, V3, ka);
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

  PMXTwoCptODE f2cpt;
  const std::vector<double> theta{CL, Q, V2, V3, ka};
  const PMXOdeIntegrator<PkBdf> integrator;
  using ode_model_t = refactor::PKODEModel<double, var, double, double, PMXTwoCptODE>;
  auto f2 = [&](std::vector<var>& amt_vec) {
    double t = t0;
    Eigen::Matrix<var, -1, 1> y = y0;
    for (int i = 0; i < 50; ++i) {
      Eigen::Matrix<var, 1, -1> yt = y.transpose();
      ode_model_t model_i(t, yt, rate, theta, f2cpt);
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
    torsten::test::test_grad(amt_vec, y1, y2, 5.e-8, 5.e-10);
  }

  {
    cmt = 3;
    Eigen::Matrix<var, -1, 1> y1 = f1(amt_vec);
    Eigen::Matrix<var, -1, 1> y2 = f2(amt_vec);
    torsten::test::test_grad(amt_vec, y1, y2, 5.e-8, 5.e-10);
  }
}

TEST_F(TorstenTwoCptModelTest, ss_bolus_grad_vs_long_run_sd) {
  y0[0] = 150;
  y0[1] = 50;
  y0[2] = 50;
  using model_t = PMXTwoCptModel<double, var, double, double>;

  int cmt = 0;
  double ii = 11.9;
  
  auto f1 = [&](std::vector<var>& amt_vec) {
    double t = t0;
    Eigen::Matrix<var, -1, 1> y = y0;
    for (int i = 0; i < 100; ++i) {
      Eigen::Matrix<var, 1, -1> yt = y.transpose();
      model_t model_i(t, yt, rate, CL, Q, V2, V3, ka);
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

  auto f2 = [&](std::vector<var>& amt_vec) {
    PMXTwoCptModel<double, double, double, double> model(t0, y0, rate, CL, Q, V2, V3, ka);
    return model.solve(amt_vec[0], rate[cmt - 1], ii, cmt);
  };

  std::vector<var> amt_vec{300.00};

  {
    cmt = 1;
    Eigen::Matrix<var, -1, 1> y1 = f1(amt_vec);
    Eigen::Matrix<var, -1, 1> y2 = f2(amt_vec);
    torsten::test::test_grad(amt_vec, y1, y2, 1.e-11, 1.e-10);
  }

  {
    cmt = 2;
    Eigen::Matrix<var, -1, 1> y1 = f1(amt_vec);
    Eigen::Matrix<var, -1, 1> y2 = f2(amt_vec);
    torsten::test::test_grad(amt_vec, y1, y2, 5.e-12, 1.e-11);
  }

  {
    cmt = 3;
    Eigen::Matrix<var, -1, 1> y1 = f1(amt_vec);
    Eigen::Matrix<var, -1, 1> y2 = f2(amt_vec);
    torsten::test::test_grad(amt_vec, y1, y2, 1.e-11, 1.e-10);
  }
}

TEST_F(TorstenTwoCptModelTest, ss_const_infusion_grad_by_long_run_sd_vs_bdf_result) {
  y0[0] = 150;
  y0[1] = 50;
  using model_t = refactor::PMXTwoCptModel<double, double, var, double>;

  int cmt = 0;
  
  auto f1 = [&](std::vector<var>& rate_vec) {
    model_t model(t0, y0, rate_vec, CL, Q, V2, V3, ka);
    double t_next = 1.0e3;
    return model.solve(t_next);
  };

  PMXTwoCptODE f2cpt;
  const std::vector<double> theta{CL, Q, V2, V3, ka};
  const PMXOdeIntegrator<PkBdf> integrator;
  using ode_model_t = refactor::PKODEModel<double, double, var, double, PMXTwoCptODE>;
  auto f2 = [&](std::vector<var>& rate_vec) {
    ode_model_t model(t0, y0, rate_vec, theta, f2cpt);
    double t_next = 1.0e3;
    return model.solve(t_next, integrator);
  };
  
  std::vector<var> rate_vec(3, 0.0);

  {
    cmt = 1;
    rate_vec[cmt - 1] = 300.0;
    Eigen::Matrix<var, -1, 1> y1 = f1(rate_vec);
    Eigen::Matrix<var, -1, 1> y2 = f2(rate_vec);
    torsten::test::test_grad(rate_vec, y1, y2, 5.e-9, 1.e-10);
    rate_vec[cmt - 1] = 0.0;
  }

  {
    cmt = 2;
    rate_vec[cmt - 1] = 300.0;
    Eigen::Matrix<var, -1, 1> y1 = f1(rate_vec);
    Eigen::Matrix<var, -1, 1> y2 = f2(rate_vec);
    torsten::test::test_grad(rate_vec, y1, y2, 5.e-8, 1.e-10);
  }

  {
    cmt = 3;
    rate_vec[cmt - 1] = 300.0;
    Eigen::Matrix<var, -1, 1> y1 = f1(rate_vec);
    Eigen::Matrix<var, -1, 1> y2 = f2(rate_vec);
    torsten::test::test_grad(rate_vec, y1, y2, 5.e-8, 1.e-10);
  }
}

TEST_F(TorstenTwoCptModelTest, ss_const_infusion_grad_vs_long_run_sd) {
  y0[0] = 150;
  y0[1] = 50;
  using model_t = refactor::PMXTwoCptModel<double, double, var, double>;
  using ss_model_t = refactor::PMXTwoCptModel<double, double, double, double>;

  int cmt = 0;
  
  auto f1 = [&](std::vector<var>& rate_vec) {
    model_t model(t0, y0, rate_vec, CL, Q, V2, V3, ka);
    double t_next = 1.0e3;
    return model.solve(t_next);
  };

  auto f2 = [&](std::vector<var>& rate_vec) {
    ss_model_t model(t0, y0, rate, CL, Q, V2, V3, ka);
    double amt = 0.0;
    double t_next = 1.0e2;
    double ii = 0.0;
    return model.solve(amt, rate_vec[cmt - 1], ii, cmt);
  };

  std::vector<var> rate_vec(3, 0.0);

  {
    cmt = 1;
    rate_vec[cmt - 1] = 300.0;
    Eigen::Matrix<var, -1, 1> y1 = f1(rate_vec);
    Eigen::Matrix<var, -1, 1> y2 = f2(rate_vec);
    std::vector<var> par{rate_vec[cmt - 1]};
    torsten::test::test_grad(par, y1, y2, 1.e-12, 1.e-12);
    rate_vec[cmt - 1] = 0.0;
  }

  {
    cmt = 2;
    rate_vec[cmt - 1] = 300.0;
    Eigen::Matrix<var, -1, 1> y1 = f1(rate_vec);
    Eigen::Matrix<var, -1, 1> y2 = f2(rate_vec);
    std::vector<var> par{rate_vec[cmt - 1]};
    torsten::test::test_grad(par, y1, y2, 1.e-12, 1.e-12);
    rate_vec[cmt - 1] = 0.0;
  }

  {
    cmt = 3;
    rate_vec[cmt - 1] = 300.0;
    Eigen::Matrix<var, -1, 1> y1 = f1(rate_vec);
    Eigen::Matrix<var, -1, 1> y2 = f2(rate_vec);
    std::vector<var> par{rate_vec[cmt - 1]};
    torsten::test::test_grad(par, y1, y2, 1.e-12, 1.e-12);
    rate_vec[cmt - 1] = 0.0;
  }
}

TEST_F(TorstenTwoCptModelTest, ss_bolus) {
  rate[0] = 0;
  rate[1] = 0;
  rate[2] = 0;
  y0[0] = 150;
  y0[1] = 250;
  y0[2] = 350;
  ts[0] = 10.0;
  ts.resize(1);
  var CLv = to_var(CL);
  var Qv  = to_var(Q);
  var V2v = to_var(V2);
  var V3v = to_var(V3);
  var kav = to_var(ka);
  std::vector<stan::math::var> rate_var{to_var(rate)};
  using model_t = PMXTwoCptModel<double, double, var, var>;
  model_t model(t0, y0, rate_var, CLv, Qv, V2v, V3v, kav);
  std::vector<var> theta{model.par()};

  std::cout.precision(12);

  double amt = 1800;
  int cmt = 1;
  double ii = 12.5;
  
  auto y1 = model.solve(amt, rate_var[cmt - 1], ii, cmt);
  EXPECT_FLOAT_EQ(y1(0).val(), 0.00055062433);
  EXPECT_FLOAT_EQ(y1(1).val(), 738.870108248);
  EXPECT_FLOAT_EQ(y1(2).val(), 763.661542764);

  std::vector<double> g1, g2;
  stan::math::set_zero_all_adjoints();
  y1(0).grad(theta, g1);
  EXPECT_FLOAT_EQ(g1[0], 0.0);
  EXPECT_FLOAT_EQ(g1[1], 0.0);
  EXPECT_FLOAT_EQ(g1[2], 0.0);
  EXPECT_FLOAT_EQ(g1[3], 0.0);
  EXPECT_FLOAT_EQ(g1[4], -0.0068828063);
  stan::math::set_zero_all_adjoints();
  y1(1).grad(theta, g1);
  EXPECT_FLOAT_EQ(g1[0], -106.766204214);
  EXPECT_FLOAT_EQ(g1[1], 1.71716124388 );
  EXPECT_FLOAT_EQ(g1[2], 11.6448699701 );
  EXPECT_FLOAT_EQ(g1[3], 1.25702756718 );
  EXPECT_FLOAT_EQ(g1[4], -33.4223234266);
  stan::math::set_zero_all_adjoints();
  y1(2).grad(theta, g1);
  EXPECT_FLOAT_EQ(g1[0],  -97.6726925784);
  EXPECT_FLOAT_EQ(g1[1],  -2.69921403074);
  EXPECT_FLOAT_EQ(g1[2], 1.70090555664  );
  EXPECT_FLOAT_EQ(g1[3], 13.0890353445  );
  EXPECT_FLOAT_EQ(g1[4],  -34.1903160204);

  cmt = 2;
  y1 = model.solve(amt, rate_var[cmt - 1], ii, cmt);
  
  EXPECT_FLOAT_EQ(y1(0).val(), 0.0);
  EXPECT_FLOAT_EQ(y1(1).val(), 700.955088612);
  EXPECT_FLOAT_EQ(y1(2).val(), 724.611476748);

  stan::math::set_zero_all_adjoints();
  y1(0).grad(theta, g1);
  EXPECT_FLOAT_EQ(g1[0], 0.0);
  EXPECT_FLOAT_EQ(g1[1], 0.0);
  EXPECT_FLOAT_EQ(g1[2], 0.0);
  EXPECT_FLOAT_EQ(g1[3], 0.0);
  EXPECT_FLOAT_EQ(g1[4], 0.0);
  stan::math::set_zero_all_adjoints();
  y1(1).grad(theta, g1);
  EXPECT_FLOAT_EQ(g1[0], -104.758709533);
  EXPECT_FLOAT_EQ(g1[1], 1.53235878944 );
  EXPECT_FLOAT_EQ(g1[2], 11.2582785156 );
  EXPECT_FLOAT_EQ(g1[3], 1.48598239969 );
  EXPECT_FLOAT_EQ(g1[4], 0);
  stan::math::set_zero_all_adjoints();
  y1(2).grad(theta, g1);
  EXPECT_FLOAT_EQ(g1[0], -96.2552712812);
  EXPECT_FLOAT_EQ(g1[1], -2.69328253921);
  EXPECT_FLOAT_EQ(g1[2], 1.83656017604 );
  EXPECT_FLOAT_EQ(g1[3], 12.7291401404 );
  EXPECT_FLOAT_EQ(g1[4], 0);

  cmt = 3;
  y1 = model.solve(amt, rate_var[cmt - 1], ii, cmt);
  
  EXPECT_FLOAT_EQ(y1(0).val(), 0.0);
  EXPECT_FLOAT_EQ(y1(1).val(), 828.127401998);
  EXPECT_FLOAT_EQ(y1(2).val(), 856.228976486);

  stan::math::set_zero_all_adjoints();
  y1(0).grad(theta, g1);
  EXPECT_FLOAT_EQ(g1[0], 0.0);
  EXPECT_FLOAT_EQ(g1[1], 0.0);
  EXPECT_FLOAT_EQ(g1[2], 0.0);
  EXPECT_FLOAT_EQ(g1[3], 0.0);
  EXPECT_FLOAT_EQ(g1[4], 0.0);
  stan::math::set_zero_all_adjoints();
  y1(1).grad(theta, g1);
  EXPECT_FLOAT_EQ(g1[0], -110.006024321);
  EXPECT_FLOAT_EQ(g1[1], -3.07803718766);
  EXPECT_FLOAT_EQ(g1[2], 12.4505184404 );
  EXPECT_FLOAT_EQ(g1[3], 2.71719727476 );
  EXPECT_FLOAT_EQ(g1[4], 0.0);
  stan::math::set_zero_all_adjoints();
  y1(2).grad(theta, g1);
  EXPECT_FLOAT_EQ(g1[0], -99.5058577806);
  EXPECT_FLOAT_EQ(g1[1], -8.28726650911);
  EXPECT_FLOAT_EQ(g1[2], 1.30023459973 );
  EXPECT_FLOAT_EQ(g1[3], 16.044046744  );
  EXPECT_FLOAT_EQ(g1[4], 0.0);
}

TEST_F(TorstenTwoCptModelTest, ss_multi_trunc_infusion) {
  rate[0] = 1200;
  rate[1] = 1100;
  rate[2] = 800;
  y0[0] = 150;
  y0[1] = 250;
  y0[2] = 350;
  ts[0] = 10.0;
  ts.resize(1);
  var CLv = to_var(CL);
  var Qv  = to_var(Q);
  var V2v = to_var(V2);
  var V3v = to_var(V3);
  var kav = to_var(ka);
  std::vector<stan::math::var> rate_var{to_var(rate)};
  using model_t = PMXTwoCptModel<double, double, var, var>;
  model_t model(t0, y0, rate_var, CLv, Qv, V2v, V3v, kav);
  std::vector<var> theta(model.par());

  std::cout.precision(12);

  double amt = 1800;
  int cmt = 1;
  double ii = 12.5;
  std::vector<double> g1, g2;

  auto y1 = model.solve(amt, rate_var[cmt - 1], ii, cmt);
  EXPECT_FLOAT_EQ(y1(0).val(), 0.00154469934961);
  EXPECT_FLOAT_EQ(y1(1).val(), 774.104413167);
  EXPECT_FLOAT_EQ(y1(2).val(), 799.879747792);

  stan::math::set_zero_all_adjoints();
  y1(0).grad(theta, g1);
  EXPECT_FLOAT_EQ(g1[0], 0.0);
  EXPECT_FLOAT_EQ(g1[1], 0.0);
  EXPECT_FLOAT_EQ(g1[2], 0.0);
  EXPECT_FLOAT_EQ(g1[3], 0.0);
  EXPECT_FLOAT_EQ(g1[4], -0.0178200945892);
  stan::math::set_zero_all_adjoints();
  y1(1).grad(theta, g1);
  EXPECT_FLOAT_EQ(g1[0], -108.530709442);
  EXPECT_FLOAT_EQ(g1[1], 1.88466941803 );
  EXPECT_FLOAT_EQ(g1[2], 11.9988317957 );
  EXPECT_FLOAT_EQ(g1[3], 1.0375686723  );
  EXPECT_FLOAT_EQ(g1[4], -35.1732711317);
  stan::math::set_zero_all_adjoints();
  y1(2).grad(theta, g1);
  EXPECT_FLOAT_EQ(g1[0], -98.8831379736);
  EXPECT_FLOAT_EQ(g1[1], -2.69526994326);
  EXPECT_FLOAT_EQ(g1[2], 1.56750688457 );
  EXPECT_FLOAT_EQ(g1[3], 13.4128341054 );
  EXPECT_FLOAT_EQ(g1[4], -35.6693548262);

  cmt = 2;
  y1 = model.solve(amt, rate_var[cmt - 1], ii, cmt);
  
  EXPECT_FLOAT_EQ(y1(0).val(), 0.0);
  EXPECT_FLOAT_EQ(y1(1).val(), 737.454228957);
  EXPECT_FLOAT_EQ(y1(2).val(), 762.268194176);

  stan::math::set_zero_all_adjoints();
  y1(0).grad(theta, g1);
  EXPECT_FLOAT_EQ(g1[0], 0.0);
  EXPECT_FLOAT_EQ(g1[1], 0.0);
  EXPECT_FLOAT_EQ(g1[2], 0.0);
  EXPECT_FLOAT_EQ(g1[3], 0.0);
  EXPECT_FLOAT_EQ(g1[4], 0.0);
  stan::math::set_zero_all_adjoints();
  y1(1).grad(theta, g1);
  EXPECT_FLOAT_EQ(g1[0], -106.758381233);
  EXPECT_FLOAT_EQ(g1[1], 1.71458102804 );
  EXPECT_FLOAT_EQ(g1[2], 11.6337986397 );
  EXPECT_FLOAT_EQ(g1[3], 1.26959503378 );
  EXPECT_FLOAT_EQ(g1[4], 0);
  stan::math::set_zero_all_adjoints();
  y1(2).grad(theta, g1);
  EXPECT_FLOAT_EQ(g1[0], -97.6910811146);
  EXPECT_FLOAT_EQ(g1[1], -2.70645214819);
  EXPECT_FLOAT_EQ(g1[2], 1.71099901438 );
  EXPECT_FLOAT_EQ(g1[3], 13.0830221449 );
  EXPECT_FLOAT_EQ(g1[4], 0.0);

  cmt = 3;
  y1 = model.solve(amt, rate_var[cmt - 1], ii, cmt);
  
  EXPECT_FLOAT_EQ(y1(0).val(), 0.0);
  EXPECT_FLOAT_EQ(y1(1).val(), 888.053512871);
  EXPECT_FLOAT_EQ(y1(2).val(), 918.312019204);

  stan::math::set_zero_all_adjoints();
  y1(0).grad(theta, g1);
  EXPECT_FLOAT_EQ(g1[0], 0.0);
  EXPECT_FLOAT_EQ(g1[1], 0.0);
  EXPECT_FLOAT_EQ(g1[2], 0.0);
  EXPECT_FLOAT_EQ(g1[3], 0.0);
  EXPECT_FLOAT_EQ(g1[4], 0.0);
  stan::math::set_zero_all_adjoints();
  y1(1).grad(theta, g1);
  EXPECT_FLOAT_EQ(g1[0], -112.218811747);
  EXPECT_FLOAT_EQ(g1[1], -3.09372149524);
  EXPECT_FLOAT_EQ(g1[2], 12.9947786493 );
  EXPECT_FLOAT_EQ(g1[3], 2.41757181997 );
  EXPECT_FLOAT_EQ(g1[4], 0.0);
  stan::math::set_zero_all_adjoints();
  y1(2).grad(theta, g1);
  EXPECT_FLOAT_EQ(g1[0], -100.772892327);
  EXPECT_FLOAT_EQ(g1[1], -8.7063180167 );
  EXPECT_FLOAT_EQ(g1[2], 1.03257281914 );
  EXPECT_FLOAT_EQ(g1[3], 16.69857146   );
  EXPECT_FLOAT_EQ(g1[4], 0.0);
}

TEST_F(TorstenTwoCptModelTest, ss_solver_const_infusion) {
  using stan::math::var;
  using stan::math::to_var;
  using refactor::PMXTwoCptModel;
  using torsten::pmx_integrate_ode_bdf;
  using stan::math::integrate_ode_bdf;
  using refactor::PMXTwoCptODE;
  using refactor::PMXOdeFunctorRateAdaptor;

  rate[0] = 1200;
  rate[1] = 1100;
  rate[2] = 800;
  y0[0] = 150;
  y0[1] = 250;
  y0[2] = 350;
  ts[0] = 10.0;
  ts.resize(1);
  var CLv = to_var(CL);
  var Qv  = to_var(Q);
  var V2v = to_var(V2);
  var V3v = to_var(V3);
  var kav = to_var(ka);
  std::vector<stan::math::var> rate_var{to_var(rate)};
  using model_t = PMXTwoCptModel<double, double, var, var>;
  model_t model(t0, y0, rate_var, CLv, Qv, V2v, V3v, kav);
  std::vector<var> theta(model.par());

  std::cout.precision(12);

  double amt = 1800;
  int cmt = 1;
  double ii = 0.0;
  std::vector<double> g1, g2;

  auto y1 = model.solve(amt, rate_var[cmt - 1], ii, cmt);
  EXPECT_FLOAT_EQ(y1(0).val(), 1000);
  EXPECT_FLOAT_EQ(y1(1).val(), 9600);
  EXPECT_FLOAT_EQ(y1(2).val(), 8400);

  stan::math::set_zero_all_adjoints();
  y1(0).grad(theta, g1);
  EXPECT_FLOAT_EQ(g1[0], 0.0);
  EXPECT_FLOAT_EQ(g1[1], 0.0);
  EXPECT_FLOAT_EQ(g1[2], 0.0);
  EXPECT_FLOAT_EQ(g1[3], 0.0);
  EXPECT_FLOAT_EQ(g1[4], -833.333333333);
  stan::math::set_zero_all_adjoints();
  y1(1).grad(theta, g1);
  EXPECT_FLOAT_EQ(g1[0], -960);
  EXPECT_FLOAT_EQ(g1[1], 2.27373675443e-14 );
  EXPECT_FLOAT_EQ(g1[2], 120 );
  EXPECT_FLOAT_EQ(g1[3], -1.81898940355e-14 );
  EXPECT_FLOAT_EQ(g1[4], 0.0);
  stan::math::set_zero_all_adjoints();
  y1(2).grad(theta, g1);
  EXPECT_FLOAT_EQ(g1[0], -840);
  EXPECT_FLOAT_EQ(g1[1], 1.70530256582e-13);
  EXPECT_FLOAT_EQ(g1[2], -7.1054273576e-14);
  EXPECT_FLOAT_EQ(g1[3], 120 );
  // EXPECT_FLOAT_EQ(g1[4], 1.02318153949e-12); //FIXME: fail on g++ but not clang
  EXPECT_NEAR(g1[4], 1.02318153949e-12, 2.e-12);

  // FIXME: check the results
  // cmt = 2;
  // y1 = model.solve(amt, rate_var[cmt - 1], ii, cmt);

  // EXPECT_FLOAT_EQ(y1(0).val(), 1000.0);
  // EXPECT_FLOAT_EQ(y1(1).val(), 9600);
  // EXPECT_FLOAT_EQ(y1(2).val(), 8400);

  // stan::math::set_zero_all_adjoints();
  // y1(0).grad(theta, g1);
  // EXPECT_FLOAT_EQ(g1[0], 0.0);
  // EXPECT_FLOAT_EQ(g1[1], 0.0);
  // EXPECT_FLOAT_EQ(g1[2], 0.0);
  // EXPECT_FLOAT_EQ(g1[3], 0.0);
  // EXPECT_FLOAT_EQ(g1[4], -833.333333333);
  // stan::math::set_zero_all_adjoints();
  // y1(1).grad(theta, g1);
  // EXPECT_FLOAT_EQ(g1[0], -960);
  // EXPECT_FLOAT_EQ(g1[1], 2.27373675443e-14);
  // EXPECT_FLOAT_EQ(g1[2], 120);
  // EXPECT_FLOAT_EQ(g1[3], -1.81898940355e-14);
  // EXPECT_FLOAT_EQ(g1[4], 0);
  // stan::math::set_zero_all_adjoints();
  // y1(2).grad(theta, g1);
  // EXPECT_FLOAT_EQ(g1[0], -840);
  // EXPECT_FLOAT_EQ(g1[1], 1.70530256582e-13);
  // EXPECT_FLOAT_EQ(g1[2], -7.1054273576e-14);
  // EXPECT_FLOAT_EQ(g1[3], 120);
  // EXPECT_FLOAT_EQ(g1[4], 1.02318153949e-12);

  cmt = 3;
  y1 = model.solve(amt, rate_var[cmt - 1], ii, cmt);
  
  EXPECT_FLOAT_EQ(y1(0).val(), 0.0);
  EXPECT_FLOAT_EQ(y1(1).val(), 6400);
  EXPECT_FLOAT_EQ(y1(2).val(), 7600);

  stan::math::set_zero_all_adjoints();
  y1(0).grad(theta, g1);
  EXPECT_FLOAT_EQ(g1[0], 0.0);
  EXPECT_FLOAT_EQ(g1[1], 0.0);
  EXPECT_FLOAT_EQ(g1[2], 0.0);
  EXPECT_FLOAT_EQ(g1[3], 0.0);
  EXPECT_FLOAT_EQ(g1[4], 0.0);
  stan::math::set_zero_all_adjoints();
  y1(1).grad(theta, g1);
  EXPECT_FLOAT_EQ(g1[0], -640);
  EXPECT_FLOAT_EQ(g1[1], -1.18559130767e-13);
  EXPECT_FLOAT_EQ(g1[2], 80);
  EXPECT_FLOAT_EQ(g1[3], 1.55913377447e-14);
  EXPECT_FLOAT_EQ(g1[4], 0);
  stan::math::set_zero_all_adjoints();
  y1(2).grad(theta, g1);
  EXPECT_FLOAT_EQ(g1[0], -560);
  EXPECT_FLOAT_EQ(g1[1], -71.4285714286);
  EXPECT_FLOAT_EQ(g1[2], 5.68434188608e-14);
  EXPECT_FLOAT_EQ(g1[3], 108.571428571);
  EXPECT_FLOAT_EQ(g1[4], 0.0);
}
