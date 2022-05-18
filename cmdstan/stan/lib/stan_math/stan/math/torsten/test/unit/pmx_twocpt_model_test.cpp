#include <stan/math.hpp>
#include <stan/math/rev/core.hpp>
#include <stan/math/torsten/dsolve/pmx_ode_integrator.hpp>
#include <test/unit/math/rev/fun/util.hpp>
#include <stan/math/torsten/test/unit/pmx_cpt_model_test_fixture.hpp>
#include <test/unit/util.hpp>
#include <gtest/gtest.h>

using stan::math::var;
using stan::math::to_var;
using stan::math::vector_v;
using stan::math::matrix_v;
using torsten::PMXTwoCptModel;
using torsten::pmx_integrate_ode_bdf;
using stan::math::integrate_ode_bdf;
using torsten::PMXTwoCptODE;
using torsten::PMXOdeFunctorRateAdaptor;
using torsten::dsolve::PMXOdeIntegrator;
using torsten::dsolve::PMXVariadicOdeSystem;
using torsten::dsolve::PMXCvodesIntegrator;
using torsten::dsolve::PMXOdeintIntegrator;

PMXTwoCptODE f0;

TEST_F(TorstenTwoCptModelTest, ka_zero) {
  y0(0) = 745;
  y0(1) = 100;
  y0(2) = 130;  
  rate[0] = 1200;
  rate[1] = 200;
  rate[2] = 300;
  ka = 0.0;
  using model_t = PMXTwoCptModel<double>;
  model_t model(CL, Q, V2, V3, ka);
  torsten::PKRec<double> y(y0);
  model.solve(y, t0, ts[0], rate);
  EXPECT_FLOAT_EQ(y(0), 865.0);
  EXPECT_FLOAT_EQ(y(1), 120.52635);
  EXPECT_FLOAT_EQ(y(2), 158.09552);
}

TEST_F(TorstenTwoCptModelTest, rate_dbl) {
  rate[0] = 1200;
  rate[1] = 200;
  rate[2] = 300;
  using model_t = PMXTwoCptModel<double>;
  model_t model(CL, Q, V2, V3, ka);
  PMXOdeFunctorRateAdaptor<PMXTwoCptODE> f1(f0);

  Eigen::VectorXd y = f1(t0, y0, msgs, model.par(),rate, x_r, x_i);
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
  using model_t = PMXTwoCptModel<var>;
  model_t model(CLv, Qv, V2v, V3v, kav);
  std::vector<stan::math::var> theta(model.par());
  PMXOdeFunctorRateAdaptor<PMXTwoCptODE> f1(f0);

  Eigen::Matrix<var, -1, 1> y = f1(t0, y0, msgs, theta, rate_var, x_r, x_i);
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
  // std::vector<var> theta{CLv, Qv, V2v, V3v, kav};
  std::vector<stan::math::var> rate_var{to_var(rate)};
  using model_t = PMXTwoCptModel<var>;
  model_t model(CLv, Qv, V2v, V3v, kav);
  PMXOdeFunctorRateAdaptor<PMXTwoCptODE> f1(f0);
  std::vector<var> theta = model.par();

  auto y1 = pmx_ode_bdf(f1, y0, t0, ts, msgs, theta, rate_var, x_r, x_i);
  torsten::PKRec<var> y2(to_var(y0));
  model.solve(y2, t0, ts[0], rate_var);

  torsten::test::test_grad(theta, y1[0], y2, 1.e-6, 1.e-6);
  torsten::test::test_grad(rate_var, y1[0], y2, 1.e-6, 1.e-7);
}

TEST_F(TorstenTwoCptModelTest, infusion_theta_grad) {
  rate[0] = 1100;
  rate[1] = 810;
  rate[2] = 200;
  y0[0] = 150;
  y0[1] = 50.0;
  y0[2] = 50.0;

  double dt = 2.5;
  
  auto f1 = [&](const std::vector<double>& pars) {
    using model_t = PMXTwoCptModel<double>;
    model_t model(pars[0], pars[1], pars[2], pars[3], pars[4]);
    torsten::PKRec<double> y(y0);
    model.solve(y, t0, dt, rate);
    return y;
  };
  auto f2 = [&](const std::vector<double>& pars) {
    using model_t = PMXTwoCptModel<double>;
    model_t model(pars[0], pars[1], pars[2], pars[3], pars[4]);
    torsten::PKRec<double> y(y0);
    model.solve_analytical(y, t0, dt, rate);
    return y;
  };
  auto f3 = [&](const std::vector<var>& pars) {
    using model_t = PMXTwoCptModel<var>;
    torsten::PKRec<var> y(to_var(y0));
    model_t model(pars[0], pars[1], pars[2], pars[3], pars[4]);
    model.solve(y, t0, dt, rate);
    return y;
  };

  std::vector<double> pars{CL, Q, V2, V3, ka};
  torsten::test::test_grad(f1, f3, pars, 1.e-3, 2.e-12, 1.e-6, 2.e-10);
  torsten::test::test_grad(f2, f3, pars, 1.e-3, 3.e-12, 1.e-6, 2.e-10);
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
  using model_t = PMXTwoCptModel<double>;
  model_t model(CL, Q, V2, V3, ka);

  double ii = 12.5;
  
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
  cmt = 3; torsten::test::test_grad(f1, f2, amt_vec, 1.e-3, 1.e-16, 2.e-10, 1.e-12);

  // compare against textbook analytical sol.
  auto f3 = [&](const std::vector<var>& amt_vec) {
    return model.solve_analytical(t0, amt_vec[0], rate[cmt-1], ii, cmt);
  };
  torsten::PKRec<var> y2, y3;
  std::vector<var> amt_vec_var = stan::math::to_var(amt_vec);
  for (int i = 0; i < 3; ++i) {
    cmt = i + 1;
    y2 = f2(amt_vec_var);
    y3 = f3(amt_vec_var);
    torsten::test::test_grad(amt_vec_var, y2, y3, 1.e-12, 5.e-16);
  }
}

TEST_F(TorstenTwoCptModelTest, ss_infusion_rate_grad) {
  y0[0] = 150;
  y0[1] = 250;
  y0[2] = 350;
  ts[0] = 10.0;
  ts.resize(1);
  using model_t = PMXTwoCptModel<double>;
  model_t model(CL, Q, V2, V3, ka);

  double amt = 1100;
  double ii = 12.5;
  
  int cmt = 0;
  auto f1 = [&](const std::vector<double>& rate_vec) {
    return model.solve(t0, amt, rate_vec[0], ii, cmt);
  };
  auto f2 = [&](const std::vector<var>& rate_vec) {
    return model.solve(t0, amt, rate_vec[0], ii, cmt);
  };
  std::vector<double> rate_vec{130.0};
  cmt = 1; torsten::test::test_grad(f1, f2, rate_vec, 1.e-3, 1.e-16, 1.e-8, 1.e-12);
  cmt = 2; torsten::test::test_grad(f1, f2, rate_vec, 1.e-3, 1.e-16, 1.e-9, 1.e-12);
  cmt = 3; torsten::test::test_grad(f1, f2, rate_vec, 1.e-3, 1.e-16, 1.e-9, 1.e-12);

  // compare against textbook analytical sol.
  auto f3 = [&](const std::vector<var>& rate_vec) {
    return model.solve_analytical(t0, amt, rate_vec[0], ii, cmt);
  };
  torsten::PKRec<var> y2, y3;
  std::vector<var> rate_var = stan::math::to_var(rate_vec);
  for (int i = 0; i < 3; ++i) {
    cmt = i + 1;
    y2 = f2(rate_var);
    y3 = f3(rate_var);
    torsten::test::test_grad(rate_var, y2, y3, 1.e-12, 5.e-15);
  }
}

TEST_F(TorstenTwoCptModelTest, ss_bolus_by_long_run_sd_vs_bdf_result) {
  using torsten::PKODEModel;

  rate[0] = 0;
  rate[1] = 0;
  rate[2] = 0;
  y0[0] = 150;
  y0[1] = 55;
  y0[2] = 120;
  using model_t = PMXTwoCptModel<double>;
  model_t model(CL, Q, V2, V3, ka);

  int cmt = 0;
  const double ii = 8.5;
  
  auto f1 = [&](const std::vector<double>& amt_vec) {
    double t = t0;
    Eigen::Matrix<double, -1, 1> y = y0;
    for (int i = 0; i < 100; ++i) {
      model_t model_i(CL, Q, V2, V3, ka);
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

  PMXTwoCptODE f2cpt;
  const std::vector<double> theta{CL, Q, V2, V3, ka};
  const PMXOdeIntegrator<PMXVariadicOdeSystem, PMXCvodesIntegrator<CV_BDF, CV_STAGGERED>> integrator;
  using ode_model_t = PKODEModel<double, PMXTwoCptODE>;
  auto f2 = [&](const std::vector<double>& amt_vec) {
    double t = t0;
    Eigen::Matrix<double, -1, 1> y = y0;
    for (int i = 0; i < 100; ++i) {
      ode_model_t model_i(theta, y.size(), f2cpt);
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
  std::vector<double> amt_vec{1000.0};

  {
    cmt = 1;
    Eigen::Matrix<double, -1, 1> y1 = f1(amt_vec);
    Eigen::Matrix<double, -1, 1> y2 = f2(amt_vec);
    EXPECT_VEC_VAL_FLOAT_EQ(y1, y2);
  }

  {
    cmt = 2;
    Eigen::Matrix<double, -1, 1> y1 = f1(amt_vec);
    Eigen::Matrix<double, -1, 1> y2 = f2(amt_vec);
    EXPECT_VEC_VAL_FLOAT_EQ(y1, y2);    
  }

  {
    // this is an unlikely dosing event but we still need to test the validity
    // of the solution
    cmt = 3;
    Eigen::Matrix<double, -1, 1> y1 = f1(amt_vec);
    Eigen::Matrix<double, -1, 1> y2 = f2(amt_vec);
    EXPECT_VEC_VAL_FLOAT_EQ(y1, y2);    
  }
}

TEST_F(TorstenTwoCptModelTest, ss_infusion_grad_vs_long_run_sd) {
  y0[0] = 150;
  y0[1] = 50;
  using model_t = PMXTwoCptModel<double>;

  int cmt = 0;
  double ii = 6.0;
  double amt = 1000;
  
  auto f1 = [&](const std::vector<var>& rate_vec) {
    var t = t0;
    Eigen::Matrix<var, -1, 1> y = y0;
    var t_infus = amt/rate_vec[cmt - 1];
    const std::vector<var> rate_zero{0.0, 0.0, 0.0};
    for (int i = 0; i < 100; ++i) {
      model_t model_i(CL, Q, V2, V3, ka);
      var t_next = t + t_infus;
      model_i.solve(y, t, t_next, rate_vec);
      t = t_next;
      model_t model_j(CL, Q, V2, V3, ka);
      t_next = t + ii - t_infus;
      model_j.solve(y, t, t_next, rate_zero);
    }
    return y;
  };

  auto f2 = [&](const std::vector<var>& rate_vec) {
    PMXTwoCptModel<double> model(CL, Q, V2, V3, ka);
    return model.solve(t0, amt, rate_vec[cmt - 1], ii, cmt);
  };

  for (int i = 0; i < 3; ++i) {
    cmt = i + 1;
    std::vector<var> rate_vec(3, 0.0);
    rate_vec[i] = 300.0;
    // For SS we only compare gradient wrt the dosing compartment as the
    // non-dosing compartment doesn't enter the SS system
    Eigen::Matrix<var, -1, 1> y1 = f1(rate_vec);
    Eigen::Matrix<var, -1, 1> y2 = f2(rate_vec);
    std::vector<var> par{rate_vec[i]};
    torsten::test::test_grad(par, y1, y2, 5.e-11, 1.e-10);
  }
}

TEST_F(TorstenTwoCptModelTest, ss_infusion_by_long_run_sd_vs_bdf_result) {
  y0[0] = 150;
  y0[1] = 50;
  using model_t = torsten::PMXTwoCptModel<double>;

  int cmt = 0;
  const double ii = 11.9;
  const double amt = 1000;
  
  auto f1 = [&](const std::vector<double>& rate_vec) {
    double t = t0;
    Eigen::Matrix<double, -1, 1> y = y0;
    double t_infus = amt/rate_vec[cmt - 1];
    const std::vector<double> rate_zero{0.0, 0.0, 0.};
    for (int i = 0; i < 50; ++i) {
      model_t model_i(CL, Q, V2, V3, ka);
      double t_next = t + t_infus;
      model_i.solve(y, t, t_next, rate_vec);
      t = t_next;
      model_t model_j(CL, Q, V2, V3, ka);
      t_next = t + ii - t_infus;
      model_j.solve(y, t, t_next, rate_zero);
    }
    return y;
  };

  PMXTwoCptODE f2cpt;
  const std::vector<double> theta{CL, Q, V2, V3, ka};
  const PMXOdeIntegrator<PMXVariadicOdeSystem, PMXCvodesIntegrator<CV_BDF, CV_STAGGERED>> integrator;
  using ode_model_t = torsten::PKODEModel<double, PMXTwoCptODE>;
  auto f2 = [&](const std::vector<double>& rate_vec) {
    double t = t0;
    Eigen::Matrix<double, -1, 1> y = y0;
    double t_infus = amt/rate_vec[cmt - 1];
    const std::vector<double> rate_zero{0.0, 0.0, 0.0};
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
  
  std::vector<double> rate_vec(3, 0.0);

  {
    cmt = 1;
    rate_vec[cmt - 1] = 300.0;
    Eigen::Matrix<double, -1, 1> y1 = f1(rate_vec);
    Eigen::Matrix<double, -1, 1> y2 = f2(rate_vec);
    EXPECT_VEC_VAL_FLOAT_EQ(y1, y2);
    rate_vec[cmt - 1] = 0.0;
  }

  {
    cmt = 2;
    rate_vec[cmt - 1] = 300.0;
    Eigen::Matrix<double, -1, 1> y1 = f1(rate_vec);
    Eigen::Matrix<double, -1, 1> y2 = f2(rate_vec);
    EXPECT_VEC_VAL_FLOAT_EQ(y1, y2);
  }

  {
    cmt = 3;
    rate_vec[cmt - 1] = 300.0;
    Eigen::Matrix<double, -1, 1> y1 = f1(rate_vec);
    Eigen::Matrix<double, -1, 1> y2 = f2(rate_vec);
    EXPECT_VEC_VAL_FLOAT_EQ(y1, y2);
  }
}

TEST_F(TorstenTwoCptModelTest, ss_infusion_grad_by_long_run_sd_vs_bdf_result) {
  y0[0] = 150;
  y0[1] = 50;
  using model_t = torsten::PMXTwoCptModel<double>;

  int cmt = 0;
  const double ii = 11.9;
  const double amt = 1000;
  
  auto f1 = [&](const std::vector<var>& rate_vec) {
    var t = t0;
    Eigen::Matrix<var, -1, 1> y = y0;
    var t_infus = amt/rate_vec[cmt - 1];
    const std::vector<var> rate_zero(3, 0.0);
    for (int i = 0; i < 50; ++i) {
      model_t model_i(CL, Q, V2, V3, ka);
      var t_next = t + t_infus;
      model_i.solve(y, t, t_next, rate_vec);
      t = t_next;
      model_t model_j(CL, Q, V2, V3, ka);
      t_next = t + ii - t_infus;
      model_j.solve(y, t, t_next, rate_zero);
    }
    return y;
  };

  PMXTwoCptODE f2cpt;
  const std::vector<double> theta{CL, Q, V2, V3, ka};
  const PMXOdeIntegrator<PMXVariadicOdeSystem, PMXCvodesIntegrator<CV_BDF, CV_STAGGERED>> integrator;
  using ode_model_t = torsten::PKODEModel<double, PMXTwoCptODE>;
  auto f2 = [&](const std::vector<var>& rate_vec) {
    var t = t0;
    Eigen::Matrix<var, -1, 1> y = y0;
    var t_infus = amt/rate_vec[cmt - 1];
    const std::vector<var> rate_zero(3, 0.0);
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
  using model_t = torsten::PMXTwoCptModel<double>;

  int cmt = 0;
  const double ii = 11.9;
  
  auto f1 = [&](const std::vector<var>& amt_vec) {
    double t = t0;
    Eigen::Matrix<var, -1, 1> y = y0;
    for (int i = 0; i < 50; ++i) {
      model_t model_i(CL, Q, V2, V3, ka);
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

  PMXTwoCptODE f2cpt;
  const std::vector<double> theta{CL, Q, V2, V3, ka};
  const PMXOdeIntegrator<PMXVariadicOdeSystem, PMXCvodesIntegrator<CV_BDF, CV_STAGGERED>> integrator;
  using ode_model_t = torsten::PKODEModel<double, PMXTwoCptODE>;
  auto f2 = [&](const std::vector<var>& amt_vec) {
    double t = t0;
    Eigen::Matrix<var, -1, 1> y = y0;
    for (int i = 0; i < 50; ++i) {
      ode_model_t model_i(theta, y.size(), f2cpt);
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
  using model_t = PMXTwoCptModel<double>;

  int cmt = 0;
  double ii = 11.9;
  
  auto f1 = [&](const std::vector<var>& amt_vec) {
    double t = t0;
    Eigen::Matrix<var, -1, 1> y = y0;
    for (int i = 0; i < 100; ++i) {
      // Eigen::Matrix<var, 1, -1> yt = y.transpose();
      // model_t model_i(t, yt, rate, CL, Q, V2, V3, ka);
      // double t_next = t + ii;
      // Eigen::Matrix<var, -1, 1> ys = model_i.solve(t_next);
      // ys(cmt - 1) += amt_vec[0];
      // y = ys;
      // t = t_next;
      model_t model_i(CL, Q, V2, V3, ka);
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
    PMXTwoCptModel<double> model(CL, Q, V2, V3, ka);
    return model.solve(t0, amt_vec[0], rate[cmt - 1], ii, cmt);
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
  using model_t = torsten::PMXTwoCptModel<double>;

  int cmt = 0;
  
  auto f1 = [&](std::vector<var>& rate_vec) {
    model_t model(CL, Q, V2, V3, ka);
    double t_next = 1.0e3;
    torsten::PKRec<var> y(to_var(y0));
    model.solve(y, t0, t_next, rate_vec);
    return y;
  };

  PMXTwoCptODE f2cpt;
  const std::vector<double> theta{CL, Q, V2, V3, ka};
  const PMXOdeIntegrator<PMXVariadicOdeSystem, PMXCvodesIntegrator<CV_BDF, CV_STAGGERED>> integrator;
  using ode_model_t = torsten::PKODEModel<double, PMXTwoCptODE>;
  auto f2 = [&](std::vector<var>& rate_vec) {
    ode_model_t model(theta, y0.size(), f2cpt);
    double t_next = 1.0e3;
    torsten::PKRec<var> y(to_var(y0));
    model.solve(y, t0, t_next, rate_vec, integrator);
    return y;
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
  using model_t = torsten::PMXTwoCptModel<double>;
  using ss_model_t = torsten::PMXTwoCptModel<double>;

  int cmt = 0;
  
  auto f1 = [&](std::vector<var>& rate_vec) {
    model_t model(CL, Q, V2, V3, ka);
    double t_next = 1.0e3;
    torsten::PKRec<var> y(to_var(y0));
    model.solve(y, t0, t_next, rate_vec);
    return y;
  };

  auto f2 = [&](std::vector<var>& rate_vec) {
    ss_model_t model(CL, Q, V2, V3, ka);
    double amt = 0.0;
    double t_next = 1.0e2;
    double ii = 0.0;
    return model.solve(t0, amt, rate_vec[cmt - 1], ii, cmt);
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
  using model_t = PMXTwoCptModel<var>;
  model_t model(CLv, Qv, V2v, V3v, kav);
  std::vector<var> theta{model.par()};

  std::cout.precision(12);

  double amt = 1800;
  int cmt = 1;
  double ii = 12.5;
  
  auto y1 = model.solve(t0, amt, rate_var[cmt - 1], ii, cmt);
  EXPECT_FLOAT_EQ(y1(0).val(), 0.00055062433);
  EXPECT_FLOAT_EQ(y1(1).val(), 738.870108248);
  EXPECT_FLOAT_EQ(y1(2).val(), 763.661542764);

  std::vector<double> g1, g2;
  stan::math::set_zero_all_adjoints();
  y1(0).grad(theta, g1);
  EXPECT_NEAR(g1[0], 0.0, 1.e-18);
  EXPECT_NEAR(g1[1], 0.0, 1.e-18);
  EXPECT_NEAR(g1[2], 0.0, 1.e-18);
  EXPECT_NEAR(g1[3], 0.0, 1.e-18);
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
  y1 = model.solve(t0, amt, rate_var[cmt - 1], ii, cmt);
  
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
  y1 = model.solve(t0, amt, rate_var[cmt - 1], ii, cmt);
  
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
  using model_t = PMXTwoCptModel<var>;
  model_t model(CLv, Qv, V2v, V3v, kav);
  std::vector<var> theta(model.par());

  std::cout.precision(12);

  double amt = 1800;
  int cmt = 1;
  double ii = 12.5;
  std::vector<double> g1, g2;

  auto y1 = model.solve(t0, amt, rate_var[cmt - 1], ii, cmt);
  EXPECT_FLOAT_EQ(y1(0).val(), 0.00154469934961);
  EXPECT_FLOAT_EQ(y1(1).val(), 774.104413167);
  EXPECT_FLOAT_EQ(y1(2).val(), 799.879747792);

  stan::math::set_zero_all_adjoints();
  y1(0).grad(theta, g1);
  EXPECT_NEAR(g1[0], 0.0, 1.e-18);
  EXPECT_NEAR(g1[1], 0.0, 1.e-18);
  EXPECT_NEAR(g1[2], 0.0, 1.e-18);
  EXPECT_NEAR(g1[3], 0.0, 1.e-18);
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
  y1 = model.solve(t0, amt, rate_var[cmt - 1], ii, cmt);
  
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
  y1 = model.solve(t0, amt, rate_var[cmt - 1], ii, cmt);
  
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
  using torsten::PMXTwoCptModel;
  using torsten::pmx_integrate_ode_bdf;
  using stan::math::integrate_ode_bdf;
  using torsten::PMXTwoCptODE;

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
  using model_t = PMXTwoCptModel<var>;
  model_t model(CLv, Qv, V2v, V3v, kav);
  std::vector<var> theta(model.par());

  std::cout.precision(12);

  double amt = 1800;
  int cmt = 1;
  double ii = 0.0;
  std::vector<double> g1, g2;

  auto y1 = model.solve(t0, amt, rate_var[cmt - 1], ii, cmt);
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
  EXPECT_NEAR(g1[1], 0.0, 1.e-12 );
  EXPECT_FLOAT_EQ(g1[2], 120 );
  EXPECT_NEAR(g1[3], 0.0, 1.e-12 );
  EXPECT_NEAR(g1[4], 0.0, 5.e-11 );
  stan::math::set_zero_all_adjoints();
  y1(2).grad(theta, g1);
  EXPECT_FLOAT_EQ(g1[0], -840);
  EXPECT_NEAR(g1[1], 0.0, 1.e-12);
  EXPECT_NEAR(g1[2], 0.0, 1.e-12);
  EXPECT_FLOAT_EQ(g1[3], 120 );
  // EXPECT_FLOAT_EQ(g1[4], 1.02318153949e-12); //FIXME: fail on g++ but not clang
  EXPECT_NEAR(g1[4], 0.0, 1.e-11);

  // FIXME: check the results
  // cmt = 2;
  // y1 = model.solve(t0, amt, rate_var[cmt - 1], ii, cmt);

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
  y1 = model.solve(t0, amt, rate_var[cmt - 1], ii, cmt);
  
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
  EXPECT_NEAR(g1[1], 0.0, 1e-12);
  EXPECT_FLOAT_EQ(g1[2], 80);
  EXPECT_NEAR(g1[3], 0.0, 1e-12);
  EXPECT_FLOAT_EQ(g1[4], 0);
  stan::math::set_zero_all_adjoints();
  y1(2).grad(theta, g1);
  EXPECT_FLOAT_EQ(g1[0], -560);
  EXPECT_FLOAT_EQ(g1[1], -71.4285714286);
  EXPECT_NEAR(g1[2], 0.0, 1e-13);
  EXPECT_FLOAT_EQ(g1[3], 108.571428571);
  EXPECT_FLOAT_EQ(g1[4], 0.0);
}

TEST_F(TorstenTwoCptModelTest, long_long_ss_infusion_vs_ode) {
  double amt = 1300;
  double r = 500;
  double ii = 1.2;
  std::vector<var> par_var(to_var(par));
  PMXTwoCptModel<var> model1(par_var);
  torsten::PKODEModel<var, PMXTwoCptODE> model2(par_var, model1.ncmt(), model1.f());

  const dsolve::PMXAnalyiticalIntegrator integ1;
  const PMXOdeIntegrator<PMXVariadicOdeSystem, PMXCvodesIntegrator<CV_BDF, CV_STAGGERED>> integ2;
  torsten::PKRec<var> y1 = model1.solve(ts[0], amt, r, ii, 1, integ1);
  torsten::PKRec<var> y2 = model2.solve(ts[0], amt, r, ii, 1, integ2);
  torsten::test::test_grad(par_var, y1, y2, 1e-6, 5e-6);
}
