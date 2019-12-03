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
using refactor::PKODEModel;

TEST_F(TorstenTwoCptModelTest, ode_model_ss_bolus_vs_long_run_sd) {
  y0[0] = 150;
  y0[1] = 55;
  y0[2] = 120;

  int cmt = 0;
  const double ii = 8.5;
  PMXTwoCptODE f2cpt;
  const std::vector<double> theta{CL, Q, V2, V3, ka};
  const PMXOdeIntegrator<PkAdams> integrator_adams;
  const PMXOdeIntegrator<PkBdf> integrator_bdf;
  const PMXOdeIntegrator<PkRk45> integrator_rk45;

  auto f1 = [&](const double amt, const auto& integrator) {
    double t = t0;
    Eigen::Matrix<double, -1, 1> y = y0;
    for (int i = 0; i < 50; ++i) {
      Eigen::Matrix<double, 1, -1> yt = y.transpose();
      PKODEModel<double, double, double, double, PMXTwoCptODE> model_i(t, yt, rate, theta, f2cpt);
      double t_next = t + ii;
      Eigen::Matrix<double, -1, 1> ys = model_i.solve(t_next,integrator);
      ys(cmt - 1) += amt;
      y = ys;
      t = t_next;
    }
    // steady state solution is the end of II dosing before
    // bolus is imposed, to check that we remove the
    // bolus(added in the last iteration) from the results
    y(cmt - 1) -= amt;
    return y;
  };

  auto f2 = [&](const double amt, const auto& integrator) {
    PKODEModel<double, double, double, double, PMXTwoCptODE> model(t0, y0, rate, theta, f2cpt);
    double r = 0;
    return model.solve(amt, r, ii, cmt, integrator);
  };

  const double amt = 1000.0;

  // check each dosing on compartment, solved by different integrator
  for (int i = 0; i < 3; ++i) {
    cmt = i + 1;
    Eigen::VectorXd y1 = f1(amt, integrator_adams);
    Eigen::VectorXd y2 = f2(amt, integrator_adams);
    torsten::test::test_val(y1, y2, 1e-10, 1e-15);    
  }

  for (int i = 0; i < 3; ++i) {
    cmt = i + 1;
    Eigen::VectorXd y1 = f1(amt, integrator_bdf);
    Eigen::VectorXd y2 = f2(amt, integrator_bdf);
    torsten::test::test_val(y1, y2, 1e-10, 1e-15);    
  }

  for (int i = 0; i < 3; ++i) {
    cmt = i + 1;
    Eigen::VectorXd y1 = f1(amt, integrator_rk45);
    Eigen::VectorXd y2 = f2(amt, integrator_rk45);
    torsten::test::test_val(y1, y2, 1e-10, 1e-15);    
  }

  // cross check between integrators
  for (int i = 0; i < 3; ++i) {
    cmt = i + 1;
    Eigen::VectorXd y1 = f1(amt, integrator_bdf);
    Eigen::VectorXd y2 = f2(amt, integrator_rk45);
    torsten::test::test_val(y1, y2, 1e-7, 1e-15);    
  }

  for (int i = 0; i < 3; ++i) {
    cmt = i + 1;
    Eigen::VectorXd y1 = f1(amt, integrator_bdf);
    Eigen::VectorXd y2 = f2(amt, integrator_adams);
    torsten::test::test_val(y1, y2, 1e-7, 1e-15);    
  }
}

TEST_F(TorstenTwoCptModelTest, ode_model_ss_bolus_grad_vs_long_run_sd) {
  y0[0] = 150;
  y0[1] = 55;
  y0[2] = 120;

  int cmt = 0;
  const double ii = 8.5;
  PMXTwoCptODE f2cpt;
  const std::vector<double> theta{CL, Q, V2, V3, ka};
  const PMXOdeIntegrator<PkAdams> integrator_adams;
  const PMXOdeIntegrator<PkBdf> integrator_bdf;
  const PMXOdeIntegrator<PkRk45> integrator_rk45;

  auto f1 = [&](const var& amt, const auto& integrator) {
    double t = t0;
    Eigen::Matrix<var, -1, 1> y = y0;
    for (int i = 0; i < 50; ++i) {
      Eigen::Matrix<var, 1, -1> yt = y.transpose();
      PKODEModel<double, var, double, double, PMXTwoCptODE> model_i(t, yt, rate, theta, f2cpt);
      double t_next = t + ii;
      Eigen::Matrix<var, -1, 1> ys = model_i.solve(t_next,integrator);
      ys(cmt - 1) += amt;
      y = ys;
      t = t_next;
    }
    // steady state solution is the end of II dosing before
    // bolus is imposed, to check that we remove the
    // bolus(added in the last iteration) from the results
    y(cmt - 1) -= amt;
    return y;
  };

  auto f2 = [&](const var& amt, const auto& integrator) {
    PKODEModel<double, var, double, double, PMXTwoCptODE> model(t0, y0, rate, theta, f2cpt);
    double r = 0;
    return model.solve(amt, r, ii, cmt, integrator);
  };

  std::vector<var> amt_vec{1000.0};

  for (int i = 0; i < 3; ++i) {
    cmt = i + 1;
    Eigen::Matrix<var, -1, 1> y1 = f1(amt_vec[0], integrator_adams);
    Eigen::Matrix<var, -1, 1> y2 = f2(amt_vec[0], integrator_adams);
    torsten::test::test_grad(amt_vec, y1, y2, 1.e-8, 1.e-11);    
  }

  for (int i = 0; i < 3; ++i) {
    cmt = i + 1;
    Eigen::Matrix<var, -1, 1> y1 = f1(amt_vec[0], integrator_bdf);
    Eigen::Matrix<var, -1, 1> y2 = f2(amt_vec[0], integrator_bdf);
    torsten::test::test_grad(amt_vec, y1, y2, 1.e-8, 1.e-11);    
  }

  for (int i = 0; i < 3; ++i) {
    cmt = i + 1;
    Eigen::Matrix<var, -1, 1> y1 = f1(amt_vec[0], integrator_rk45);
    Eigen::Matrix<var, -1, 1> y2 = f2(amt_vec[0], integrator_rk45);
    torsten::test::test_grad(amt_vec, y1, y2, 1.e-8, 2.e-11);    
  }

  // cross check
  for (int i = 0; i < 3; ++i) {
    cmt = i + 1;
    Eigen::Matrix<var, -1, 1> y1 = f1(amt_vec[0], integrator_bdf);
    Eigen::Matrix<var, -1, 1> y2 = f2(amt_vec[0], integrator_adams);
    torsten::test::test_grad(amt_vec, y1, y2, 3.e-7, 5.e-10);    
  }

  for (int i = 0; i < 3; ++i) {
    cmt = i + 1;
    Eigen::Matrix<var, -1, 1> y1 = f1(amt_vec[0], integrator_bdf);
    Eigen::Matrix<var, -1, 1> y2 = f2(amt_vec[0], integrator_rk45);
    torsten::test::test_grad(amt_vec, y1, y2, 3.e-7, 5.e-10);    
  }
}

TEST_F(TorstenTwoCptModelTest, ode_model_amt_data_ss_infusion_vs_long_run_sd) {
  y0[0] = 150;
  y0[1] = 55;
  y0[2] = 120;

  int cmt = 0;
  const double amt = 1100.0;
  const double ii = 8.5;
  PMXTwoCptODE f2cpt;
  const std::vector<double> theta{CL, Q, V2, V3, ka};
  const PMXOdeIntegrator<PkAdams> integrator_adams;
  const PMXOdeIntegrator<PkBdf> integrator_bdf;
  const PMXOdeIntegrator<PkRk45> integrator_rk45;

  auto f1 = [&](std::vector<double>& rate_vec, const auto& integrator) {
    using model_t = PKODEModel<double, double, double, double, PMXTwoCptODE>;
    double t = t0;
    Eigen::Matrix<double, -1, 1> y = y0;
    double t_infus = amt/rate_vec[cmt - 1];
    const std::vector<double> rate_zero(3, 0.0);
    for (int i = 0; i < 50; ++i) {
      Eigen::Matrix<double, 1, -1> yt;
      Eigen::Matrix<double, -1, 1> ys;
      yt = y.transpose();
      model_t model_i(t, yt, rate_vec, theta, f2cpt);
      double t_next = t + t_infus;
      ys = model_i.solve(t_next, integrator);
      yt = ys.transpose();
      t = t_next;
      model_t model_j(t, yt, rate_zero, theta, f2cpt);
      t_next = t + ii - t_infus;
      ys = model_j.solve(t_next, integrator);
      y = ys;
    }
    return y;
  };

  auto f2 = [&](std::vector<double>& rate_vec, const auto& integrator) {
    PKODEModel<double, double, double, double, PMXTwoCptODE> model(t0, y0, rate, theta, f2cpt);
    return model.solve(amt, rate_vec[cmt - 1], ii, cmt, integrator);
  };

  std::vector<double> rate_vec(3, 0.0);

  for (int i = 0; i < 3; ++i) {
    cmt = i + 1;
    rate_vec[cmt - 1] = 330.0;
    Eigen::VectorXd y1 = f1(rate_vec , integrator_adams);
    Eigen::VectorXd y2 = f2(rate_vec , integrator_adams);
    torsten::test::test_val(y1, y2, 1.e-9, 1.e-12);
    rate_vec[cmt - 1] = 0.0;
  }

  for (int i = 0; i < 3; ++i) {
    cmt = i + 1;
    rate_vec[cmt - 1] = 330.0;
    Eigen::VectorXd y1 = f1(rate_vec , integrator_bdf);
    Eigen::VectorXd y2 = f2(rate_vec , integrator_bdf);
    torsten::test::test_val(y1, y2, 1.e-9, 1.e-12);
    rate_vec[cmt - 1] = 0.0;
  }

  for (int i = 0; i < 3; ++i) {
    cmt = i + 1;
    rate_vec[cmt - 1] = 330.0;
    Eigen::VectorXd y1 = f1(rate_vec , integrator_rk45);
    Eigen::VectorXd y2 = f2(rate_vec , integrator_rk45);
    torsten::test::test_val(y1, y2, 1.e-9, 1.e-12);
    rate_vec[cmt - 1] = 0.0;
  }

  for (int i = 0; i < 3; ++i) {
    cmt = i + 1;
    rate_vec[cmt - 1] = 330.0;
    Eigen::VectorXd y1 = f1(rate_vec , integrator_bdf);
    Eigen::VectorXd y2 = f2(rate_vec , integrator_rk45);
    torsten::test::test_val(y1, y2, 5.e-8, 1.e-12);
    rate_vec[cmt - 1] = 0.0;
  }
}

TEST_F(TorstenTwoCptModelTest, ode_model_amt_parm_ss_infusion_vs_long_run_sd) {
  y0[0] = 150;
  y0[1] = 55;
  y0[2] = 120;

  int cmt = 0;
  std::vector<var> amt{1100.0};
  const double ii = 8.5;
  PMXTwoCptODE f2cpt;
  const std::vector<double> theta{CL, Q, V2, V3, ka};
  const PMXOdeIntegrator<PkAdams> integrator_adams;
  const PMXOdeIntegrator<PkBdf> integrator_bdf;
  const PMXOdeIntegrator<PkRk45> integrator_rk45;

  auto f1 = [&](std::vector<double>& rate_vec, const auto& integrator) {
    using model_t = PKODEModel<var, var, double, double, PMXTwoCptODE>;
    var t = t0;
    Eigen::Matrix<var, -1, 1> y = y0;
    var t_infus = amt[0]/rate_vec[cmt - 1];
    const std::vector<double> rate_zero(3, 0.0);
    for (int i = 0; i < 50; ++i) {
      Eigen::Matrix<var, 1, -1> yt;
      Eigen::Matrix<var, -1, 1> ys;
      yt = y.transpose();
      model_t model_i(t, yt, rate_vec, theta, f2cpt);
      var t_next = t + t_infus;
      ys = model_i.solve(t_next, integrator);
      yt = ys.transpose();
      t = t_next;
      model_t model_j(t, yt, rate_zero, theta, f2cpt);
      t_next = t + ii - t_infus;
      ys = model_j.solve(t_next, integrator);
      y = ys;
    }
    return y;
  };

  auto f2 = [&](std::vector<double>& rate_vec, const auto& integrator) {
    PKODEModel<double, double, double, double, PMXTwoCptODE> model(t0, y0, rate, theta, f2cpt);
    return model.solve(amt[0], rate_vec[cmt - 1], ii, cmt, integrator);
  };

  std::vector<double> rate_vec(3, 0.0);

  for (int i = 0; i < 0; ++i) {
    cmt = i + 1;
    rate_vec[cmt - 1] = 330.0;
    Eigen::Matrix<var, -1, 1> y1 = f1(rate_vec, integrator_adams);
    Eigen::Matrix<var, -1, 1> y2 = f2(rate_vec, integrator_adams);
    torsten::test::test_grad(amt, y1, y2, 5.e-9, 5.e-12);
    rate_vec[cmt - 1] = 0.0;    
  }

  for (int i = 0; i < 0; ++i) {
    cmt = i + 1;
    rate_vec[cmt - 1] = 330.0;
    Eigen::Matrix<var, -1, 1> y1 = f1(rate_vec, integrator_bdf);
    Eigen::Matrix<var, -1, 1> y2 = f2(rate_vec, integrator_bdf);
    torsten::test::test_grad(amt, y1, y2, 5.e-9, 5.e-12);
    rate_vec[cmt - 1] = 0.0;    
  }

  for (int i = 0; i < 0; ++i) {
    cmt = i + 1;
    rate_vec[cmt - 1] = 330.0;
    Eigen::Matrix<var, -1, 1> y1 = f1(rate_vec, integrator_rk45);
    Eigen::Matrix<var, -1, 1> y2 = f2(rate_vec, integrator_rk45);
    torsten::test::test_grad(amt, y1, y2, 5.e-9, 5.e-12);
    rate_vec[cmt - 1] = 0.0;    
  }

  // cross check
  for (int i = 0; i < 0; ++i) {
    cmt = i + 1;
    rate_vec[cmt - 1] = 330.0;
    Eigen::Matrix<var, -1, 1> y1 = f1(rate_vec, integrator_bdf);
    Eigen::Matrix<var, -1, 1> y2 = f2(rate_vec, integrator_adams);
    torsten::test::test_grad(amt, y1, y2, 5.e-9, 5.e-12);
    rate_vec[cmt - 1] = 0.0;    
  }

  for (int i = 0; i < 0; ++i) {
    cmt = i + 1;
    rate_vec[cmt - 1] = 330.0;
    Eigen::Matrix<var, -1, 1> y1 = f1(rate_vec, integrator_bdf);
    Eigen::Matrix<var, -1, 1> y2 = f2(rate_vec, integrator_rk45);
    torsten::test::test_grad(amt, y1, y2, 5.e-9, 5.e-12);
    rate_vec[cmt - 1] = 0.0;    
  }
}

TEST_F(TorstenTwoCptModelTest, ode_model_amt_data_ss_infusion_grad_vs_long_run_sd) {
  y0[0] = 150;
  y0[1] = 55;
  y0[2] = 120;

  int cmt = 0;
  double amt = 1100.0;
  const double ii = 8.5;
  PMXTwoCptODE f2cpt;
  const std::vector<double> theta{CL, Q, V2, V3, ka};
  const PMXOdeIntegrator<PkAdams> integrator_adams;
  const PMXOdeIntegrator<PkBdf> integrator_bdf;
  const PMXOdeIntegrator<PkRk45> integrator_rk45;

  auto f1 = [&](const var& r, const auto& integrator) {
    using model_t = PKODEModel<var, var, var, double, PMXTwoCptODE>;
    var t = t0;
    Eigen::Matrix<var, -1, 1> y = y0;
    var t_infus = amt/r;
    const std::vector<var> rate_zero(3, 0.0);
    std::vector<var> rate_vec(3, 0.0);
    rate_vec[cmt - 1] = r;
    for (int i = 0; i < 50; ++i) {
      Eigen::Matrix<var, 1, -1> yt;
      Eigen::Matrix<var, -1, 1> ys;
      yt = y.transpose();
      model_t model_i(t, yt, rate_vec, theta, f2cpt);
      var t_next = t + t_infus;
      ys = model_i.solve(t_next, integrator);
      yt = ys.transpose();
      t = t_next;
      model_t model_j(t, yt, rate_zero, theta, f2cpt);
      t_next = t + ii - t_infus;
      ys = model_j.solve(t_next, integrator);
      y = ys;
    }
    return y;
  };

  auto f2 = [&](const var& r, const auto& integrator) {
    PKODEModel<double, double, double, double, PMXTwoCptODE> model(t0, y0, rate, theta, f2cpt);
    return model.solve(amt, r, ii, cmt, integrator);
  };

  for (int i = 0; i < 3; ++i) {
    cmt = i + 1;
    const var r = 330.0;
    std::vector<var> rvec{r};
    Eigen::Matrix<var, -1, 1> y1 = f1(r, integrator_adams);
    Eigen::Matrix<var, -1, 1> y2 = f2(r, integrator_adams);
    torsten::test::test_grad(rvec, y1, y2, 5.e-9, 3.e-12);    
  }

  for (int i = 0; i < 3; ++i) {
    cmt = i + 1;
    const var r = 330.0;
    std::vector<var> rvec{r};
    Eigen::Matrix<var, -1, 1> y1 = f1(r, integrator_bdf);
    Eigen::Matrix<var, -1, 1> y2 = f2(r, integrator_bdf);
    torsten::test::test_grad(rvec, y1, y2, 5.e-9, 3.e-12);    
  }

  for (int i = 0; i < 3; ++i) {
    cmt = i + 1;
    const var r = 330.0;
    std::vector<var> rvec{r};
    Eigen::Matrix<var, -1, 1> y1 = f1(r, integrator_rk45);
    Eigen::Matrix<var, -1, 1> y2 = f2(r, integrator_rk45);
    torsten::test::test_grad(rvec, y1, y2, 5.e-9, 3.e-12);    
  }
}

TEST_F(TorstenTwoCptModelTest, ode_model_amt_data_ss_const_infusion_grad_vs_long_run_sd) {
  y0[0] = 150;
  y0[1] = 55;
  y0[2] = 120;

  int cmt = 0;
  double amt = 0.0;
  const double ii = 0.0;
  PMXTwoCptODE f2cpt;
  const std::vector<double> theta{CL, Q, V2, V3, ka};
  const PMXOdeIntegrator<PkAdams> integrator_adams;
  const PMXOdeIntegrator<PkBdf> integrator_bdf;
  const PMXOdeIntegrator<PkRk45> integrator_rk45;

  auto f1 = [&](const var& r, const auto& integrator) {
    using model_t = PKODEModel<double, double, var, double, PMXTwoCptODE>;
    var t = t0;
    Eigen::Matrix<var, -1, 1> y = y0;
    std::vector<var> rate_vec(3, 0.0);
    rate_vec[cmt - 1] = r;
    model_t model(t0, y0, rate_vec, theta, f2cpt);
    double t_next = 5.0e2;
    return model.solve(t_next, integrator);
  };

  auto f2 = [&](const var& r, const auto& integrator) {
    PKODEModel<double, double, double, double, PMXTwoCptODE> model(t0, y0, rate, theta, f2cpt);
    return model.solve(amt, r, ii, cmt, integrator);
  };

  const var r = 330.0;
  std::vector<var> rvec{r};

  for (int i = 0; i < 3; ++i) {
    cmt = i + 1;
    Eigen::Matrix<var, -1, 1> y1 = f1(r, integrator_adams);
    Eigen::Matrix<var, -1, 1> y2 = f2(r, integrator_adams);
    torsten::test::test_grad(rvec, y1, y2, 1.e-7, 1.e-9);    
  }

  for (int i = 0; i < 3; ++i) {
    cmt = i + 1;
    Eigen::Matrix<var, -1, 1> y1 = f1(r, integrator_bdf);
    Eigen::Matrix<var, -1, 1> y2 = f2(r, integrator_bdf);
    torsten::test::test_grad(rvec, y1, y2, 3e-7, 1.e-9);    
  }

  for (int i = 0; i < 3; ++i) {
    cmt = i + 1;
    Eigen::Matrix<var, -1, 1> y1 = f1(r, integrator_rk45);
    Eigen::Matrix<var, -1, 1> y2 = f2(r, integrator_rk45);
    torsten::test::test_grad(rvec, y1, y2, 1e-7, 1.e-9);    
  }
}

TEST_F(TorstenTwoCptModelTest, ode_model_rate_param_ss_bolus_vs_long_run_sd) {
  y0[0] = 150;
  y0[1] = 55;
  y0[2] = 120;

  int cmt = 0;
  const double ii = 8.5;
  PMXTwoCptODE f2cpt;
  const std::vector<double> theta{CL, Q, V2, V3, ka};
  const PMXOdeIntegrator<PkAdams> integrator_adams;
  const PMXOdeIntegrator<PkBdf> integrator_bdf;
  const PMXOdeIntegrator<PkRk45> integrator_rk45;

  auto f1 = [&](const double amt, const auto& integrator) {
    double t = t0;
    Eigen::Matrix<double, -1, 1> y = y0;
    for (int i = 0; i < 50; ++i) {
      Eigen::Matrix<double, 1, -1> yt = y.transpose();
      PKODEModel<double, double, double, double, PMXTwoCptODE> model_i(t, yt, rate, theta, f2cpt);
      double t_next = t + ii;
      Eigen::Matrix<double, -1, 1> ys = model_i.solve(t_next,integrator);
      ys(cmt - 1) += amt;
      y = ys;
      t = t_next;
    }
    // steady state solution is the end of II dosing before
    // bolus is imposed, to check that we remove the
    // bolus(added in the last iteration) from the results
    y(cmt - 1) -= amt;
    return y;
  };

  auto f2 = [&](const double amt, const auto& integrator) {
    PKODEModel<double, double, double, double, PMXTwoCptODE> model(t0, y0, rate, theta, f2cpt);
    var r = 0.0;
    return model.solve(amt, r, ii, cmt, integrator);
  };

  const double amt = 1000.0;

  for (int i = 0; i < 3; ++i) {
    cmt = i + 1;
    Eigen::Matrix<double, -1, 1> y1 = f1(amt, integrator_adams);
    Eigen::Matrix<var, -1, 1>    y2 = f2(amt, integrator_adams);
    torsten::test::test_val(y1, stan::math::value_of(y2), 1e-10, 1e-15);    
  }

  for (int i = 0; i < 3; ++i) {
    cmt = i + 1;
    Eigen::Matrix<double, -1, 1> y1 = f1(amt, integrator_bdf);
    Eigen::Matrix<var, -1, 1>    y2 = f2(amt, integrator_bdf);
    torsten::test::test_val(y1, stan::math::value_of(y2), 1e-10, 1e-15);    
  }

  for (int i = 0; i < 3; ++i) {
    cmt = i + 1;
    Eigen::Matrix<double, -1, 1> y1 = f1(amt, integrator_rk45);
    Eigen::Matrix<var, -1, 1>    y2 = f2(amt, integrator_rk45);
    torsten::test::test_val(y1, stan::math::value_of(y2), 1e-10, 1e-15);    
  }
}

TEST_F(TorstenTwoCptModelTest, ode_model_rate_param_ss_bolus_grad_vs_long_run_sd) {
  y0[0] = 150;
  y0[1] = 55;
  y0[2] = 120;

  int cmt = 0;
  const double ii = 8.5;
  PMXTwoCptODE f2cpt;
  const std::vector<double> theta{CL, Q, V2, V3, ka};
  const PMXOdeIntegrator<PkAdams> integrator_adams;
  const PMXOdeIntegrator<PkBdf> integrator_bdf;
  const PMXOdeIntegrator<PkRk45> integrator_rk45;

  auto f1 = [&](const var& amt, const auto& integrator) {
    double t = t0;
    Eigen::Matrix<var, -1, 1> y = y0;
    for (int i = 0; i < 50; ++i) {
      Eigen::Matrix<var, 1, -1> yt = y.transpose();
      PKODEModel<double, var, double, double, PMXTwoCptODE> model_i(t, yt, rate, theta, f2cpt);
      double t_next = t + ii;
      Eigen::Matrix<var, -1, 1> ys = model_i.solve(t_next,integrator);
      ys(cmt - 1) += amt;
      y = ys;
      t = t_next;
    }
    // steady state solution is the end of II dosing before
    // bolus is imposed, to check that we remove the
    // bolus(added in the last iteration) from the results
    y(cmt - 1) -= amt;
    return y;
  };

  auto f2 = [&](const var& amt, const auto& integrator) {
    PKODEModel<double, double, double, double, PMXTwoCptODE> model(t0, y0, rate, theta, f2cpt);
    var r = 0.0;
    return model.solve(amt, r, ii, cmt, integrator);
  };

  std::vector<var> amt_vec{1000.0};

  for (int i = 0; i < 3; ++i) {
    cmt = i + 1;
    Eigen::Matrix<var, -1, 1> y1 = f1(amt_vec[0], integrator_adams);
    Eigen::Matrix<var, -1, 1> y2 = f2(amt_vec[0], integrator_adams);
    torsten::test::test_grad(amt_vec, y1, y2, 5.5e-9, 1e-11);    
  }

  for (int i = 0; i < 3; ++i) {
    cmt = i + 1;
    Eigen::Matrix<var, -1, 1> y1 = f1(amt_vec[0], integrator_bdf);
    Eigen::Matrix<var, -1, 1> y2 = f2(amt_vec[0], integrator_bdf);
    torsten::test::test_grad(amt_vec, y1, y2, 5.5e-9, 1e-11);    
  }

  for (int i = 0; i < 3; ++i) {
    cmt = i + 1;
    Eigen::Matrix<var, -1, 1> y1 = f1(amt_vec[0], integrator_rk45);
    Eigen::Matrix<var, -1, 1> y2 = f2(amt_vec[0], integrator_rk45);
    torsten::test::test_grad(amt_vec, y1, y2, 5.5e-9, 2e-11);    
  }
}

TEST_F(TorstenTwoCptModelTest, ode_model_amt_param_ss_infusion_grad_vs_long_run_sd) {
  y0[0] = 150;
  y0[1] = 55;
  y0[2] = 120;

  int cmt = 0;
  var amt = 1100.0;
  const double ii = 8.5;
  PMXTwoCptODE f2cpt;
  const std::vector<double> theta{CL, Q, V2, V3, ka};
  const PMXOdeIntegrator<PkAdams> integrator_adams;
  const PMXOdeIntegrator<PkBdf> integrator_bdf;
  const PMXOdeIntegrator<PkRk45> integrator_rk45;

  auto f1 = [&](const var& r, const auto& integrator) {
    using model_t = PKODEModel<var, var, var, double, PMXTwoCptODE>;
    var t = t0;
    Eigen::Matrix<var, -1, 1> y = y0;
    var t_infus = amt/r;
    const std::vector<var> rate_zero(3, 0.0);
    std::vector<var> rate_vec(3, 0.0);
    rate_vec[cmt - 1] = r;
    for (int i = 0; i < 50; ++i) {
      Eigen::Matrix<var, 1, -1> yt;
      Eigen::Matrix<var, -1, 1> ys;
      yt = y.transpose();
      model_t model_i(t, yt, rate_vec, theta, f2cpt);
      var t_next = t + t_infus;
      ys = model_i.solve(t_next, integrator);
      yt = ys.transpose();
      t = t_next;
      model_t model_j(t, yt, rate_zero, theta, f2cpt);
      t_next = t + ii - t_infus;
      ys = model_j.solve(t_next, integrator);
      y = ys;
    }
    return y;
  };

  auto f2 = [&](const var& r, const auto& integrator) {
    PKODEModel<double, double, double, double, PMXTwoCptODE> model(t0, y0, rate, theta, f2cpt);
    return model.solve(amt, r, ii, cmt, integrator);
  };

  const var r = 330.0;
  std::vector<var> params{amt, r};

  for (int i = 0; i < 3; ++i) {
    cmt = i + 1;
    Eigen::Matrix<var, -1, 1> y1 = f1(r, integrator_adams);
    Eigen::Matrix<var, -1, 1> y2 = f2(r, integrator_adams);
    torsten::test::test_grad(params, y1, y2, 5.e-9, 5.e-12);    
  }

  for (int i = 0; i < 3; ++i) {
    cmt = i + 1;
    Eigen::Matrix<var, -1, 1> y1 = f1(r, integrator_bdf);
    Eigen::Matrix<var, -1, 1> y2 = f2(r, integrator_bdf);
    torsten::test::test_grad(params, y1, y2, 5.e-9, 5.e-12);    
  }

  for (int i = 0; i < 3; ++i) {
    cmt = i + 1;
    Eigen::Matrix<var, -1, 1> y1 = f1(r, integrator_rk45);
    Eigen::Matrix<var, -1, 1> y2 = f2(r, integrator_rk45);
    torsten::test::test_grad(params, y1, y2, 5.e-9, 5.e-12);    
  }
}

TEST_F(TorstenTwoCptModelTest, ode_model_amt_param_ss_const_infusion_grad_vs_long_run_sd) {
  y0[0] = 150;
  y0[1] = 55;
  y0[2] = 120;

  int cmt = 0;
  var amt = 0.0;
  const double ii = 0.0;
  PMXTwoCptODE f2cpt;
  const std::vector<double> theta{CL, Q, V2, V3, ka};
  const PMXOdeIntegrator<PkAdams> integrator_adams;
  const PMXOdeIntegrator<PkBdf> integrator_bdf;
  const PMXOdeIntegrator<PkRk45> integrator_rk45;

  auto f1 = [&](const var& r, const auto& integrator) {
    using model_t = PKODEModel<double, double, var, double, PMXTwoCptODE>;
    var t = t0;
    Eigen::Matrix<var, -1, 1> y = y0;
    var t_infus = amt/r;
    std::vector<var> rate_vec(3, 0.0);
    rate_vec[cmt - 1] = r;
    model_t model(t0, y0, rate_vec, theta, f2cpt);
    double t_next = 5.0e2;
    return model.solve(t_next, integrator);
  };

  auto f2 = [&](const var& r, const auto& integrator) {
    PKODEModel<double, double, double, double, PMXTwoCptODE> model(t0, y0, rate, theta, f2cpt);
    return model.solve(amt, r, ii, cmt, integrator);
  };

  const var r = 330.0;
  std::vector<var> rvec{r};

  for (int i = 0; i < 3; ++i) {
    cmt = i + 1;
    Eigen::Matrix<var, -1, 1> y1 = f1(r, integrator_adams);
    Eigen::Matrix<var, -1, 1> y2 = f2(r, integrator_adams);
    torsten::test::test_grad(rvec, y1, y2, 3.e-7, 1.e-9);    
  }

  for (int i = 0; i < 3; ++i) {
    cmt = i + 1;
    Eigen::Matrix<var, -1, 1> y1 = f1(r, integrator_bdf);
    Eigen::Matrix<var, -1, 1> y2 = f2(r, integrator_bdf);
    torsten::test::test_grad(rvec, y1, y2, 3.e-7, 1.e-9);    
  }

  for (int i = 0; i < 3; ++i) {
    cmt = i + 1;
    Eigen::Matrix<var, -1, 1> y1 = f1(r, integrator_rk45);
    Eigen::Matrix<var, -1, 1> y2 = f2(r, integrator_rk45);
    torsten::test::test_grad(rvec, y1, y2, 3.e-7, 1.e-9);    
  }
}

TEST_F(TorstenTwoCptModelTest, ode_model_ss_bolus_theta_grad_vs_long_run_sd) {
  y0[0] = 150;
  y0[1] = 55;
  y0[2] = 120;

  int cmt = 0;
  const double ii = 8.5;
  PMXTwoCptODE f2cpt;
  double amt = 1000.0;
  std::vector<var> params{CL, Q, V2, V3, ka, amt};
  std::vector<var> theta(params.begin(), params.begin() + 5);
  const PMXOdeIntegrator<PkAdams> integrator_adams;
  const PMXOdeIntegrator<PkBdf> integrator_bdf;
  const PMXOdeIntegrator<PkRk45> integrator_rk45;

  auto f1 = [&](const auto& integrator) {
    double t = t0;
    Eigen::Matrix<var, -1, 1> y = y0;
    for (int i = 0; i < 50; ++i) {
      Eigen::Matrix<var, 1, -1> yt = y.transpose();
      PKODEModel<double, var, double, var, PMXTwoCptODE> model_i(t, yt, rate, theta, f2cpt);
      double t_next = t + ii;
      Eigen::Matrix<var, -1, 1> ys = model_i.solve(t_next,integrator);
      ys(cmt - 1) += amt;
      y = ys;
      t = t_next;
    }
    // steady state solution is the end of II dosing before
    // bolus is imposed, to check that we remove the
    // bolus(added in the last iteration) from the results
    y(cmt - 1) -= amt;
    return y;
  };

  auto f2 = [&](const auto& integrator) {
    PKODEModel<double, double, double, var, PMXTwoCptODE> model(t0, y0, rate, theta, f2cpt);
    double r = 0;
    return model.solve(amt, r, ii, cmt, integrator);
  };


  auto f3 = [&](const auto& integrator) {
    double t = t0;
    Eigen::Matrix<var, -1, 1> y = y0;
    for (int i = 0; i < 50; ++i) {
      Eigen::Matrix<var, 1, -1> yt = y.transpose();
      PKODEModel<double, var, double, var, PMXTwoCptODE> model_i(t, yt, rate, theta, f2cpt);
      double t_next = t + ii;
      Eigen::Matrix<var, -1, 1> ys = model_i.solve(t_next,integrator);
      ys(cmt - 1) += params.back();
      y = ys;
      t = t_next;
    }
    // steady state solution is the end of II dosing before
    // bolus is imposed, to check that we remove the
    // bolus(added in the last iteration) from the results
    y(cmt - 1) -= params.back();
    return y;
  };

  auto f4 = [&](const auto& integrator) {
    PKODEModel<double, double, double, var, PMXTwoCptODE> model(t0, y0, rate, theta, f2cpt);
    double r = 0;
    return model.solve(params.back(), r, ii, cmt, integrator);
  };

  for (int i = 0; i < 3; ++i) {
    cmt = i + 1;
    Eigen::Matrix<var, -1, 1> y1 = f1(integrator_adams);
    Eigen::Matrix<var, -1, 1> y2 = f2(integrator_adams);
    Eigen::Matrix<var, -1, 1> y3 = f3(integrator_adams);
    Eigen::Matrix<var, -1, 1> y4 = f4(integrator_adams);
    torsten::test::test_grad(params, y1, y2, 8e-9, 2e-8);    
    torsten::test::test_grad(params, y3, y4, 8e-9, 2e-8);    
  }

  for (int i = 0; i < 3; ++i) {
    cmt = i + 1;
    Eigen::Matrix<var, -1, 1> y1 = f1(integrator_bdf);
    Eigen::Matrix<var, -1, 1> y2 = f2(integrator_bdf);
    Eigen::Matrix<var, -1, 1> y3 = f3(integrator_bdf);
    Eigen::Matrix<var, -1, 1> y4 = f4(integrator_bdf);
    torsten::test::test_grad(params, y1, y2, 8e-9, 2e-8);    
    torsten::test::test_grad(params, y3, y4, 8e-9, 2e-8);    
  }

  for (int i = 0; i < 3; ++i) {
    cmt = i + 1;
    Eigen::Matrix<var, -1, 1> y1 = f1(integrator_rk45);
    Eigen::Matrix<var, -1, 1> y2 = f2(integrator_rk45);
    Eigen::Matrix<var, -1, 1> y3 = f3(integrator_rk45);
    Eigen::Matrix<var, -1, 1> y4 = f4(integrator_rk45);
    torsten::test::test_grad(params, y1, y2, 2e-8, 2e-8);    
    torsten::test::test_grad(params, y3, y4, 2e-8, 2e-8);    
  }
}

TEST_F(TorstenTwoCptModelTest, ode_model_ss_infusion_theta_grad_vs_long_run_sd) {
  y0[0] = 150;
  y0[1] = 55;
  y0[2] = 120;

  int cmt = 0;
  double amt = 1100.0;
  const double ii = 8.5;
  PMXTwoCptODE f2cpt;
  std::vector<var> params{CL, Q, V2, V3, ka, amt, 330};
  std::vector<var> theta(params.begin(), params.begin() + 5);
  const PMXOdeIntegrator<PkAdams> integrator_adams;
  const PMXOdeIntegrator<PkBdf> integrator_bdf;
  const PMXOdeIntegrator<PkRk45> integrator_rk45;

  auto f1 = [&](const auto& integrator) {
    using model_t = PKODEModel<var, var, var, var, PMXTwoCptODE>;
    var& r = params.back();
    var t = t0;
    Eigen::Matrix<var, -1, 1> y = y0;
    var t_infus = amt/r;
    const std::vector<var> rate_zero(3, 0.0);
    std::vector<var> rate_vec(3, 0.0);
    rate_vec[cmt - 1] = r;
    for (int i = 0; i < 50; ++i) {
      Eigen::Matrix<var, 1, -1> yt;
      Eigen::Matrix<var, -1, 1> ys;
      yt = y.transpose();
      model_t model_i(t, yt, rate_vec, theta, f2cpt);
      var t_next = t + t_infus;
      ys = model_i.solve(t_next, integrator);
      yt = ys.transpose();
      t = t_next;
      model_t model_j(t, yt, rate_zero, theta, f2cpt);
      t_next = t + ii - t_infus;
      ys = model_j.solve(t_next, integrator);
      y = ys;
    }
    return y;
  };

  auto f2 = [&](const auto& integrator) {
    PKODEModel<double, double, double, var, PMXTwoCptODE> model(t0, y0, rate, theta, f2cpt);
    var& r = params.back();
    return model.solve(amt, r, ii, cmt, integrator);
  };

  auto f3 = [&](const auto& integrator) {
    using model_t = PKODEModel<var, var, var, var, PMXTwoCptODE>;
    var& r = params.back();
    var t = t0;
    Eigen::Matrix<var, -1, 1> y = y0;
    var t_infus = params[5]/r;
    const std::vector<var> rate_zero(3, 0.0);
    std::vector<var> rate_vec(3, 0.0);
    rate_vec[cmt - 1] = r;
    for (int i = 0; i < 50; ++i) {
      Eigen::Matrix<var, 1, -1> yt;
      Eigen::Matrix<var, -1, 1> ys;
      yt = y.transpose();
      model_t model_i(t, yt, rate_vec, theta, f2cpt);
      var t_next = t + t_infus;
      ys = model_i.solve(t_next, integrator);
      yt = ys.transpose();
      t = t_next;
      model_t model_j(t, yt, rate_zero, theta, f2cpt);
      t_next = t + ii - t_infus;
      ys = model_j.solve(t_next, integrator);
      y = ys;
    }
    return y;
  };

  auto f4 = [&](const auto& integrator) {
    PKODEModel<double, double, double, var, PMXTwoCptODE> model(t0, y0, rate, theta, f2cpt);
    var& r = params.back();
    return model.solve(params[5], r, ii, cmt, integrator);
  };

  for (int i = 0; i < 3; ++i) {
    cmt = i + 1;
    Eigen::Matrix<var, -1, 1> y1 = f1(integrator_adams);
    Eigen::Matrix<var, -1, 1> y2 = f2(integrator_adams);
    Eigen::Matrix<var, -1, 1> y3 = f3(integrator_adams);
    Eigen::Matrix<var, -1, 1> y4 = f4(integrator_adams);
    torsten::test::test_grad(params, y1, y2, 5.e-9, 1.e-8);
    torsten::test::test_grad(params, y3, y4, 5.e-9, 1.e-8);    
  }

  for (int i = 0; i < 3; ++i) {
    cmt = i + 1;
    Eigen::Matrix<var, -1, 1> y1 = f1(integrator_bdf);
    Eigen::Matrix<var, -1, 1> y2 = f2(integrator_bdf);
    Eigen::Matrix<var, -1, 1> y3 = f3(integrator_bdf);
    Eigen::Matrix<var, -1, 1> y4 = f4(integrator_bdf);
    torsten::test::test_grad(params, y1, y2, 5.e-9, 1.e-8);
    torsten::test::test_grad(params, y3, y4, 5.e-9, 1.e-8);    
  }

  for (int i = 0; i < 3; ++i) {
    cmt = i + 1;
    Eigen::Matrix<var, -1, 1> y1 = f1(integrator_rk45);
    Eigen::Matrix<var, -1, 1> y2 = f2(integrator_rk45);
    Eigen::Matrix<var, -1, 1> y3 = f3(integrator_rk45);
    Eigen::Matrix<var, -1, 1> y4 = f4(integrator_rk45);
    torsten::test::test_grad(params, y1, y2, 5.e-9, 1.e-8);
    torsten::test::test_grad(params, y3, y4, 5.e-9, 1.e-8);    
  }
}

TEST_F(TorstenTwoCptModelTest, ode_model_const_infusion_theta_grad_vs_long_run_sd) {
  y0[0] = 150;
  y0[1] = 55;
  y0[2] = 120;

  int cmt = 0;
  double amt = 0.0;
  const double ii = 0.0;
  PMXTwoCptODE f2cpt;
  std::vector<var> params{CL, Q, V2, V3, ka, 330.0};
  std::vector<var> theta(params.begin(), params.begin() + 5);
  const PMXOdeIntegrator<PkBdf> integrator;

  auto f1 = [&]() {
    using model_t = PKODEModel<double, double, var, var, PMXTwoCptODE>;
    var t = t0;
    Eigen::Matrix<var, -1, 1> y = y0;
    std::vector<var> rate_vec(3, 0.0);
    rate_vec[cmt - 1] = params.back();
    model_t model(t0, y0, rate_vec, theta, f2cpt);
    double t_next = 5.0e2;
    return model.solve(t_next, integrator);
  };

  auto f2 = [&]() {
    PKODEModel<double, double, double, var, PMXTwoCptODE> model(t0, y0, rate, theta, f2cpt);
    return model.solve(amt, params.back(), ii, cmt, integrator);
  };

  {
    cmt = 1;
    Eigen::Matrix<var, -1, 1> y1 = f1();
    Eigen::Matrix<var, -1, 1> y2 = f2();
    torsten::test::test_grad(params, y1, y2, 5.e-8, 5.e-8);
  }

  {
    cmt = 2;
    Eigen::Matrix<var, -1, 1> y1 = f1();
    Eigen::Matrix<var, -1, 1> y2 = f2();
    torsten::test::test_grad(params, y1, y2, 5.e-7, 5.e-7);
  }

  {
    cmt = 3;
    Eigen::Matrix<var, -1, 1> y1 = f1();
    Eigen::Matrix<var, -1, 1> y2 = f2();
    torsten::test::test_grad(params, y1, y2, 1.e-7, 1.e-7);
  }
}
