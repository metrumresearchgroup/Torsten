#include <stan/math/rev/fun/typedefs.hpp>
#include <stan/math/torsten/events_manager.hpp>
#include <stan/math/torsten/torsten_def.hpp>
#include <stan/math/torsten/pmx_twocpt_model.hpp>
#include <stan/math/torsten/pmx_ode_model.hpp>
#include <stan/math/torsten/test/unit/pmx_twocpt_test_fixture.hpp>
#include <stan/math/torsten/test/unit/test_util.hpp>
#include <gtest/gtest.h>
#include <vector>

using std::vector;
using Eigen::Matrix;
using Eigen::Dynamic;
using stan::math::var;
using torsten::PKRec;
using torsten::EventsManager;
using torsten::PMXTwoCptModel;
using torsten::PKODEModel;

TEST_F(TorstenTwoCptTest, model_solve_d_data_only) {
  using model_t = PMXTwoCptModel<double, double, double, double>;

  const int ncmt = model_t::Ncmt;
  torsten::PKRec<double> init(ncmt);
  init << 0, 100, 0;

  double t = time[1], t1 = t + 0.1;
  std::vector<double> rate{1000., 0., 0.};

  model_t model(t, init, rate, pMatrix[0]);
  
  Eigen::VectorXd sol1 = model.solve(t1);
  Eigen::VectorXd sol2 = model.solve_d(t1);
  
  torsten::test::test_val(sol1, sol2);
}

TEST_F(TorstenTwoCptTest, model_solve_d_init_var) {
  
  using stan::math::vector_v;
  using stan::math::var;
  using model_t = PMXTwoCptModel<double, var, double, double>;

  const int ncmt = model_t::Ncmt;
  torsten::PKRec<var> init(ncmt);
  init << 200, 100, 0;

  double t = time[1], t1 = t + 0.1;
  std::vector<double> rate{400., 1000., 0.};

  model_t model(t, init, rate, pMatrix[0]);

  std::vector<var> vars(model.vars(t1));
  EXPECT_EQ(vars.size(), init.size());

  vector_v sol1 = model.solve(t1);
  Eigen::VectorXd sol2_d = model.solve_d(t1);
  vector_v sol2 = torsten::mpi::precomputed_gradients(sol2_d, vars);
  torsten::test::test_grad(vars, sol1, sol2, 1.E-8, 1.E-5);
}

TEST_F(TorstenTwoCptTest, model_solve_d_rate_var) {
  
  using stan::math::vector_v;
  using stan::math::var;
  using model_t = PMXTwoCptModel<double, double, var, double>;

  const int ncmt = model_t::Ncmt;
  torsten::PKRec<double> init(ncmt);
  init << 200, 100, 0;

  double t = time[1], t1 = t + 0.1;
  std::vector<var> rate{400., 1000., 0.};

  model_t model(t, init, rate, pMatrix[0]);

  std::vector<var> vars(model.vars(t1));
  EXPECT_EQ(vars.size(), rate.size());

  vector_v sol1 = model.solve(t1);
  Eigen::VectorXd sol2_d = model.solve_d(t1);
  vector_v sol2 = torsten::mpi::precomputed_gradients(sol2_d, vars);
  torsten::test::test_grad(vars, sol1, sol2, 1.E-8, 1.E-5);
}

TEST_F(TorstenTwoCptTest, model_solve_d_par_var) {
  using stan::math::vector_v;
  using stan::math::var;
  using model_t = PMXTwoCptModel<double, double, double, var>;

  const int ncmt = model_t::Ncmt;
  torsten::PKRec<double> init(ncmt);
  init << 200, 100, 0;

  double t = time[1], t1 = t + 0.1;
  std::vector<double> rate{400., 1000., 0.};

  std::vector<var> pars(stan::math::to_var(pMatrix[0]));
  model_t model(t, init, rate, pars);

  std::vector<var> vars(model.vars(t1));
  EXPECT_EQ(vars.size(), pars.size());

  vector_v sol1 = model.solve(t1);
  Eigen::VectorXd sol2_d = model.solve_d(t1);
  vector_v sol2 = torsten::mpi::precomputed_gradients(sol2_d, vars);
  torsten::test::test_grad(vars, sol1, sol2, 1.E-8, 1.E-5);
}

TEST_F(TorstenTwoCptTest, model_solve_d_rate_par_var) {
  
  using stan::math::vector_v;
  using stan::math::var;
  using model_t = PMXTwoCptModel<double, double, var, var>;

  const int ncmt = model_t::Ncmt;
  torsten::PKRec<double> init(ncmt);
  init << 200, 100, 0;

  double t = time[1], t1 = t + 0.1;
  std::vector<var> rate{400., 1000., 0.};

  std::vector<var> pars(stan::math::to_var(pMatrix[0]));
  model_t model(t, init, rate, pars);

  std::vector<var> vars(model.vars(t1));
  EXPECT_EQ(vars.size(), pars.size() + rate.size());

  vector_v sol1 = model.solve(t1);
  Eigen::VectorXd sol2_d = model.solve_d(t1);
  vector_v sol2 = torsten::mpi::precomputed_gradients(sol2_d, vars);
  torsten::test::test_grad(vars, sol1, sol2, 1.E-8, 1.E-5);
}

TEST_F(TorstenTwoCptTest, model_solve_d_data_only_ss) {
  
  using model_t = PMXTwoCptModel<double, double, double, double>;

  const int ncmt = model_t::Ncmt;
  torsten::PKRec<double> init(ncmt);
  init << 0, 100, 0;

  double t = time[1];
  std::vector<double> rate{1000., 0., 0.};

  model_t model(t, init, rate, pMatrix[0]);
  
  double a = 1500, r = 100, ii = 18.0;
  int cmt = 1;
  Eigen::VectorXd sol1 = model.solve(a, r, ii, cmt);
  Eigen::VectorXd sol2 = model.solve_d(a, r, ii, cmt);
  torsten::test::test_val(sol1, sol2);
}

TEST_F(TorstenTwoCptTest, model_solve_d_init_var_ss) {
  
  using stan::math::var;
  using stan::math::vector_v;
  using model_t = PMXTwoCptModel<double, var, double, double>;

  const int ncmt = model_t::Ncmt;
  torsten::PKRec<var> init(ncmt);
  init << 0, 100, 0;

  double t = time[1];
  std::vector<double> rate{1000., 0., 0.};

  model_t model(t, init, rate, pMatrix[0]);
  
  double a = 1500, r = 100, ii = 18.0;
  int cmt = 1;

  // initial condition should not affect steady-state solution
  std::vector<var> vars(model.vars(a, r, ii));
  EXPECT_EQ(vars.size(), 0);

  Eigen::VectorXd sol1 = model.solve(a, r, ii, cmt);
  Eigen::VectorXd sol2 = model.solve_d(a, r, ii, cmt);
  EXPECT_EQ(sol2.size(), sol1.size());
  torsten::test::test_val(sol1, sol2);
}

TEST_F(TorstenTwoCptTest, model_solve_d_rate_var_ss) {
  
  using stan::math::var;
  using stan::math::vector_v;
  using model_t = PMXTwoCptModel<double, double, var, double>;

  const int ncmt = model_t::Ncmt;
  torsten::PKRec<double> init(ncmt);
  init << 0, 100, 0;

  double t = time[1];
  std::vector<var> rate{1000., 0., 0.};

  model_t model(t, init, rate, pMatrix[0]);
  
  double a = 1500, r = 100, ii = 18.0;
  int cmt = 1;

  // rate in model constructor should not affect steady-state solution
  std::vector<var> vars(model.vars(a, r, ii));
  EXPECT_EQ(vars.size(), 0);

  Eigen::VectorXd sol1 = model.solve(a, r, ii, cmt);
  Eigen::VectorXd sol2 = model.solve_d(a, r, ii, cmt);
  EXPECT_EQ(sol2.size(), sol1.size());
  torsten::test::test_val(sol1, sol2);
}

TEST_F(TorstenTwoCptTest, model_solve_d_par_var_ss) {
  
  using stan::math::var;
  using stan::math::vector_v;
  using model_t = PMXTwoCptModel<double, double, double, var>;

  const int ncmt = model_t::Ncmt;
  torsten::PKRec<double> init(ncmt);
  init << 0, 100, 0;

  double t = time[1];
  std::vector<double> rate{1000., 0., 0.};

  std::vector<var> pars(stan::math::to_var(pMatrix[0]));

  model_t model(t, init, rate, pars);
  
  double a = 1500, r = 100, ii = 18.0;
  int cmt = 1;

  // rate in model constructor should not affect steady-state solution
  std::vector<var> vars(model.vars(a, r, ii));
  EXPECT_EQ(vars.size(), pars.size());

  vector_v sol1 = model.solve(a, r, ii, cmt);
  Eigen::VectorXd sol2_d = model.solve_d(a, r, ii, cmt);
  vector_v sol2 = torsten::mpi::precomputed_gradients(sol2_d, vars);
  
  torsten::test::test_grad(vars, pars, sol1, sol2, 1.E-8, 1.E-5);
}

TEST_F(TorstenTwoCptTest, model_solve_d_amt_par_var_ss) {
  
  using stan::math::var;
  using stan::math::vector_v;
  using model_t = PMXTwoCptModel<double, double, double, var>;

  const int ncmt = model_t::Ncmt;
  torsten::PKRec<double> init(ncmt);
  init << 0, 100, 0;

  double t = time[1];
  std::vector<double> rate{1000., 0., 0.};

  std::vector<var> pars(stan::math::to_var(pMatrix[0]));

  model_t model(t, init, rate, pars);
  
  var a = 1500;
  double r = 100, ii = 18.0;
  int cmt = 1;

  // rate in model constructor should not affect steady-state solution
  std::vector<var> vars(model.vars(a, r, ii));
  EXPECT_EQ(vars.size(), pars.size() + 1);

  vector_v sol1 = model.solve(a, r, ii, cmt);
  Eigen::VectorXd sol2_d = model.solve_d(a, r, ii, cmt);
  vector_v sol2 = torsten::mpi::precomputed_gradients(sol2_d, vars);
  
  torsten::test::test_grad(vars, sol1, sol2, 1.E-8, 1.E-5);
}

TEST_F(TorstenTwoCptTest, model_solve_d_amt_ii_var_ss) {
  
  using stan::math::var;
  using stan::math::vector_v;
  using model_t = PMXTwoCptModel<double, double, double, double>;

  const int ncmt = model_t::Ncmt;
  torsten::PKRec<double> init(ncmt);
  init << 0, 100, 0;

  double t = time[1];
  std::vector<double> rate{1000., 0., 0.};

  model_t model(t, init, rate, pMatrix[0]);
  
  var a = 1500, ii = 18.0;
  double r = 100;
  int cmt = 1;

  // rate in model constructor should not affect steady-state solution
  std::vector<var> vars(model.vars(a, r, ii));
  EXPECT_EQ(vars.size(), 2);

  vector_v sol1 = model.solve(a, r, ii, cmt);
  Eigen::VectorXd sol2_d = model.solve_d(a, r, ii, cmt);
  vector_v sol2 = torsten::mpi::precomputed_gradients(sol2_d, vars);
  
  torsten::test::test_grad(vars, sol1, sol2, 1.E-8, 1.E-5);
}

TEST_F(TorstenTwoCptTest, rk45_model_solve_d_data_only) {
  

  
  using torsten::PMXTwoCptODE;
  using torsten::PMXOdeIntegrator;

  const int ncmt = PMXTwoCptModel<double, double, double, double>::Ncmt;
  PMXTwoCptODE f;

  using model_t = PKODEModel<double, double, double, double, PMXTwoCptODE>;
  PMXOdeIntegrator<torsten::StanRk45> integrator;

  torsten::PKRec<double> init(ncmt);
  init << 0, 100, 0;

  double t = time[1], t1 = t + 0.1;
  std::vector<double> rate{1000., 0., 0.};

  model_t model(t, init, rate, pMatrix[0], f);
  
  Eigen::VectorXd sol1 = model.solve(t1, integrator);
  Eigen::VectorXd sol2 = model.solve_d(t1, integrator);

  torsten::test::test_val(sol1, sol2);
}

TEST_F(TorstenTwoCptTest, rk45_model_solve_d_init_var) {
  
  
  using torsten::PMXTwoCptODE;
  using torsten::PMXOdeIntegrator;
  using stan::math::var;
  using stan::math::vector_v;

  const int ncmt = PMXTwoCptModel<double, double, double, double>::Ncmt;
  PMXTwoCptODE f;

  using model_t = PKODEModel<double, var, double, double, PMXTwoCptODE>;
  PMXOdeIntegrator<torsten::StanRk45> integrator;

  torsten::PKRec<var> init(ncmt);
  init << 0, 100, 0;

  double t = time[1], t1 = t + 0.1;
  std::vector<double> rate{1000., 0., 0.};

  model_t model(t, init, rate, pMatrix[0], f);
  std::vector<var> vars(model.vars(t1));
  EXPECT_EQ(vars.size(), init.size());
  
  vector_v sol1 = model.solve(t1, integrator);
  Eigen::VectorXd sol2_d = model.solve_d(t1, integrator);
  vector_v sol2 = torsten::mpi::precomputed_gradients(sol2_d, vars);
  
  // vars and init should be pointing to the same @c vari
  torsten::test::test_grad(vars, init, sol1, sol2, 1.E-8, 1.E-5);
}

TEST_F(TorstenTwoCptTest, rk45_model_solve_d_rate_var) {
  
  using torsten::PMXTwoCptODE;
  using torsten::PMXOdeIntegrator;
  using stan::math::var;
  using stan::math::vector_v;

  const int ncmt = PMXTwoCptModel<double, double, double, double>::Ncmt;
  PMXTwoCptODE f;

  using model_t = PKODEModel<double, double, var, double, PMXTwoCptODE>;
  PMXOdeIntegrator<torsten::StanRk45> integrator;

  torsten::PKRec<double> init(ncmt);
  init << 0, 100, 0;

  double t = time[1], t1 = t + 0.1;
  std::vector<var> rate{1000., 0., 0.};

  std::vector<double>& par = pMatrix[0];
  model_t model(t, init, rate, par, f);
  std::vector<var> vars(model.vars(t1));
  std::vector<var> par_rate;
  par_rate.insert(par_rate.end(), par.begin(), par.end());
  par_rate.insert(par_rate.end(), rate.begin(), rate.end());
  EXPECT_EQ(vars.size(), par_rate.size());
  
  vector_v sol1 = model.solve(t1, integrator);
  Eigen::VectorXd sol2_d = model.solve_d(t1, integrator);
  vector_v sol2 = torsten::mpi::precomputed_gradients(sol2_d, vars);
  
  // vars and init should be pointing to the same @c vari
  torsten::test::test_grad(vars, par_rate, sol1, sol2, 1.E-8, 1.E-5);
}

TEST_F(TorstenTwoCptTest, rk45_model_solve_d_par_var) {
  
  using torsten::PMXTwoCptODE;
  using torsten::PMXOdeIntegrator;
  using stan::math::var;
  using stan::math::vector_v;

  const int ncmt = PMXTwoCptModel<double, double, double, double>::Ncmt;
  PMXTwoCptODE f;

  using model_t = PKODEModel<double, double, double, var, PMXTwoCptODE>;
  PMXOdeIntegrator<torsten::StanRk45> integrator;

  torsten::PKRec<double> init(ncmt);
  init << 0, 100, 0;

  double t = time[1], t1 = t + 0.1;
  std::vector<double> rate{1000., 0., 0.};

  std::vector<var> par(stan::math::to_var(pMatrix[0]));
  model_t model(t, init, rate, par, f);
  std::vector<var> vars(model.vars(t1));
  EXPECT_EQ(vars.size(), par.size());
  
  vector_v sol1 = model.solve(t1, integrator);
  Eigen::VectorXd sol2_d = model.solve_d(t1, integrator);
  vector_v sol2 = torsten::mpi::precomputed_gradients(sol2_d, vars);
  
  // vars and init should be pointing to the same @c vari
  torsten::test::test_grad(vars, par, sol1, sol2, 1.E-8, 1.E-5);
}

TEST_F(TorstenTwoCptTest, rk45_model_solve_d_data_only_ss) {
  
  using torsten::PMXTwoCptODE;
  using torsten::PMXOdeIntegrator;

  const int ncmt = PMXTwoCptModel<double, double, double, double>::Ncmt;
  PMXTwoCptODE f;

  using model_t = PKODEModel<double, double, double, double, PMXTwoCptODE>;
  PMXOdeIntegrator<torsten::StanRk45> integrator;

  torsten::PKRec<double> init(ncmt);
  init << 0, 100, 0;

  double t = time[1];
  std::vector<double> rate{1000., 0., 0.};

  double a = 1500, r = 100, ii = 18.0;
  int cmt = 1;

  model_t model(t, init, rate, pMatrix[0], f);
  
  Eigen::VectorXd sol1 = model.solve(a, r, ii, cmt, integrator);
  Eigen::VectorXd sol2 = model.solve_d(a, r, ii, cmt, integrator);
  
  torsten::test::test_val(sol1, sol2);
}

TEST_F(TorstenTwoCptTest, rk45_model_solve_d_init_var_ss) {
  
  using torsten::PMXTwoCptODE;
  using torsten::PMXOdeIntegrator;
  using stan::math::var;
  using stan::math::vector_v;

  const int ncmt = PMXTwoCptModel<double, double, double, double>::Ncmt;
  PMXTwoCptODE f;

  using model_t = PKODEModel<double, var, double, double, PMXTwoCptODE>;
  PMXOdeIntegrator<torsten::StanRk45> integrator;

  torsten::PKRec<var> init(ncmt);
  init << 0, 100, 0;

  double t = time[1];
  std::vector<double> rate{1000., 0., 0.};

  double a = 1500, r = 100, ii = 18.0;
  int cmt = 1;

  model_t model(t, init, rate, pMatrix[0], f);
  std::vector<var> vars(model.vars(a, r, ii));

  // type of @c init should not affect steady state solution type.
  EXPECT_EQ(vars.size(), 0);
  
  Eigen::VectorXd sol1 = model.solve(a, r, ii, cmt, integrator);
  Eigen::VectorXd sol2 = model.solve_d(a, r, ii, cmt, integrator);
  
  torsten::test::test_val(sol1, sol2);
}

TEST_F(TorstenTwoCptTest, rk45_model_solve_d_par_var_ss) {
  
  using torsten::PMXTwoCptODE;
  using torsten::PMXOdeIntegrator;
  using stan::math::var;
  using stan::math::vector_v;

  const int ncmt = PMXTwoCptModel<double, double, double, double>::Ncmt;
  PMXTwoCptODE f;

  using model_t = PKODEModel<double, double, double, var, PMXTwoCptODE>;
  PMXOdeIntegrator<torsten::StanRk45> integrator;

  torsten::PKRec<double> init(ncmt);
  init << 0, 100, 0;

  double t = time[1];
  std::vector<double> rate{1000., 3000., 0.};

  double a = 1500, r = 100, ii = 18.0;
  int cmt = 1;

  std::vector<var> par(stan::math::to_var(pMatrix[0]));
  model_t model(t, init, rate, par, f);
  std::vector<var> vars(model.vars(a, r, ii));

  EXPECT_EQ(vars.size(), par.size());
  
  vector_v sol1 = model.solve(a, r, ii, cmt, integrator);
  Eigen::VectorXd sol2_d = model.solve_d(a, r, ii, cmt, integrator);
  vector_v sol2 = torsten::mpi::precomputed_gradients(sol2_d, vars);
  
  // vars and init should be pointing to the same @c vari
  torsten::test::test_grad(vars, par, sol1, sol2, 1.E-8, 1.E-5);
}

// PkBdf
TEST_F(TorstenTwoCptTest, PkBdf_model_solve_d_data_only) {
  
  using torsten::PMXTwoCptODE;
  using torsten::PMXOdeIntegrator;

  const int ncmt = PMXTwoCptModel<double, double, double, double>::Ncmt;
  PMXTwoCptODE f;

  using model_t = PKODEModel<double, double, double, double, PMXTwoCptODE>;
  PMXOdeIntegrator<torsten::PkBdf> integrator;

  torsten::PKRec<double> init(ncmt);
  init << 0, 100, 0;

  double t = time[1], t1 = t + 0.1;
  std::vector<double> rate{1000., 0., 0.};

  model_t model(t, init, rate, pMatrix[0], f);
  
  Eigen::VectorXd sol1 = model.solve(t1, integrator);
  Eigen::VectorXd sol2 = model.solve_d(t1, integrator);

  torsten::test::test_val(sol1, sol2);
}

TEST_F(TorstenTwoCptTest, PkBdf_model_solve_d_init_var) {
  
  using torsten::PMXTwoCptODE;
  using torsten::PMXOdeIntegrator;
  using stan::math::var;
  using stan::math::vector_v;

  const int ncmt = PMXTwoCptModel<double, double, double, double>::Ncmt;
  PMXTwoCptODE f;

  using model_t = PKODEModel<double, var, double, double, PMXTwoCptODE>;
  PMXOdeIntegrator<torsten::PkBdf> integrator;

  torsten::PKRec<var> init(ncmt);
  init << 0, 100, 0;

  double t = time[1], t1 = t + 0.1;
  std::vector<double> rate{1000., 0., 0.};

  model_t model(t, init, rate, pMatrix[0], f);
  std::vector<var> vars(model.vars(t1));
  EXPECT_EQ(vars.size(), init.size());
  
  vector_v sol1 = model.solve(t1, integrator);
  Eigen::VectorXd sol2_d = model.solve_d(t1, integrator);
  vector_v sol2 = torsten::mpi::precomputed_gradients(sol2_d, vars);
  
  // vars and init should be pointing to the same @c vari
  torsten::test::test_grad(vars, init, sol1, sol2, 1.E-8, 1.E-5);
}

TEST_F(TorstenTwoCptTest, PkBdf_model_solve_d_rate_var) {
  
  using torsten::PMXTwoCptODE;
  using torsten::PMXOdeIntegrator;
  using stan::math::var;
  using stan::math::vector_v;

  const int ncmt = PMXTwoCptModel<double, double, double, double>::Ncmt;
  PMXTwoCptODE f;

  using model_t = PKODEModel<double, double, var, double, PMXTwoCptODE>;
  PMXOdeIntegrator<torsten::PkBdf> integrator;

  torsten::PKRec<double> init(ncmt);
  init << 0, 100, 0;

  double t = time[1], t1 = t + 0.1;
  std::vector<var> rate{1000., 0., 0.};

  std::vector<double>& par = pMatrix[0];
  model_t model(t, init, rate, par, f);
  std::vector<var> vars(model.vars(t1));
  std::vector<var> par_rate;
  par_rate.insert(par_rate.end(), par.begin(), par.end());
  par_rate.insert(par_rate.end(), rate.begin(), rate.end());
  EXPECT_EQ(vars.size(), par_rate.size());
  
  vector_v sol1 = model.solve(t1, integrator);
  Eigen::VectorXd sol2_d = model.solve_d(t1, integrator);
  vector_v sol2 = torsten::mpi::precomputed_gradients(sol2_d, vars);
  
  // vars and init should be pointing to the same @c vari
  torsten::test::test_grad(vars, par_rate, sol1, sol2, 1.E-8, 1.E-5);
}

TEST_F(TorstenTwoCptTest, PkBdf_model_solve_d_par_var) {
  
  using torsten::PMXTwoCptODE;
  using torsten::PMXOdeIntegrator;
  using stan::math::var;
  using stan::math::vector_v;

  const int ncmt = PMXTwoCptModel<double, double, double, double>::Ncmt;
  PMXTwoCptODE f;

  using model_t = PKODEModel<double, double, double, var, PMXTwoCptODE>;
  PMXOdeIntegrator<torsten::PkBdf> integrator;

  torsten::PKRec<double> init(ncmt);
  init << 0, 100, 0;

  double t = time[1], t1 = t + 0.1;
  std::vector<double> rate{1000., 0., 0.};

  std::vector<var> par(stan::math::to_var(pMatrix[0]));
  model_t model(t, init, rate, par, f);
  std::vector<var> vars(model.vars(t1));
  EXPECT_EQ(vars.size(), par.size());
  
  vector_v sol1 = model.solve(t1, integrator);
  Eigen::VectorXd sol2_d = model.solve_d(t1, integrator);
  vector_v sol2 = torsten::mpi::precomputed_gradients(sol2_d, vars);
  
  // vars and init should be pointing to the same @c vari
  torsten::test::test_grad(vars, par, sol1, sol2, 1.E-8, 1.E-5);
}

TEST_F(TorstenTwoCptTest, PkBdf_model_solve_d_data_only_ss) {
  
  using torsten::PMXTwoCptODE;
  using torsten::PMXOdeIntegrator;

  const int ncmt = PMXTwoCptModel<double, double, double, double>::Ncmt;
  PMXTwoCptODE f;

  using model_t = PKODEModel<double, double, double, double, PMXTwoCptODE>;
  PMXOdeIntegrator<torsten::PkBdf> integrator;

  torsten::PKRec<double> init(ncmt);
  init << 0, 100, 0;

  double t = time[1];
  std::vector<double> rate{1000., 0., 0.};

  double a = 1500, r = 100, ii = 18.0;
  int cmt = 1;

  model_t model(t, init, rate, pMatrix[0], f);
  
  Eigen::VectorXd sol1 = model.solve(a, r, ii, cmt, integrator);
  Eigen::VectorXd sol2 = model.solve_d(a, r, ii, cmt, integrator);
  
  torsten::test::test_val(sol1, sol2);
}

TEST_F(TorstenTwoCptTest, PkBdf_model_solve_d_init_var_ss) {
  
  using torsten::PMXTwoCptODE;
  using torsten::PMXOdeIntegrator;
  using stan::math::var;
  using stan::math::vector_v;

  const int ncmt = PMXTwoCptModel<double, double, double, double>::Ncmt;
  PMXTwoCptODE f;

  using model_t = PKODEModel<double, var, double, double, PMXTwoCptODE>;
  PMXOdeIntegrator<torsten::PkBdf> integrator;

  torsten::PKRec<var> init(ncmt);
  init << 0, 100, 0;

  double t = time[1];
  std::vector<double> rate{1000., 0., 0.};

  double a = 1500, r = 100, ii = 18.0;
  int cmt = 1;

  model_t model(t, init, rate, pMatrix[0], f);
  std::vector<var> vars(model.vars(a, r, ii));

  // type of @c init should not affect steady state solution type.
  EXPECT_EQ(vars.size(), 0);
  
  Eigen::VectorXd sol1 = model.solve(a, r, ii, cmt, integrator);
  Eigen::VectorXd sol2 = model.solve_d(a, r, ii, cmt, integrator);
  
  torsten::test::test_val(sol1, sol2);
}

TEST_F(TorstenTwoCptTest, PkBdf_model_solve_d_par_var_ss) {
  
  using torsten::PMXTwoCptODE;
  using torsten::PMXOdeIntegrator;
  using stan::math::var;
  using stan::math::vector_v;

  const int ncmt = PMXTwoCptModel<double, double, double, double>::Ncmt;
  PMXTwoCptODE f;

  using model_t = PKODEModel<double, double, double, var, PMXTwoCptODE>;
  PMXOdeIntegrator<torsten::PkBdf> integrator;

  torsten::PKRec<double> init(ncmt);
  init << 0, 100, 0;

  double t = time[1];
  std::vector<double> rate{1000., 3000., 0.};

  double a = 1500, r = 100, ii = 18.0;
  int cmt = 1;

  std::vector<var> par(stan::math::to_var(pMatrix[0]));
  model_t model(t, init, rate, par, f);
  std::vector<var> vars(model.vars(a, r, ii));

  EXPECT_EQ(vars.size(), par.size());
  
  vector_v sol1 = model.solve(a, r, ii, cmt, integrator);
  Eigen::VectorXd sol2_d = model.solve_d(a, r, ii, cmt, integrator);
  vector_v sol2 = torsten::mpi::precomputed_gradients(sol2_d, vars);
  
  // vars and init should be pointing to the same @c vari
  torsten::test::test_grad(vars, par, sol1, sol2, 1.E-8, 1.E-5);
}
