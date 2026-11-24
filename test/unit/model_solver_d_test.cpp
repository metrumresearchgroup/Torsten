#include <stan/math/rev/core/typedefs.hpp>
#include <stan/math/torsten/ev_manager.hpp>
#include <stan/math/torsten/torsten_def.hpp>
#include <stan/math/torsten/dsolve/pk_vars.hpp>
#include <stan/math/torsten/pmx_twocpt_model.hpp>
#include <stan/math/torsten/model_solve_d.hpp>
#include <stan/math/torsten/pmx_ode_model.hpp>
#include <stan/math/torsten/test/unit/pmx_twocpt_test_fixture.hpp>
#include <stan/math/torsten/test/unit/test_macros.hpp>
#include <gtest/gtest.h>
#include <vector>

using std::vector;
using Eigen::Matrix;
using Eigen::Dynamic;
using stan::math::var;
using stan::math::to_var;
using stan::math::vector_v;
using torsten::PKRec;
using torsten::PMXTwoCptModel;
using torsten::PKODEModel;
using torsten::PMXTwoCptODE;
using torsten::dsolve::pk_vars;
using torsten::pmx_model_vars;
using torsten::dsolve::PMXOdeIntegrator;
using torsten::dsolve::PMXVariadicOdeSystem;
using torsten::dsolve::PMXCvodesIntegrator;
using torsten::dsolve::PMXOdeintIntegrator;

using integ_t = torsten::dsolve::PMXAnalyiticalIntegrator;
using torsten::dsolve::odeint_scheme_rk45;

TEST_F(TorstenTwoCptTest, model_solve_d_data_only) {
  using model_t = PMXTwoCptModel<double>;

  const int ncmt = model_t::Ncmt;
  torsten::PKRec<double> init(ncmt);
  init << 0, 100, 0;

  double t = time[1], t1 = t + 0.1;
  std::vector<double> rate{1000., 0., 0.};

  model_t model(pMatrix[0]);
  
  torsten::PKRec<double> y1(init), y2(init);
  model.solve(y1, t, t1, rate);
  y2 = model_solve_d(model, init, t, t1, rate, integ_t());
  
  torsten::test::test_val(y1, y2);
}

TEST_F(TorstenTwoCptTest, model_solve_d_init_var) {
  using model_t = PMXTwoCptModel<double>;

  const int ncmt = model_t::Ncmt;
  torsten::PKRec<var> init(ncmt);
  init << 200, 100, 0;

  double t = time[1], t1 = t + 0.1;
  std::vector<double> rate{400., 1000., 0.};

  model_t model(pMatrix[0]);

  std::vector<var> vars(pk_vars(t1, init, rate, pMatrix[0]));
  EXPECT_EQ(vars.size(), init.size());

  vector_v sol1(to_var(init));
  model.solve(sol1, t, t1, rate);
  Eigen::VectorXd sol2_d;
  sol2_d = model_solve_d(model, init, t, t1, rate, integ_t());
  vector_v sol2 = torsten::mpi::precomputed_gradients(sol2_d, vars);
  torsten::test::test_grad(vars, sol1, sol2, 1.E-8, 1.E-5);
}

TEST_F(TorstenTwoCptTest, model_solve_d_rate_var) {
  using model_t = PMXTwoCptModel<double>;

  const int ncmt = model_t::Ncmt;
  torsten::PKRec<double> init(ncmt);
  init << 200, 100, 0;

  double t = time[1], t1 = t + 0.1;
  std::vector<var> rate{400., 1000., 0.};

  model_t model(pMatrix[0]);

  torsten::PKRec<var> y0(to_var(init));
  std::vector<var> vars(pk_vars(t1, init, rate, pMatrix[0]));
  EXPECT_EQ(vars.size(), rate.size());

  torsten::PKRec<var> sol1(y0);
  model.solve(sol1, t, t1, rate);
  Eigen::VectorXd sol2_d;
  sol2_d = model_solve_d(model, init, t, t1, rate, integ_t());
  vector_v sol2 = torsten::mpi::precomputed_gradients(sol2_d, vars);
  torsten::test::test_grad(vars, sol1, sol2, 1.E-8, 1.E-5);
}

TEST_F(TorstenTwoCptTest, model_solve_d_par_var) {
  using model_t = PMXTwoCptModel<var>;

  const int ncmt = model_t::Ncmt;
  torsten::PKRec<double> init(ncmt);
  init << 200, 100, 0;

  double t = time[1], t1 = t + 0.1;
  std::vector<double> rate{400., 1000., 0.};

  std::vector<var> pars(stan::math::to_var(pMatrix[0]));
  model_t model(pars);

  torsten::PKRec<var> y0(to_var(init));
  std::vector<var> vars(pk_vars(t1, init, rate, pars));
  EXPECT_EQ(vars.size(), pars.size());

  torsten::PKRec<var> sol1(y0);
  model.solve(sol1, t, t1, rate);
  Eigen::VectorXd sol2_d;
  sol2_d = model_solve_d(model, init, t, t1, rate, integ_t());
  vector_v sol2 = torsten::mpi::precomputed_gradients(sol2_d, vars);
  torsten::test::test_grad(vars, sol1, sol2, 1.E-8, 1.E-5);
}

// TEST_F(TorstenTwoCptTest, model_solve_d_rate_par_var) {
//   using model_t = PMXTwoCptModel<var>;

//   const int ncmt = model_t::Ncmt;
//   torsten::PKRec<double> init(ncmt);
//   init << 200, 100, 0;

//   double t = time[1], t1 = t + 0.1;
//   std::vector<var> rate{400., 1000., 0.};

//   std::vector<var> pars(stan::math::to_var(pMatrix[0]));
//   model_t model(pars);

//   std::vector<var> vars(model.vars(t1));
//   EXPECT_EQ(vars.size(), pars.size() + rate.size());

//   vector_v sol1(to_var(init));
//   model.solve(sol1, t, t1, rate, integ_t());
//   Eigen::VectorXd sol2_d = model_solve_d(model, init, t, t1, rate, integ_t());
//   vector_v sol2 = torsten::mpi::precomputed_gradients(sol2_d, vars);
//   torsten::test::test_grad(vars, sol1, sol2, 1.E-8, 1.E-5);
// }

TEST_F(TorstenTwoCptTest, model_solve_d_data_only_ss) {
  using model_t = PMXTwoCptModel<double>;
  const int ncmt = model_t::Ncmt;
  torsten::PKRec<double> init(ncmt);
  init << 0, 100, 0;

  double t = time[1];
  std::vector<double> rate{1000., 0., 0.};

  model_t model(pMatrix[0]);
  
  double a = 1500, r = 100, ii = 18.0;
  int cmt = 1;
  Eigen::VectorXd sol1 = model.solve(t, a, r, ii, cmt);
  Eigen::VectorXd sol2 = model_solve_d(model, t, a, r, ii, cmt, integ_t());
  torsten::test::test_val(sol1, sol2);
}

TEST_F(TorstenTwoCptTest, model_solve_d_init_var_ss) {
  using model_t = PMXTwoCptModel<double>;

  const int ncmt = model_t::Ncmt;
  torsten::PKRec<var> init(ncmt);
  init << 0, 100, 0;

  double t = time[1];
  std::vector<double> rate{1000., 0., 0.};

  model_t model(pMatrix[0]);
  
  double a = 1500, r = 100, ii = 18.0;
  int cmt = 1;

  // initial condition should not affect steady-state solution
  std::vector<var> vars(pk_vars(a, r, ii, model.par()));
  EXPECT_EQ(vars.size(), 0);

  Eigen::VectorXd sol1 = model.solve(t, a, r, ii, cmt);
  Eigen::VectorXd sol2 = model_solve_d(model, t, a, r, ii, cmt, integ_t());
  EXPECT_EQ(sol2.size(), sol1.size());
  torsten::test::test_val(sol1, sol2);
}

TEST_F(TorstenTwoCptTest, model_solve_d_rate_var_ss) {
  using model_t = PMXTwoCptModel<double>;

  const int ncmt = model_t::Ncmt;
  torsten::PKRec<double> init(ncmt);
  init << 0, 100, 0;

  double t = time[1];
  std::vector<var> rate{1000., 0., 0.};

  model_t model(pMatrix[0]);
  
  double a = 1500, r = 100, ii = 18.0;
  int cmt = 1;

  // rate in model constructor should not affect steady-state solution
  std::vector<var> vars(pk_vars(a, r, ii, model.par()));
  EXPECT_EQ(vars.size(), 0);

  Eigen::VectorXd sol1 = model.solve(t, a, r, ii, cmt);
  Eigen::VectorXd sol2 = model_solve_d(model, t, a, r, ii, cmt, integ_t());
  EXPECT_EQ(sol2.size(), sol1.size());
  torsten::test::test_val(sol1, sol2);
}

TEST_F(TorstenTwoCptTest, model_solve_d_par_var_ss) {
  
  using stan::math::var;
  using stan::math::vector_v;
  using model_t = PMXTwoCptModel<var>;

  const int ncmt = model_t::Ncmt;
  torsten::PKRec<double> init(ncmt);
  init << 0, 100, 0;

  double t = time[1];
  std::vector<double> rate{1000., 0., 0.};

  std::vector<var> pars(stan::math::to_var(pMatrix[0]));

  model_t model(pars);
  
  double a = 1500, r = 100, ii = 18.0;
  int cmt = 1;

  // rate in model constructor should not affect steady-state solution
  std::vector<var> vars(pk_vars(a, r, ii, model.par()));
  EXPECT_EQ(vars.size(), pars.size());

  vector_v sol1 = model.solve(t, a, r, ii, cmt);
  Eigen::VectorXd sol2_d = model_solve_d(model, t, a, r, ii, cmt, integ_t());
  vector_v sol2 = torsten::mpi::precomputed_gradients(sol2_d, vars);
  
  torsten::test::test_grad(vars, pars, sol1, sol2, 1.E-8, 1.E-5);
}

TEST_F(TorstenTwoCptTest, model_solve_d_amt_par_var_ss) {
  using model_t = PMXTwoCptModel<var>;

  const int ncmt = model_t::Ncmt;
  torsten::PKRec<double> init(ncmt);
  init << 0, 100, 0;

  double t = time[1];
  std::vector<double> rate{1000., 0., 0.};

  std::vector<var> pars(stan::math::to_var(pMatrix[0]));

  model_t model(pars);
  
  var a = 1500;
  double r = 100, ii = 18.0;
  int cmt = 1;

  // rate in model constructor should not affect steady-state solution
  std::vector<var> vars(pk_vars(a, r, ii, model.par()));
  EXPECT_EQ(vars.size(), pars.size() + 1);

  vector_v sol1 = model.solve(t, a, r, ii, cmt);
  Eigen::VectorXd sol2_d = model_solve_d(model, t, a, r, ii, cmt, integ_t());
  vector_v sol2 = torsten::mpi::precomputed_gradients(sol2_d, vars);
  
  torsten::test::test_grad(vars, sol1, sol2, 1.E-8, 1.E-5);
}

TEST_F(TorstenTwoCptTest, model_solve_d_amt_ii_var_ss) {
  
  using stan::math::var;
  using stan::math::vector_v;
  using model_t = PMXTwoCptModel<double>;

  const int ncmt = model_t::Ncmt;
  torsten::PKRec<double> init(ncmt);
  init << 0, 100, 0;

  double t = time[1];
  std::vector<double> rate{1000., 0., 0.};

  model_t model(pMatrix[0]);
  
  var a = 1500, ii = 18.0;
  double r = 100;
  int cmt = 1;

  // rate in model constructor should not affect steady-state solution
  std::vector<var> vars(pk_vars(a, r, ii, model.par()));
  EXPECT_EQ(vars.size(), 2);

  vector_v sol1 = model.solve(t, a, r, ii, cmt);
  Eigen::VectorXd sol2_d = model_solve_d(model, t, a, r, ii, cmt, integ_t());
  vector_v sol2 = torsten::mpi::precomputed_gradients(sol2_d, vars);
  
  torsten::test::test_grad(vars, sol1, sol2, 1.E-8, 1.E-5);
}

TEST_F(TorstenTwoCptTest, rk45_model_solve_d_data_only) {
  const int ncmt = PMXTwoCptModel<double>::Ncmt;
  PMXTwoCptODE f;

  using model_t = PKODEModel<double, PMXTwoCptODE>;
  PMXOdeIntegrator<PMXVariadicOdeSystem, PMXOdeintIntegrator<odeint_scheme_rk45>> integrator;

  torsten::PKRec<double> init(ncmt);
  init << 0, 100, 0;

  double t = time[1], t1 = t + 0.1;
  std::vector<double> rate{1000., 0., 0.};

  model_t model(pMatrix[0], ncmt, f);
  
  Eigen::VectorXd sol1(init);
  model.solve(sol1, t, t1, rate, integrator);
  Eigen::VectorXd sol2 = model_solve_d(model, init, t, t1, rate, integrator, ncmt, f);

  torsten::test::test_val(sol1, sol2);
}

TEST_F(TorstenTwoCptTest, rk45_model_solve_d_init_var) {
  const int ncmt = PMXTwoCptModel<double>::Ncmt;
  PMXTwoCptODE f;

  using model_t = PKODEModel<double, PMXTwoCptODE>;
  PMXOdeIntegrator<PMXVariadicOdeSystem, PMXOdeintIntegrator<odeint_scheme_rk45>> integrator;

  torsten::PKRec<var> init(ncmt);
  init << 0, 100, 0;

  double t = time[1], t1 = t + 0.1;
  std::vector<double> rate{1000., 0., 0.};

  std::vector<double>& par = pMatrix[0];
  model_t model(pMatrix[0], ncmt, f);
  std::vector<var> vars(pmx_model_vars<model_t>::vars(t1, init, rate, par));
  EXPECT_EQ(vars.size(), init.size());
  
  vector_v sol1(init);
  model.solve(sol1, t, t1, rate, integrator);
  Eigen::VectorXd sol2_d = model_solve_d(model, init, t, t1, rate, integrator, ncmt, f);
  vector_v sol2 = torsten::mpi::precomputed_gradients(sol2_d, vars);
  
  // vars and init should be pointing to the same @c vari
  torsten::test::test_grad(vars, init, sol1, sol2, 1.E-8, 1.E-5);
}

TEST_F(TorstenTwoCptTest, rk45_model_solve_d_rate_var) {
  const int ncmt = PMXTwoCptModel<double>::Ncmt;
  PMXTwoCptODE f;

  using model_t = PKODEModel<double, PMXTwoCptODE>;
  PMXOdeIntegrator<PMXVariadicOdeSystem, PMXOdeintIntegrator<odeint_scheme_rk45>> integrator;

  torsten::PKRec<double> init(ncmt);
  init << 0, 100, 0;

  double t = time[1], t1 = t + 0.1;
  std::vector<var> rate{1000., 0., 0.};

  std::vector<double>& par = pMatrix[0];
  model_t model(par, ncmt, f);
  std::vector<var> vars(pmx_model_vars<model_t>::vars(t1, init, rate, par));
  EXPECT_EQ(vars.size(), rate.size());
  
  vector_v sol1(to_var(init));
  model.solve(sol1, t, t1, rate, integrator);
  Eigen::VectorXd sol2_d = model_solve_d(model, init, t, t1, rate, integrator, ncmt, f);
  vector_v sol2 = torsten::mpi::precomputed_gradients(sol2_d, vars);
  
  // vars and init should be pointing to the same @c vari
  torsten::test::test_grad(vars, rate, sol1, sol2, 1.E-8, 1.E-5);
}

TEST_F(TorstenTwoCptTest, rk45_model_solve_d_par_var) {
  const int ncmt = PMXTwoCptModel<double>::Ncmt;
  PMXTwoCptODE f;

  using model_t = PKODEModel<var, PMXTwoCptODE>;
  PMXOdeIntegrator<PMXVariadicOdeSystem, PMXOdeintIntegrator<odeint_scheme_rk45>> integrator;

  torsten::PKRec<double> init(ncmt);
  init << 0, 100, 0;

  double t = time[1], t1 = t + 0.1;
  std::vector<double> rate{1000., 0., 0.};

  std::vector<var> par(stan::math::to_var(pMatrix[0]));
  model_t model(par, ncmt, f);
  std::vector<var> vars(pmx_model_vars<model_t>::vars(t1, init, rate, par));
  EXPECT_EQ(vars.size(), par.size());
  
  vector_v sol1(to_var(init));
  model.solve(sol1, t, t1, rate, integrator);
  Eigen::VectorXd sol2_d = model_solve_d(model, init, t, t1, rate, integrator, ncmt, f);
  vector_v sol2 = torsten::mpi::precomputed_gradients(sol2_d, vars);
  
  // vars and init should be pointing to the same @c vari
  torsten::test::test_grad(vars, par, sol1, sol2, 1.E-8, 1.E-5);
}

TEST_F(TorstenTwoCptTest, rk45_model_solve_d_data_only_ss) {
  const int ncmt = PMXTwoCptModel<double>::Ncmt;
  PMXTwoCptODE f;

  using model_t = PKODEModel<double, PMXTwoCptODE>;
  PMXOdeIntegrator<PMXVariadicOdeSystem, PMXOdeintIntegrator<odeint_scheme_rk45>> integrator;

  torsten::PKRec<double> init(ncmt);
  init << 0, 100, 0;

  double t = time[1];
  std::vector<double> rate{1000., 0., 0.};

  double a = 1500, r = 100, ii = 18.0;
  int cmt = 1;

  model_t model(pMatrix[0], ncmt, f);
  
  Eigen::VectorXd sol1 = model.solve(t, a, r, ii, cmt, integrator);
  Eigen::VectorXd sol2 = model_solve_d(model, t, a, r, ii, cmt, integrator, ncmt, f);
  
  torsten::test::test_val(sol1, sol2);
}

TEST_F(TorstenTwoCptTest, rk45_model_solve_d_init_var_ss) {
  const int ncmt = PMXTwoCptModel<double>::Ncmt;
  PMXTwoCptODE f;

  using model_t = PKODEModel<double, PMXTwoCptODE>;
  PMXOdeIntegrator<PMXVariadicOdeSystem, PMXOdeintIntegrator<odeint_scheme_rk45>> integrator;

  torsten::PKRec<var> init(ncmt);
  init << 0, 100, 0;

  double t = time[1];
  std::vector<double> rate{1000., 0., 0.};

  double a = 1500, r = 100, ii = 18.0;
  int cmt = 1;

  model_t model(pMatrix[0], ncmt, f);
  std::vector<var> vars(pk_vars(a, r, ii, model.par()));

  // type of @c init should not affect steady state solution type.
  EXPECT_EQ(vars.size(), 0);
  
  Eigen::VectorXd sol1 = model.solve(t, a, r, ii, cmt, integrator);
  Eigen::VectorXd sol2 = model_solve_d(model, t, a, r, ii, cmt, integrator, ncmt, f);
  
  torsten::test::test_val(sol1, sol2);
}

TEST_F(TorstenTwoCptTest, rk45_model_solve_d_par_var_ss) {
  
  using torsten::PMXTwoCptODE;
  using torsten::dsolve::PMXOdeIntegrator;
  using stan::math::var;
  using stan::math::vector_v;

  const int ncmt = PMXTwoCptModel<double>::Ncmt;
  PMXTwoCptODE f;

  using model_t = PKODEModel<var, PMXTwoCptODE>;
  PMXOdeIntegrator<PMXVariadicOdeSystem, PMXOdeintIntegrator<odeint_scheme_rk45>> integrator;

  torsten::PKRec<double> init(ncmt);
  init << 0, 100, 0;

  double t = time[1];
  std::vector<double> rate{1000., 3000., 0.};

  double a = 1500, r = 100, ii = 18.0;
  int cmt = 1;

  std::vector<var> par(stan::math::to_var(pMatrix[0]));
  model_t model(par, ncmt, f);
  std::vector<var> vars(pk_vars(a, r, ii, model.par()));

  EXPECT_EQ(vars.size(), par.size());
  
  vector_v sol1 = model.solve(t, a, r, ii, cmt, integrator);
  Eigen::VectorXd sol2_d = model_solve_d(model, t, a, r, ii, cmt, integrator, ncmt, f);
  vector_v sol2 = torsten::mpi::precomputed_gradients(sol2_d, vars);
  
  // vars and init should be pointing to the same @c vari
  torsten::test::test_grad(vars, par, sol1, sol2, 1.E-8, 1.E-5);
}

// PkBdf
TEST_F(TorstenTwoCptTest, PkBdf_model_solve_d_data_only) {
  const int ncmt = PMXTwoCptModel<double>::Ncmt;
  PMXTwoCptODE f;

  using model_t = PKODEModel<double, PMXTwoCptODE>;
  PMXOdeIntegrator<PMXVariadicOdeSystem, PMXCvodesIntegrator<CV_BDF, CV_STAGGERED>> integrator;

  torsten::PKRec<double> init(ncmt);
  init << 0, 100, 0;

  double t = time[1], t1 = t + 0.1;
  std::vector<double> rate{1000., 0., 0.};

  model_t model(pMatrix[0], ncmt, f);
  
  Eigen::VectorXd sol1(init);
  model.solve(sol1, t, t1, rate, integrator);
  Eigen::VectorXd sol2 = model_solve_d(model, init, t, t1, rate, integrator, ncmt, f);

  torsten::test::test_val(sol1, sol2);
}

TEST_F(TorstenTwoCptTest, PkBdf_model_solve_d_init_var) {
  const int ncmt = PMXTwoCptModel<double>::Ncmt;
  PMXTwoCptODE f;

  using model_t = PKODEModel<double, PMXTwoCptODE>;
  PMXOdeIntegrator<PMXVariadicOdeSystem, PMXCvodesIntegrator<CV_BDF, CV_STAGGERED>> integrator;

  torsten::PKRec<var> init(ncmt);
  init << 0, 100, 0;

  double t = time[1], t1 = t + 0.1;
  std::vector<double> rate{1000., 0., 0.};

  std::vector<double>& par = pMatrix[0];
  model_t model(pMatrix[0], ncmt, f);
  std::vector<var> vars(pmx_model_vars<model_t>::vars(t1, init, rate, par));
  EXPECT_EQ(vars.size(), init.size());
  
  vector_v sol1(init);
  model.solve(sol1, t, t1, rate, integrator);
  Eigen::VectorXd sol2_d = model_solve_d(model, init, t, t1, rate, integrator, ncmt, f);
  vector_v sol2 = torsten::mpi::precomputed_gradients(sol2_d, vars);
  
  // vars and init should be pointing to the same @c vari
  torsten::test::test_grad(vars, init, sol1, sol2, 1.E-8, 1.E-5);
}

TEST_F(TorstenTwoCptTest, PkBdf_model_solve_d_rate_var) {
  const int ncmt = PMXTwoCptModel<double>::Ncmt;
  PMXTwoCptODE f;

  using model_t = PKODEModel<double, PMXTwoCptODE>;
  PMXOdeIntegrator<PMXVariadicOdeSystem, PMXCvodesIntegrator<CV_BDF, CV_STAGGERED>> integrator;

  torsten::PKRec<double> init(ncmt);
  init << 0, 100, 0;

  double t = time[1], t1 = t + 0.1;
  std::vector<var> rate{1000., 0., 0.};

  std::vector<double>& par = pMatrix[0];
  model_t model(par, ncmt, f);
  std::vector<var> vars(pmx_model_vars<model_t>::vars(t1, init, rate, par));
  EXPECT_EQ(vars.size(), rate.size());
  
  vector_v sol1(to_var(init));
  model.solve(sol1, t, t1, rate, integrator);
  Eigen::VectorXd sol2_d = model_solve_d(model, init, t, t1, rate, integrator, ncmt, f);
  vector_v sol2 = torsten::mpi::precomputed_gradients(sol2_d, vars);
  
  // vars and init should be pointing to the same @c vari
  torsten::test::test_grad(vars, rate, sol1, sol2, 1.E-8, 1.E-5);
}

TEST_F(TorstenTwoCptTest, PkBdf_model_solve_d_par_var) {
  
  using torsten::PMXTwoCptODE;
  using torsten::dsolve::PMXOdeIntegrator;
  using stan::math::var;
  using stan::math::vector_v;

  const int ncmt = PMXTwoCptModel<double>::Ncmt;
  PMXTwoCptODE f;

  using model_t = PKODEModel<var, PMXTwoCptODE>;
  PMXOdeIntegrator<PMXVariadicOdeSystem, PMXCvodesIntegrator<CV_BDF, CV_STAGGERED>> integrator;

  torsten::PKRec<double> init(ncmt);
  init << 0, 100, 0;

  double t = time[1], t1 = t + 0.1;
  std::vector<double> rate{1000., 0., 0.};

  std::vector<var> par(stan::math::to_var(pMatrix[0]));
  model_t model(par, ncmt, f);
  std::vector<var> vars(pmx_model_vars<model_t>::vars(t1, init, rate, par));
  EXPECT_EQ(vars.size(), par.size());
  
  vector_v sol1(to_var(init));
  model.solve(sol1, t, t1, rate, integrator);
  Eigen::VectorXd sol2_d = model_solve_d(model, init, t, t1, rate, integrator, ncmt, f);
  vector_v sol2 = torsten::mpi::precomputed_gradients(sol2_d, vars);
  
  // vars and init should be pointing to the same @c vari
  torsten::test::test_grad(vars, par, sol1, sol2, 1.E-8, 1.E-5);
}

TEST_F(TorstenTwoCptTest, PkBdf_model_solve_d_data_only_ss) {
  
  using torsten::PMXTwoCptODE;
  using torsten::dsolve::PMXOdeIntegrator;

  const int ncmt = PMXTwoCptModel<double>::Ncmt;
  PMXTwoCptODE f;

  using model_t = PKODEModel<double, PMXTwoCptODE>;
  PMXOdeIntegrator<PMXVariadicOdeSystem, PMXCvodesIntegrator<CV_BDF, CV_STAGGERED>> integrator;

  torsten::PKRec<double> init(ncmt);
  init << 0, 100, 0;

  double t = time[1];
  std::vector<double> rate{1000., 0., 0.};

  double a = 1500, r = 100, ii = 18.0;
  int cmt = 1;

  model_t model(pMatrix[0], ncmt, f);
  
  Eigen::VectorXd sol1 = model.solve(t, a, r, ii, cmt, integrator);
  Eigen::VectorXd sol2 = model_solve_d(model, t, a, r, ii, cmt, integrator, ncmt, f);
  
  torsten::test::test_val(sol1, sol2);
}

TEST_F(TorstenTwoCptTest, PkBdf_model_solve_d_init_var_ss) {
  
  using torsten::PMXTwoCptODE;
  using torsten::dsolve::PMXOdeIntegrator;
  using stan::math::var;
  using stan::math::vector_v;

  const int ncmt = PMXTwoCptModel<double>::Ncmt;
  PMXTwoCptODE f;

  using model_t = PKODEModel<double, PMXTwoCptODE>;
  PMXOdeIntegrator<PMXVariadicOdeSystem, PMXCvodesIntegrator<CV_BDF, CV_STAGGERED>> integrator;

  torsten::PKRec<var> init(ncmt);
  init << 0, 100, 0;

  double t = time[1];
  std::vector<double> rate{1000., 0., 0.};

  double a = 1500, r = 100, ii = 18.0;
  int cmt = 1;

  model_t model(pMatrix[0], ncmt, f);
  std::vector<var> vars(pk_vars(a, r, ii, model.par()));

  // type of @c init should not affect steady state solution type.
  EXPECT_EQ(vars.size(), 0);
  
  Eigen::VectorXd sol1 = model.solve(t, a, r, ii, cmt, integrator);
  Eigen::VectorXd sol2 = model_solve_d(model, t, a, r, ii, cmt, integrator, ncmt, f);
  
  torsten::test::test_val(sol1, sol2);
}

TEST_F(TorstenTwoCptTest, PkBdf_model_solve_d_par_var_ss) {
  const int ncmt = PMXTwoCptModel<double>::Ncmt;
  PMXTwoCptODE f;

  using model_t = PKODEModel<var, PMXTwoCptODE>;
  PMXOdeIntegrator<PMXVariadicOdeSystem, PMXCvodesIntegrator<CV_BDF, CV_STAGGERED>> integrator;

  torsten::PKRec<double> init(ncmt);
  init << 0, 100, 0;

  double t = time[1];
  std::vector<double> rate{1000., 3000., 0.};

  double a = 1500, r = 100, ii = 18.0;
  int cmt = 1;

  std::vector<var> par(stan::math::to_var(pMatrix[0]));
  model_t model(par, ncmt, f);
  std::vector<var> vars(pk_vars(a, r, ii, model.par()));

  EXPECT_EQ(vars.size(), par.size());
  
  vector_v sol1 = model.solve(t, a, r, ii, cmt, integrator);
  Eigen::VectorXd sol2_d = model_solve_d(model, t, a, r, ii, cmt, integrator, ncmt, f);
  vector_v sol2 = torsten::mpi::precomputed_gradients(sol2_d, vars);
  
  // vars and init should be pointing to the same @c vari
  torsten::test::test_grad(vars, par, sol1, sol2, 1.E-8, 1.E-5);
}
