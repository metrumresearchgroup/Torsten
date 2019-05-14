#include <stan/lang/ast_def.cpp>
#include <stan/lang/generator.hpp>
#include <gtest/gtest.h>
#include <boost/variant/polymorphic_get.hpp>
#include <cmath>

#include <fstream>
#include <iostream>
#include <sstream>
#include <string>
#include <set>
#include <vector>

using stan::lang::idx;
using stan::lang::uni_idx;
using stan::lang::omni_idx;
using stan::lang::expression;
using stan::lang::int_literal;
using stan::lang::function_signatures;
using stan::lang::function_arg_type;
using stan::lang::expr_type;
using stan::lang::base_expr_type;
using stan::lang::ill_formed_type;
using stan::lang::void_type;
using stan::lang::double_type;
using stan::lang::int_type;
using stan::lang::vector_type;
using stan::lang::row_vector_type;
using stan::lang::matrix_type;
using std::vector;

TEST(langAst, univariate_integral) {
  using stan::lang::univariate_integral_control;
  using stan::lang::variable;
  using stan::lang::expr_type;
  using stan::lang::expression;

  univariate_integral_control so;

  std::string integration_function_name = "bar";
  std::string system_function_name = "foo";

  variable y0("y0_var_name");
  y0.set_type(double_type(), 1);  // plain old vector

  variable t0("t0_var_name");
  t0.set_type(double_type(), 0);  // double

  variable t1("t0_var_name");
  t1.set_type(double_type(), 0);  // double

  variable theta("theta_var_name");
  theta.set_type(double_type(), 1);

  variable x("x_var_name");
  x.set_type(double_type(), 1);

  variable x_int("x_int_var_name");
  x.set_type(int_type(), 1);

  univariate_integral_control so2(integration_function_name,
                                  system_function_name,
                                  t0, t1, theta, x, x_int);

  // dumb test to make sure we at least get the right types back
  EXPECT_EQ(integration_function_name, so2.integration_function_name_);
  EXPECT_EQ(system_function_name, so2.system_function_name_);
  EXPECT_EQ(t0.type_, so2.t0_.expression_type());
  EXPECT_EQ(t1.type_, so2.t1_.expression_type());
  EXPECT_EQ(theta.type_, so2.theta_.expression_type());
  EXPECT_EQ(x.type_, so2.x_r_.expression_type());
  EXPECT_EQ(x_int.type_, so2.x_i_.expression_type());

  expression e2(so2);
  EXPECT_EQ(expr_type(double_type()), e2.expression_type());
}

TEST(langAst, generalOdeModel_control) {
  using stan::lang::generalOdeModel_control;
  using stan::lang::variable;
  using stan::lang::expr_type;
  using stan::lang::expression;

  generalOdeModel_control so;

  std::string integration_function_name = "bar";
  std::string system_function_name = "foo";

  variable nCmt("nCmt_var_name");
  nCmt.set_type(int_type(), 0);  // int

  variable time("time_var_name");
  time.set_type(double_type(), 1);  // time of events

  variable amt("amt_var_name");
  amt.set_type(double_type(), 1);  // Amount of events

  variable rate("rate_var_name");
  rate.set_type(double_type(), 1);  // rate at events

  variable ii("ii_var_name");
  ii.set_type(double_type(), 1);  // interdose interval

  variable evid("evid_var_name");
  evid.set_type(int_type(), 1);  // event id

  variable cmt("cmt_var_name");
  cmt.set_type(int_type(), 1);  // compartment nb. of events

  variable addl("addl_var_name");
  addl.set_type(int_type(), 1);  // Number of additional doses

  variable ss("ss_var_name");
  ss.set_type(int_type(), 1);  // steady state flag

  variable pmatrix("pmatrix_var_name");
  pmatrix.set_type(double_type(), 2);  // linear ode system matrix

  variable biovar("biovar_var_name");
  biovar.set_type(double_type(), 2);  // biovariability

  variable tlag("tlag_var_name");
  tlag.set_type(double_type(), 2);  // lag time

  variable rtol("rtol_var_name");
  rtol.set_type(double_type(), 0);  // rel tol

  variable atol("atol_var_name");
  atol.set_type(double_type(), 0);  // abs tol

  variable maxstep("maxstep_var_name");
  maxstep.set_type(int_type(), 0);  // max num step

  generalOdeModel_control so2(integration_function_name,
                              system_function_name,
                              nCmt, time, amt, rate, ii,
                              evid, cmt, addl, ss, pmatrix,
                              biovar, tlag, rtol, atol, maxstep);

  // dumb test to make sure we at least get the right types back
  EXPECT_EQ(integration_function_name , so2.integration_function_name_);
  EXPECT_EQ(system_function_name      , so2.system_function_name_);
  EXPECT_EQ(nCmt.type_                , so2.nCmt_.expression_type());
  EXPECT_EQ(time.type_                , so2.time_.expression_type());
  EXPECT_EQ(amt.type_                 , so2.amt_.expression_type());
  EXPECT_EQ(rate.type_                , so2.rate_.expression_type());
  EXPECT_EQ(ii.type_                  , so2.ii_.expression_type());
  EXPECT_EQ(evid.type_                , so2.evid_.expression_type());
  EXPECT_EQ(cmt.type_                 , so2.addl_.expression_type());
  EXPECT_EQ(ss.type_                  , so2.ss_.expression_type());
  EXPECT_EQ(pmatrix.type_             , so2.pMatrix_.expression_type());
  EXPECT_EQ(biovar.type_              , so2.biovar_.expression_type());
  EXPECT_EQ(tlag.type_                , so2.tlag_.expression_type());
  EXPECT_EQ(rtol.type_                , so2.rel_tol_.expression_type());
  EXPECT_EQ(atol.type_                , so2.abs_tol_.expression_type());
  EXPECT_EQ(maxstep.type_             , so2.max_num_steps_.expression_type());

  expression e2(so2);
  EXPECT_EQ(expr_type(matrix_type(), 0), e2.expression_type());
}

TEST(langAst, generalOdeModel) {
  using stan::lang::generalOdeModel;
  using stan::lang::variable;
  using stan::lang::expr_type;
  using stan::lang::expression;

  generalOdeModel so;

  std::string integration_function_name = "bar";
  std::string system_function_name = "foo";

  variable nCmt("nCmt_var_name");
  nCmt.set_type(int_type(), 0);  // int

  variable time("time_var_name");
  time.set_type(double_type(), 1);  // time of events

  variable amt("amt_var_name");
  amt.set_type(double_type(), 1);  // Amount of events

  variable rate("rate_var_name");
  rate.set_type(double_type(), 1);  // rate at events

  variable ii("ii_var_name");
  ii.set_type(double_type(), 1);  // interdose interval

  variable evid("evid_var_name");
  evid.set_type(int_type(), 1);  // event id

  variable cmt("cmt_var_name");
  cmt.set_type(int_type(), 1);  // compartment nb. of events

  variable addl("addl_var_name");
  addl.set_type(int_type(), 1);  // Number of additional doses

  variable ss("ss_var_name");
  ss.set_type(int_type(), 1);  // steady state flag

  variable pmatrix("pmatrix_var_name");
  pmatrix.set_type(double_type(), 2);  // linear ode system matrix

  variable biovar("biovar_var_name");
  biovar.set_type(double_type(), 2);  // biovariability

  variable tlag("tlag_var_name");
  tlag.set_type(double_type(), 2);  // lag time

  generalOdeModel so2(integration_function_name,
                              system_function_name,
                              nCmt, time, amt, rate, ii,
                              evid, cmt, addl, ss, pmatrix,
                              biovar, tlag);

  // dumb test to make sure we at least get the right types back
  EXPECT_EQ(integration_function_name , so2.integration_function_name_);
  EXPECT_EQ(system_function_name      , so2.system_function_name_);
  EXPECT_EQ(nCmt.type_                , so2.nCmt_.expression_type());
  EXPECT_EQ(time.type_                , so2.time_.expression_type());
  EXPECT_EQ(amt.type_                 , so2.amt_.expression_type());
  EXPECT_EQ(rate.type_                , so2.rate_.expression_type());
  EXPECT_EQ(ii.type_                  , so2.ii_.expression_type());
  EXPECT_EQ(evid.type_                , so2.evid_.expression_type());
  EXPECT_EQ(cmt.type_                 , so2.addl_.expression_type());
  EXPECT_EQ(ss.type_                  , so2.ss_.expression_type());
  EXPECT_EQ(pmatrix.type_             , so2.pMatrix_.expression_type());
  EXPECT_EQ(biovar.type_              , so2.biovar_.expression_type());
  EXPECT_EQ(tlag.type_                , so2.tlag_.expression_type());

  expression e2(so2);
  EXPECT_EQ(expr_type(matrix_type(), 0), e2.expression_type());
}
