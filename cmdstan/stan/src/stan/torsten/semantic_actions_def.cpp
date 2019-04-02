#ifndef STAN_LANG_TORSTEN_GRAMMARS_SEMANTIC_ACTIONS_DEF_CPP
#define STAN_LANG_TORSTEN_GRAMMARS_SEMANTIC_ACTIONS_DEF_CPP

#include <stan/io/program_reader.hpp>
#include <stan/lang/ast.hpp>
#include <stan/lang/grammars/iterator_typedefs.hpp>
#include <stan/lang/grammars/semantic_actions.hpp>
#include <stan/torsten/semantic_actions.hpp>
#include <boost/algorithm/string.hpp>
#include <boost/format.hpp>
#include <boost/spirit/include/qi.hpp>
#include <boost/variant/apply_visitor.hpp>
#include <boost/variant/recursive_variant.hpp>
#include <cstddef>
#include <limits>
#include <climits>
#include <iomanip>
#include <iostream>
#include <map>
#include <set>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

/**********************************
   univariate_integral
**********************************/

template <class T>
void validate_univariate_integral(const T& univar_fun,
                                    const variable_map& var_map,
                                  bool& pass,
                                  std::ostream& error_msgs) {
  pass = true;
  expr_type sys_result_type(double_type(), 0);
  std::vector<function_arg_type> sys_arg_types;

  if (univar_fun.integration_function_name_ == "univariate_integral_rk45"
      || univar_fun.integration_function_name_ == "univariate_integral_bdf"){
    sys_arg_types.push_back(function_arg_type(expr_type(double_type(), 0)));
    sys_arg_types.push_back(function_arg_type(expr_type(double_type(), 1)));
    sys_arg_types.push_back(function_arg_type(expr_type(double_type(), 1)));
    sys_arg_types.push_back(function_arg_type(expr_type(int_type(), 1)));
  }
  function_signature_t system_signature(sys_result_type, sys_arg_types);
  if (!function_signatures::instance()
      .is_defined(univar_fun.system_function_name_, system_signature)) {
    error_msgs << "first argument to "
               << univar_fun.integration_function_name_
               << " must be the name of a function with signature"
               << " (real, real[], real[], int[]) : real ";
    pass = false;
  }
  // test regular argument types
  if (univar_fun.t0_.expression_type() != expr_type(double_type(), 0)) {
    error_msgs << "second argument to "
               << univar_fun.integration_function_name_
               << " must have type real for time limit;"
               << " found type="
               << univar_fun.t0_.expression_type()
               << ". ";
    pass = false;
  }
  if (univar_fun.t1_.expression_type() != expr_type(double_type(), 0)) {
    error_msgs << "third argument to "
               << univar_fun.integration_function_name_
               << " must have type real for time limit;"
               << " found type="
               << univar_fun.t1_.expression_type()
               << ". ";
    pass = false;
  }
  if (univar_fun.theta_.expression_type() != expr_type(double_type(), 1)) {
    error_msgs << "fourth argument to "
               << univar_fun.integration_function_name_
               << " must have type real[] for parameters;"
               << " found type="
               << univar_fun.theta_.expression_type()
               << ". ";
    pass = false;
  }
  if (univar_fun.x_r_.expression_type() != expr_type(double_type(), 1)) {
    error_msgs << "fifth argument to "
               << univar_fun.integration_function_name_
               << " must have type real[] for real data; found type="
               << univar_fun.x_r_.expression_type()
               << ". ";
    pass = false;
  }
  if (univar_fun.x_i_.expression_type() != expr_type(int_type(), 1)) {
    error_msgs << "sixth argument to "
               << univar_fun.integration_function_name_
               << " must have type int[] for integer data; found type="
               << univar_fun.x_i_.expression_type()
               << ". ";
    pass = false;
  }
}

void validate_univariate_integral_control::operator()(
                                                      const univariate_integral_control& univar_fun,
                                                      const variable_map& var_map,
                                                      bool& pass,
                                                      std::ostream& error_msgs) const {
  validate_univariate_integral(univar_fun, var_map, pass, error_msgs);
}
boost::phoenix::function<validate_univariate_integral_control>
validate_univariate_integral_control_f;

bool data_only_expression::operator()(const univariate_integral_control& x)
  const {
  return boost::apply_visitor(*this, x.t0_.expr_)
    && boost::apply_visitor(*this, x.t1_.expr_)
    && boost::apply_visitor(*this, x.theta_.expr_);
}

template void assign_lhs::operator()(expression&,
                                     const univariate_integral_control&)
  const;

/**********************************
   generalOdeModel
**********************************/

template <class T>
void validate_generalOdeModel_non_control_args(const T& ode_fun,
                              const variable_map& var_map,
                              bool& pass,
                              std::ostream& error_msgs) {
  pass = true;

  expr_type sys_result_type(double_type(), 1);
  std::vector<function_arg_type> sys_arg_types;
  std::string expected_signature;

  // build expected function argument type for generalOdeModel
  if (ode_fun.integration_function_name_ == "generalOdeModel_rk45"
      || ode_fun.integration_function_name_ == "generalOdeModel_bdf") {
    sys_arg_types.push_back(function_arg_type(expr_type(double_type(),
                                                        0)));  // t0
    sys_arg_types.push_back(function_arg_type(expr_type(double_type(),
                                                        1)));  // y
    sys_arg_types.push_back(function_arg_type(expr_type(double_type(),
                                                        1)));  // theta
    sys_arg_types.push_back(function_arg_type(expr_type(double_type(),
                                                        1)));  // x_r
    sys_arg_types.push_back(function_arg_type(expr_type(int_type(),
                                                        1)));  // x_i
    expected_signature = "(real, real[], real[], real[], int[]) : real[]";
  }

  // build expected function argument type for mixOdeModel
  if ((ode_fun.integration_function_name_ == "mixOde1CptModel_rk45"
       || ode_fun.integration_function_name_ == "mixOde1CptModel_bdf")
      || (ode_fun.integration_function_name_ == "mixOde2CptModel_rk45"
          || ode_fun.integration_function_name_ == "mixOde2CptModel_bdf")) {
    sys_arg_types.push_back(function_arg_type(expr_type(double_type(),
                                                        0)));  // t0
    sys_arg_types.push_back(function_arg_type(expr_type(double_type(),
                                                        1)));  // y
    sys_arg_types.push_back(function_arg_type(expr_type(double_type(),
                                                        1)));  // y_PK
    sys_arg_types.push_back(function_arg_type(expr_type(double_type(),
                                                        1)));  // theta
    sys_arg_types.push_back(function_arg_type(expr_type(double_type(),
                                                        1)));  // x_r
    sys_arg_types.push_back(function_arg_type(expr_type(int_type(),
                                                        1)));  // x_i
    expected_signature = "(real, real[], real[], real[], real[], int[]) : real[]";  // NOLINT
  }

  function_signature_t system_signature(sys_result_type, sys_arg_types);

  // test function argument type
  if (!function_signatures::instance()
      .is_defined(ode_fun.system_function_name_, system_signature)) {
    error_msgs << "first argument to"
               << ode_fun.integration_function_name_
               << " must be a function with signature "
               << expected_signature << " ";
    pass = false;
  }

  // test regular argument types
  if (!ode_fun.nCmt_.expression_type().type().is_int_type()) {
    error_msgs << "second argument to "
               << ode_fun.integration_function_name_
               << " must be type int"
               << " for nCmt (number of compartments)"
               << "; found type="
               << ode_fun.nCmt_.expression_type()
               << ". ";
    pass = false;
  }
  if (ode_fun.time_.expression_type() != expr_type(double_type(), 1)) {
    error_msgs << "third argument to "
               << ode_fun.integration_function_name_
               << " must be type real[]"
               << "for time"
               << "; found type="
               << ode_fun.time_.expression_type()
               << ". ";
    pass = false;
  }
  if (ode_fun.amt_.expression_type() != expr_type(double_type(), 1)) {
    error_msgs << "fourth argument to "
               << ode_fun.integration_function_name_
               << " must be type real[]"
               << "for amount"
               << "; found type="
               << ode_fun.amt_.expression_type()
               << ". ";
    pass = false;
  }
  if (ode_fun.rate_.expression_type() != expr_type(double_type(), 1)) {
    error_msgs << "fifth argument to "
               << ode_fun.integration_function_name_
               << " must be type real[]"
               << "for rate"
               << "; found type="
               << ode_fun.rate_.expression_type()
               << ". ";
    pass = false;
  }
  if (ode_fun.ii_.expression_type() != expr_type(double_type(), 1)) {
    error_msgs << "sixth argument to "
               << ode_fun.integration_function_name_
               << " must be type real[]"
               << "for inter-dose interval"
               << "; found type="
               << ode_fun.ii_.expression_type()
               << ". ";
    pass = false;
  }
  if (ode_fun.evid_.expression_type() != expr_type(int_type(), 1)) {
    error_msgs << "seventh argument to "
               << ode_fun.integration_function_name_
               << " must be type int[]"
               << "for evid (event ID)"
               << "; found type="
               << ode_fun.evid_.expression_type()
               << ". ";
    pass = false;
  }
  if (ode_fun.cmt_.expression_type() != expr_type(int_type(), 1)) {
    error_msgs << "eighth argument to "
               << ode_fun.integration_function_name_
               << " must be type int[]"
               << "for cmt (compartment)"
               << "; found type="
               << ode_fun.cmt_.expression_type()
               << ". ";
    pass = false;
  }
  if (ode_fun.addl_.expression_type() != expr_type(int_type(), 1)) {
    error_msgs << "ninth argument to "
               << ode_fun.integration_function_name_
               << " must be type int[]"
               << "for addl (additional dose)"
               << "; found type="
               << ode_fun.addl_.expression_type()
               << ". ";
    pass = false;
  }
  if (ode_fun.ss_.expression_type() != expr_type(int_type(), 1)) {
    error_msgs << "tenth argument to "
               << ode_fun.integration_function_name_
               << " must be type int[]"
               << "for ss (steady state)"
               << "; found type="
               << ode_fun.ss_.expression_type()
               << ". ";
    pass = false;
  }
  if ((ode_fun.pMatrix_.expression_type() != expr_type(double_type(), 2))
      && (ode_fun.pMatrix_.expression_type() !=
          expr_type(double_type(), 1))) {
    error_msgs << "eleventh argument to "
               << ode_fun.integration_function_name_
               << " must be type real[ ] or real[ , ]"
               << " for the ODE parameters"
               << "; found type="
               << ode_fun.pMatrix_.expression_type()
               << ". ";
    pass = false;
  }
  if ((ode_fun.biovar_.expression_type() != expr_type(double_type(), 2))
      && (ode_fun.biovar_.expression_type() !=
          expr_type(double_type(), 1))) {
    error_msgs << "twelth argument to "
               << ode_fun.integration_function_name_
               << " must be type real[ ] or real[ , ]"
               << " for the bio-variability"
               << "; found type="
               << ode_fun.biovar_.expression_type()
               << ". ";
    pass = false;
  }
  if ((ode_fun.tlag_.expression_type() != expr_type(double_type(), 2))
      && (ode_fun.tlag_.expression_type() != expr_type(double_type(), 1))) {
    error_msgs << "thirteenth argument to "
               << ode_fun.integration_function_name_
               << " must be type real[ ] or real[ , ]"
               << " for the lag times"
               << "; found type="
               << ode_fun.tlag_.expression_type()
               << ". ";
    pass = false;
  }

  // test data-only variables do not have parameters (int locals OK)
  if (has_var(ode_fun.nCmt_, var_map)) {
    error_msgs << "second argument to "
               << ode_fun.integration_function_name_
               << " for number of compartments"
               << " must be data only and not reference parameters";
    pass = false;
  }
  if (has_var(ode_fun.evid_, var_map)) {
    error_msgs << "seventh argument to "
               << ode_fun.integration_function_name_
               << " for event ID (evid)"
               << " must be data only and not reference parameters";
    pass = false;
  }
  if (has_var(ode_fun.cmt_, var_map)) {
    error_msgs << "eighth argument to "
               << ode_fun.integration_function_name_
               << " for compartment number (cmt)"
               << " must be data only and not reference parameters";
    pass = false;
  }
  if (has_var(ode_fun.addl_, var_map)) {
    error_msgs << "ninth argument to "
               << ode_fun.integration_function_name_
               << " for additional dose (addl)"
               << " must be data only and not reference parameters";
    pass = false;
  }
  if (has_var(ode_fun.ss_, var_map)) {
    error_msgs << "tenth argument to "
               << ode_fun.integration_function_name_
               << " for steady state approximation (ss)"
               << " must be data only and not reference parameters";
    pass = false;
  }
}

template <class T>
void validate_generalOdeModel_control_args(const T& ode_fun,
                                           const variable_map& var_map,
                                           bool& pass,
                                           std::ostream& error_msgs) {
  if (!ode_fun.rel_tol_.expression_type().is_primitive()) {
    error_msgs << "fourteenth argument to "
               << ode_fun.integration_function_name_
               << " must be type real or int"
               << " for relative tolerance"
               << "; found type="
               << ode_fun.rel_tol_.expression_type()
               << ". ";
    pass = false;
  }
  if (!ode_fun.abs_tol_.expression_type().is_primitive()) {
    error_msgs << "fifthteenth argument to "
               << ode_fun.integration_function_name_
               << " must be type real or int"
               << " for absolute tolerance"
               << "; found type="
               << ode_fun.abs_tol_.expression_type()
               << ". ";
    pass = false;
  }
  if (!ode_fun.max_num_steps_.expression_type().is_primitive()) {
    error_msgs << "sixteenth argument to "
               << ode_fun.integration_function_name_
               << " must be type real or int"
               << " for maximum number of steps"
               << "; found type="
               << ode_fun.max_num_steps_.expression_type()
               << ". ";
    pass = false;
  }

  // test data-only variables do not have parameters (int locals OK)
  if (has_var(ode_fun.rel_tol_, var_map)) {
    error_msgs << "fourteenth argument to "
               << ode_fun.integration_function_name_
               << " for relative tolerance"
               << " must be data only and not reference parameters";
    pass = false;
  }
  if (has_var(ode_fun.abs_tol_, var_map)) {
    error_msgs << "fifthteenth argument to "
               << ode_fun.integration_function_name_
               << " for absolute tolerance"
               << " must be data only and not reference parameters";
    pass = false;
  }
  if (has_var(ode_fun.max_num_steps_, var_map)) {
    error_msgs << "sixteenth argument to "
               << ode_fun.integration_function_name_
               << " for maximum number of steps"
               << " must be data only and not reference parameters";
    pass = false;
  }
}

/**********************************
   pmx_solve_group
**********************************/
template <class T>
void validate_pmx_solve_group_non_control_args(const T& ode_fun,
                              const variable_map& var_map,
                              bool& pass,
                              std::ostream& error_msgs) {
  pass = true;

  expr_type sys_result_type(double_type(), 1);
  std::vector<function_arg_type> sys_arg_types;
  std::string expected_signature;

  // build expected function argument type for generalOdeModel
  if (ode_fun.integration_function_name_ == "pmx_solve_group_rk45"
      || ode_fun.integration_function_name_ == "pmx_solve_group_adams"
      || ode_fun.integration_function_name_ == "pmx_solve_group_bdf") {
    sys_arg_types.push_back(function_arg_type(expr_type(double_type(),
                                                        0)));  // t0
    sys_arg_types.push_back(function_arg_type(expr_type(double_type(),
                                                        1)));  // y
    sys_arg_types.push_back(function_arg_type(expr_type(double_type(),
                                                        1)));  // theta
    sys_arg_types.push_back(function_arg_type(expr_type(double_type(),
                                                        1)));  // x_r
    sys_arg_types.push_back(function_arg_type(expr_type(int_type(),
                                                        1)));  // x_i
    expected_signature = "(real, real[], real[], real[], int[]) : real[]";
  }

  // build expected function argument type for mixOdeModel
  if ((ode_fun.integration_function_name_ == "mixOde1CptModel_rk45"
       || ode_fun.integration_function_name_ == "mixOde1CptModel_bdf")
      || (ode_fun.integration_function_name_ == "mixOde2CptModel_rk45"
          || ode_fun.integration_function_name_ == "mixOde2CptModel_bdf")) {
    sys_arg_types.push_back(function_arg_type(expr_type(double_type(),
                                                        0)));  // t0
    sys_arg_types.push_back(function_arg_type(expr_type(double_type(),
                                                        1)));  // y
    sys_arg_types.push_back(function_arg_type(expr_type(double_type(),
                                                        1)));  // y_PK
    sys_arg_types.push_back(function_arg_type(expr_type(double_type(),
                                                        1)));  // theta
    sys_arg_types.push_back(function_arg_type(expr_type(double_type(),
                                                        1)));  // x_r
    sys_arg_types.push_back(function_arg_type(expr_type(int_type(),
                                                        1)));  // x_i
    expected_signature = "(real, real[], real[], real[], real[], int[]) : real[]";  // NOLINT
  }

  function_signature_t system_signature(sys_result_type, sys_arg_types);

  // test function argument type
  if (!function_signatures::instance()
      .is_defined(ode_fun.system_function_name_, system_signature)) {
    error_msgs << "1st argument to"
               << ode_fun.integration_function_name_
               << " must be a function with signature "
               << expected_signature << " ";
    pass = false;
  }

  // test regular argument types
  if (!ode_fun.nCmt_.expression_type().type().is_int_type()) {
    error_msgs << "2nd argument to "
               << ode_fun.integration_function_name_
               << " must be type int"
               << " for nCmt (number of compartments)"
               << "; found type="
               << ode_fun.nCmt_.expression_type()
               << ". ";
    pass = false;
  }
  if (ode_fun.len_.expression_type() != expr_type(int_type(), 1)) {
    error_msgs << "3rd argument to "
               << ode_fun.integration_function_name_
               << " must be type int[]"
               << " for len"
               << "; found type="
               << ode_fun.time_.expression_type()
               << ". ";
    pass = false;
  }
  if (ode_fun.time_.expression_type() != expr_type(double_type(), 1)) {
    error_msgs << "4th argument to "
               << ode_fun.integration_function_name_
               << " must be type real[]"
               << "for time"
               << "; found type="
               << ode_fun.time_.expression_type()
               << ". ";
    pass = false;
  }
  if (ode_fun.amt_.expression_type() != expr_type(double_type(), 1)) {
    error_msgs << "5th argument to "
               << ode_fun.integration_function_name_
               << " must be type real[]"
               << "for amount"
               << "; found type="
               << ode_fun.amt_.expression_type()
               << ". ";
    pass = false;
  }
  if (ode_fun.rate_.expression_type() != expr_type(double_type(), 1)) {
    error_msgs << "6th argument to "
               << ode_fun.integration_function_name_
               << " must be type real[]"
               << "for rate"
               << "; found type="
               << ode_fun.rate_.expression_type()
               << ". ";
    pass = false;
  }
  if (ode_fun.ii_.expression_type() != expr_type(double_type(), 1)) {
    error_msgs << "7th argument to "
               << ode_fun.integration_function_name_
               << " must be type real[]"
               << "for inter-dose interval"
               << "; found type="
               << ode_fun.ii_.expression_type()
               << ". ";
    pass = false;
  }
  if (ode_fun.evid_.expression_type() != expr_type(int_type(), 1)) {
    error_msgs << "8th argument to "
               << ode_fun.integration_function_name_
               << " must be type int[]"
               << "for evid (event ID)"
               << "; found type="
               << ode_fun.evid_.expression_type()
               << ". ";
    pass = false;
  }
  if (ode_fun.cmt_.expression_type() != expr_type(int_type(), 1)) {
    error_msgs << "9th argument to "
               << ode_fun.integration_function_name_
               << " must be type int[]"
               << "for cmt (compartment)"
               << "; found type="
               << ode_fun.cmt_.expression_type()
               << ". ";
    pass = false;
  }
  if (ode_fun.addl_.expression_type() != expr_type(int_type(), 1)) {
    error_msgs << "10th argument to "
               << ode_fun.integration_function_name_
               << " must be type int[]"
               << "for addl (additional dose)"
               << "; found type="
               << ode_fun.addl_.expression_type()
               << ". ";
    pass = false;
  }
  if (ode_fun.ss_.expression_type() != expr_type(int_type(), 1)) {
    error_msgs << "11th argument to "
               << ode_fun.integration_function_name_
               << " must be type int[]"
               << "for ss (steady state)"
               << "; found type="
               << ode_fun.ss_.expression_type()
               << ". ";
    pass = false;
  }
  if ((ode_fun.pMatrix_.expression_type() != expr_type(double_type(), 2))
      && (ode_fun.pMatrix_.expression_type() !=
          expr_type(double_type(), 2))) {
    error_msgs << "12th argument to "
               << ode_fun.integration_function_name_
               << " must be type real[ ] or real[ , ]"
               << " for the ODE parameters"
               << "; found type="
               << ode_fun.pMatrix_.expression_type()
               << ". ";
    pass = false;
  }
  if ((ode_fun.biovar_.expression_type() != expr_type(double_type(), 2))
      && (ode_fun.biovar_.expression_type() !=
          expr_type(double_type(), 2))) {
    error_msgs << "13th argument to "
               << ode_fun.integration_function_name_
               << " must be type real[ ] or real[ , ]"
               << " for the bio-variability"
               << "; found type="
               << ode_fun.biovar_.expression_type()
               << ". ";
    pass = false;
  }
  if ((ode_fun.tlag_.expression_type() != expr_type(double_type(), 2))
      && (ode_fun.tlag_.expression_type() != expr_type(double_type(), 2))) {
    error_msgs << "14th argument to "
               << ode_fun.integration_function_name_
               << " must be type real[ ] or real[ , ]"
               << " for the lag times"
               << "; found type="
               << ode_fun.tlag_.expression_type()
               << ". ";
    pass = false;
  }

  // test data-only variables do not have parameters (int locals OK)
  if (has_var(ode_fun.nCmt_, var_map)) {
    error_msgs << "second argument to "
               << ode_fun.integration_function_name_
               << " for number of compartments"
               << " must be data only and not reference parameters";
    pass = false;
  }
  if (has_var(ode_fun.evid_, var_map)) {
    error_msgs << "seventh argument to "
               << ode_fun.integration_function_name_
               << " for event ID (evid)"
               << " must be data only and not reference parameters";
    pass = false;
  }
  if (has_var(ode_fun.cmt_, var_map)) {
    error_msgs << "eighth argument to "
               << ode_fun.integration_function_name_
               << " for compartment number (cmt)"
               << " must be data only and not reference parameters";
    pass = false;
  }
  if (has_var(ode_fun.addl_, var_map)) {
    error_msgs << "ninth argument to "
               << ode_fun.integration_function_name_
               << " for additional dose (addl)"
               << " must be data only and not reference parameters";
    pass = false;
  }
  if (has_var(ode_fun.ss_, var_map)) {
    error_msgs << "tenth argument to "
               << ode_fun.integration_function_name_
               << " for steady state approximation (ss)"
               << " must be data only and not reference parameters";
    pass = false;
  }
}


void validate_generalOdeModel_control::operator()(
                      const generalOdeModel_control& ode_fun,
                      const variable_map& var_map,
                      bool& pass,
                      std::ostream& error_msgs) const {
  validate_generalOdeModel_non_control_args(ode_fun, var_map, pass, error_msgs);
  validate_generalOdeModel_control_args(ode_fun, var_map, pass, error_msgs);
}
boost::phoenix::function<validate_generalOdeModel_control>
validate_generalOdeModel_control_f;

bool data_only_expression::operator()(const generalOdeModel_control& x)
  const {
  return ((((((boost::apply_visitor(*this, x.time_.expr_)
               && boost::apply_visitor(*this, x.amt_.expr_)))
             && boost::apply_visitor(*this, x.rate_.expr_)
             && boost::apply_visitor(*this, x.ii_.expr_))
            && boost::apply_visitor(*this, x.pMatrix_.expr_))
           && boost::apply_visitor(*this, x.biovar_.expr_))
          && boost::apply_visitor(*this, x.tlag_.expr_));
}  // include all arguments with a template type

template void assign_lhs::operator()(expression&,
                                     const generalOdeModel_control&)
  const;

void validate_generalOdeModel::operator()(
                      const generalOdeModel& ode_fun,
                      const variable_map& var_map,
                      bool& pass,
                      std::ostream& error_msgs) const {
  validate_generalOdeModel_non_control_args(ode_fun, var_map, pass, error_msgs);
}
boost::phoenix::function<validate_generalOdeModel>
validate_generalOdeModel_f;

bool data_only_expression::operator()(const generalOdeModel& x)
  const {
  return ((((((boost::apply_visitor(*this, x.time_.expr_)
               && boost::apply_visitor(*this, x.amt_.expr_)))
             && boost::apply_visitor(*this, x.rate_.expr_)
             && boost::apply_visitor(*this, x.ii_.expr_))
            && boost::apply_visitor(*this, x.pMatrix_.expr_))
           && boost::apply_visitor(*this, x.biovar_.expr_))
          && boost::apply_visitor(*this, x.tlag_.expr_));
}

template void assign_lhs::operator()(expression&,
                                     const generalOdeModel&)
  const;

// pop pk
void validate_pmx_solve_group::operator()(
                      const pmx_solve_group& ode_fun,
                      const variable_map& var_map,
                      bool& pass,
                      std::ostream& error_msgs) const {
  validate_pmx_solve_group_non_control_args(ode_fun, var_map, pass, error_msgs);
}
boost::phoenix::function<validate_pmx_solve_group>
validate_pmx_solve_group_f;

bool data_only_expression::operator()(const pmx_solve_group& x)
  const {
  return ((((((boost::apply_visitor(*this, x.time_.expr_)
               && boost::apply_visitor(*this, x.amt_.expr_)))
             && boost::apply_visitor(*this, x.rate_.expr_)
             && boost::apply_visitor(*this, x.ii_.expr_))
            && boost::apply_visitor(*this, x.pMatrix_.expr_))
           && boost::apply_visitor(*this, x.biovar_.expr_))
          && boost::apply_visitor(*this, x.tlag_.expr_));
}

template void assign_lhs::operator()(expression&,
                                     const pmx_solve_group&)
  const;

#endif
