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

struct torsten_types {
  static const double_type t_dbl;
  static const int_type t_int;
  static const vector_type t_vec;
  static const bare_expr_type t_dbl_1;
  static const bare_expr_type t_int_1;
  static const bare_expr_type t_dbl_2;
  static const bare_expr_type t_int_2;
  const bare_expr_type sys_result_type;

  torsten_types(bare_expr_type return_t) :
    sys_result_type(return_t)
  {}
};

const double_type torsten_types::t_dbl = double_type();
const int_type torsten_types::t_int = int_type();
const vector_type torsten_types::t_vec = vector_type();
const bare_expr_type torsten_types::t_dbl_1 = bare_array_type(t_dbl, 1);
const bare_expr_type torsten_types::t_int_1 = bare_array_type(t_int, 1);
const bare_expr_type torsten_types::t_dbl_2 = bare_array_type(t_dbl, 2);
const bare_expr_type torsten_types::t_int_2 = bare_array_type(t_int, 2);

/**********************************
   univariate_integral
**********************************/

template <class T>
void validate_univariate_integral(const T& univar_fun,
                                  const variable_map& var_map,
                                  bool& pass,
                                  std::ostream& error_msgs) {
pass = true;
torsten_types pmx_t(torsten_types::t_dbl);
  std::vector<bare_expr_type> sys_arg_types;

  if (univar_fun.integration_function_name_ == "univariate_integral_rk45"
      || univar_fun.integration_function_name_ == "univariate_integral_bdf"){
    sys_arg_types.push_back(torsten_types::t_dbl);
    sys_arg_types.push_back(torsten_types::t_dbl_1);
    sys_arg_types.push_back(torsten_types::t_dbl_1);
    sys_arg_types.push_back(torsten_types::t_int_1);
  }
  function_signature_t system_signature(pmx_t.sys_result_type, sys_arg_types);
  if (!function_signatures::instance()
      .is_defined(univar_fun.system_function_name_, system_signature)) {
    error_msgs << "first argument to "
               << univar_fun.integration_function_name_
               << " must be the name of a function with signature"
               << " (real, real[], real[], int[]) : real ";
    pass = false;
  }
  // test regular argument types
  if (univar_fun.t0_.bare_type() != torsten_types::t_dbl) {
    error_msgs << "second argument to "
               << univar_fun.integration_function_name_
               << " must have type real for time limit;"
               << " found type="
               << univar_fun.t0_.bare_type()
               << ". ";
    pass = false;
  }
  if (univar_fun.t1_.bare_type() != torsten_types::t_dbl) {
    error_msgs << "third argument to "
               << univar_fun.integration_function_name_
               << " must have type real for time limit;"
               << " found type="
               << univar_fun.t1_.bare_type()
               << ". ";
    pass = false;
  }
  if (univar_fun.theta_.bare_type() != torsten_types::t_dbl_1) {
    error_msgs << "fourth argument to "
               << univar_fun.integration_function_name_
               << " must have type real[] for parameters;"
               << " found type="
               << univar_fun.theta_.bare_type()
               << ". ";
    pass = false;
  }
  if (univar_fun.x_r_.bare_type() != torsten_types::t_dbl_1) {
    error_msgs << "fifth argument to "
               << univar_fun.integration_function_name_
               << " must have type real[] for real data; found type="
               << univar_fun.x_r_.bare_type()
               << ". ";
    pass = false;
  }
  if (univar_fun.x_i_.bare_type() != torsten_types::t_int_1) {
    error_msgs << "sixth argument to "
               << univar_fun.integration_function_name_
               << " must have type int[] for integer data; found type="
               << univar_fun.x_i_.bare_type()
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

  torsten_types pmx_t(torsten_types::t_dbl_1);
  std::vector<bare_expr_type> sys_arg_types;
  std::string expected_signature;

  // build expected function argument type for generalOdeModel
  if (ode_fun.integration_function_name_ == "generalOdeModel_rk45"
      || ode_fun.integration_function_name_ == "generalOdeModel_bdf"
      || ode_fun.integration_function_name_ == "pmx_solve_adams"
      || ode_fun.integration_function_name_ == "pmx_solve_bdf"
      || ode_fun.integration_function_name_ == "pmx_solve_rk45"
      ) {
    sys_arg_types.push_back(torsten_types::t_dbl);  // t0
    sys_arg_types.push_back(torsten_types::t_dbl_1);  // y
    sys_arg_types.push_back(torsten_types::t_dbl_1);  // theta
    sys_arg_types.push_back(torsten_types::t_dbl_1);  // x_r
    sys_arg_types.push_back(torsten_types::t_int_1);  // x_i
    expected_signature = "(real, real[], real[], real[], int[]) : real[]";
  }

  // build expected function argument type for mixOdeModel
  if (ode_fun.integration_function_name_ == "mixOde1CptModel_rk45"
      || ode_fun.integration_function_name_ == "mixOde1CptModel_bdf"
      || ode_fun.integration_function_name_ == "mixOde2CptModel_rk45"
      || ode_fun.integration_function_name_ == "mixOde2CptModel_bdf"
      || ode_fun.integration_function_name_ == "pmx_solve_onecpt_bdf"
      || ode_fun.integration_function_name_ == "pmx_solve_onecpt_rk45"
      || ode_fun.integration_function_name_ == "pmx_solve_twocpt_bdf"
      || ode_fun.integration_function_name_ == "pmx_solve_twocpt_rk45"
      ) {
    sys_arg_types.push_back(torsten_types::t_dbl);  // t0
    sys_arg_types.push_back(torsten_types::t_dbl_1);  // y
    sys_arg_types.push_back(torsten_types::t_dbl_1);  // y_PK
    sys_arg_types.push_back(torsten_types::t_dbl_1);  // theta
    sys_arg_types.push_back(torsten_types::t_dbl_1);  // x_r
    sys_arg_types.push_back(torsten_types::t_int_1);  // x_i
    expected_signature = "(real, real[], real[], real[], real[], int[]) : real[]";  // NOLINT
  }

  function_signature_t system_signature(pmx_t.sys_result_type, sys_arg_types);

  // test function argument type
  if (!function_signatures::instance()
      .is_defined(ode_fun.system_function_name_, system_signature)) {
    error_msgs << "1st argument to "
               << ode_fun.integration_function_name_
               << " must be a function with signature "
               << expected_signature << " ";
    pass = false;
  }

  // test regular argument types
  if (!ode_fun.nCmt_.bare_type().is_int_type()) {
    error_msgs << "2nd argument to "
               << ode_fun.integration_function_name_
               << " must be type int"
               << " for number of compartments"
               << "; found type="
               << ode_fun.nCmt_.bare_type()
               << ". ";
    pass = false;
  }
  if (ode_fun.time_.bare_type() != torsten_types::t_dbl_1) {
    error_msgs << "3rd argument to "
               << ode_fun.integration_function_name_
               << " must be type real[]"
               << " for time"
               << "; found type="
               << ode_fun.time_.bare_type()
               << ". ";
    pass = false;
  }
  if (ode_fun.amt_.bare_type() != torsten_types::t_dbl_1) {
    error_msgs << "4th argument to "
               << ode_fun.integration_function_name_
               << " must be type real[]"
               << " for amount"
               << "; found type="
               << ode_fun.amt_.bare_type()
               << ". ";
    pass = false;
  }
  if (ode_fun.rate_.bare_type() != torsten_types::t_dbl_1) {
    error_msgs << "5th argument to "
               << ode_fun.integration_function_name_
               << " must be type real[]"
               << " for rate"
               << "; found type="
               << ode_fun.rate_.bare_type()
               << ". ";
    pass = false;
  }
  if (ode_fun.ii_.bare_type() != torsten_types::t_dbl_1) {
    error_msgs << "6th argument to "
               << ode_fun.integration_function_name_
               << " must be type real[]"
               << " for inter-dose interval"
               << "; found type="
               << ode_fun.ii_.bare_type()
               << ". ";
    pass = false;
  }
  if (ode_fun.evid_.bare_type() != torsten_types::t_int_1) {
    error_msgs << "7th argument to "
               << ode_fun.integration_function_name_
               << " must be type int[]"
               << " for event ID"
               << "; found type="
               << ode_fun.evid_.bare_type()
               << ". ";
    pass = false;
  }
  if (ode_fun.cmt_.bare_type() != torsten_types::t_int_1) {
    error_msgs << "8th argument to "
               << ode_fun.integration_function_name_
               << " must be type int[]"
               << " for compartment ID"
               << "; found type="
               << ode_fun.cmt_.bare_type()
               << ". ";
    pass = false;
  }
  if (ode_fun.addl_.bare_type() != torsten_types::t_int_1) {
    error_msgs << "9th argument to "
               << ode_fun.integration_function_name_
               << " must be type int[]"
               << " for number of additional doses"
               << "; found type="
               << ode_fun.addl_.bare_type()
               << ". ";
    pass = false;
  }
  if (ode_fun.ss_.bare_type() != torsten_types::t_int_1) {
    error_msgs << "10th argument to "
               << ode_fun.integration_function_name_
               << " must be type int[]"
               << " for steady state flags"
               << "; found type="
               << ode_fun.ss_.bare_type()
               << ". ";
    pass = false;
  }
  if ((ode_fun.pMatrix_.bare_type() != torsten_types::t_dbl_2)
      && (ode_fun.pMatrix_.bare_type() != torsten_types::t_dbl_1)) {
    error_msgs << "11th argument to "
               << ode_fun.integration_function_name_
               << " must be type real[] or real[ , ]"
               << " for ODE parameters"
               << "; found type="
               << ode_fun.pMatrix_.bare_type()
               << ". ";
    pass = false;
  }
  if ((ode_fun.biovar_.bare_type() != torsten_types::t_dbl_2)
      && (ode_fun.biovar_.bare_type() != torsten_types::t_dbl_1)) {
    error_msgs << "12th argument to "
               << ode_fun.integration_function_name_
               << " must be type real[] or real[ , ]"
               << " for bioavailability"
               << "; found type="
               << ode_fun.biovar_.bare_type()
               << ". ";
    pass = false;
  }
  if ((ode_fun.tlag_.bare_type() != torsten_types::t_dbl_2)
      && (ode_fun.tlag_.bare_type() != torsten_types::t_dbl_1)) {
    error_msgs << "13th argument to "
               << ode_fun.integration_function_name_
               << " must be type real[] or real[ , ]"
               << " for lag times"
               << "; found type="
               << ode_fun.tlag_.bare_type()
               << ". ";
    pass = false;
  }

  // test data-only variables do not have parameters (int locals OK)
  if (has_var(ode_fun.nCmt_, var_map)) {
    error_msgs << "2nd argument to "
               << ode_fun.integration_function_name_
               << " for number of compartments"
               << " must be data only and not reference parameters";
    pass = false;
  }
  if (has_var(ode_fun.evid_, var_map)) {
    error_msgs << "3rd argument to "
               << ode_fun.integration_function_name_
               << " for event ID"
               << " must be data only and not reference parameters";
    pass = false;
  }
  if (has_var(ode_fun.cmt_, var_map)) {
    error_msgs << "8th argument to "
               << ode_fun.integration_function_name_
               << " for compartment ID"
               << " must be data only and not reference parameters";
    pass = false;
  }
  if (has_var(ode_fun.addl_, var_map)) {
    error_msgs << "9th argument to "
               << ode_fun.integration_function_name_
               << " for number of additional doses"
               << " must be data only and not reference parameters";
    pass = false;
  }
  if (has_var(ode_fun.ss_, var_map)) {
    error_msgs << "10th argument to "
               << ode_fun.integration_function_name_
               << " for steady state flags"
               << " must be data only and not reference parameters";
    pass = false;
  }
}

template <class T>
void validate_generalOdeModel_control_args(const T& ode_fun,
                                           const variable_map& var_map,
                                           bool& pass,
                                           std::ostream& error_msgs) {
  if (!ode_fun.rel_tol_.bare_type().is_primitive()) {
    error_msgs << "14th argument to "
               << ode_fun.integration_function_name_
               << " must be type real or int"
               << " for relative tolerance"
               << "; found type="
               << ode_fun.rel_tol_.bare_type()
               << ". ";
    pass = false;
  }
  if (!ode_fun.abs_tol_.bare_type().is_primitive()) {
    error_msgs << "15th argument to "
               << ode_fun.integration_function_name_
               << " must be type real or int"
               << " for absolute tolerance"
               << "; found type="
               << ode_fun.abs_tol_.bare_type()
               << ". ";
    pass = false;
  }
  if (!ode_fun.max_num_steps_.bare_type().is_primitive()) {
    error_msgs << "16th argument to "
               << ode_fun.integration_function_name_
               << " must be type real or int"
               << " for maximum number of steps"
               << "; found type="
               << ode_fun.max_num_steps_.bare_type()
               << ". ";
    pass = false;
  }

  // test data-only variables do not have parameters (int locals OK)
  if (has_var(ode_fun.rel_tol_, var_map)) {
    error_msgs << "14th argument to "
               << ode_fun.integration_function_name_
               << " for relative tolerance"
               << " must be data only and not reference parameters";
    pass = false;
  }
  if (has_var(ode_fun.abs_tol_, var_map)) {
    error_msgs << "15th argument to "
               << ode_fun.integration_function_name_
               << " for absolute tolerance"
               << " must be data only and not reference parameters";
    pass = false;
  }
  if (has_var(ode_fun.max_num_steps_, var_map)) {
    error_msgs << "16th argument to "
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

  torsten_types pmx_t(torsten_types::t_dbl_1);
  std::vector<bare_expr_type> sys_arg_types;
  std::string expected_signature;

  // build expected function argument type for generalOdeModel
  if (   ode_fun.integration_function_name_ == "pmx_solve_group_rk45"
      || ode_fun.integration_function_name_ == "pmx_solve_group_adams"
      || ode_fun.integration_function_name_ == "pmx_solve_group_bdf") {
    sys_arg_types.push_back(torsten_types::t_dbl);  // t0
    sys_arg_types.push_back(torsten_types::t_dbl_1);  // y
    sys_arg_types.push_back(torsten_types::t_dbl_1);  // theta
    sys_arg_types.push_back(torsten_types::t_dbl_1);  // x_r
    sys_arg_types.push_back(torsten_types::t_int_1);  // x_i
    expected_signature = "(real, real[], real[], real[], int[]) : real[]";
  }

  // build expected function argument type for mixOdeModel
  if (ode_fun.integration_function_name_ == "pmx_solve_group_onecpt_rk45"
      || ode_fun.integration_function_name_ == "pmx_solve_group_onecpt_bdf"
      || ode_fun.integration_function_name_ == "pmx_solve_group_twocpt_rk45"
      || ode_fun.integration_function_name_ == "pmx_solve_group_twocpt_bdf") {
    sys_arg_types.push_back(torsten_types::t_dbl);  // t0
    sys_arg_types.push_back(torsten_types::t_dbl_1);  // y
    sys_arg_types.push_back(torsten_types::t_dbl_1);  // y_PK
    sys_arg_types.push_back(torsten_types::t_dbl_1);  // theta
    sys_arg_types.push_back(torsten_types::t_dbl_1);  // x_r
    sys_arg_types.push_back(torsten_types::t_int_1);  // x_i
    expected_signature = "(real, real[], real[], real[], real[], int[]) : real[]";  // NOLINT
  }

  function_signature_t system_signature(pmx_t.sys_result_type, sys_arg_types);

  // test function argument type
  if (!function_signatures::instance()
      .is_defined(ode_fun.system_function_name_, system_signature)) {
    error_msgs << "1st argument to "
               << ode_fun.integration_function_name_
               << " must be a function with signature "
               << expected_signature << " ";
    pass = false;
  }

  // test regular argument types
  if (!ode_fun.nCmt_.bare_type().is_int_type()) {
    error_msgs << "2nd argument to "
               << ode_fun.integration_function_name_
               << " must be type int"
               << " for nCmt (number of compartments)"
               << "; found type="
               << ode_fun.nCmt_.bare_type()
               << ". ";
    pass = false;
  }
  if (ode_fun.len_.bare_type() != torsten_types::t_int_1) {
    error_msgs << "3rd argument to "
               << ode_fun.integration_function_name_
               << " must be type int[]"
               << " for len"
               << "; found type="
               << ode_fun.time_.bare_type()
               << ". ";
    pass = false;
  }
  if (ode_fun.time_.bare_type() != torsten_types::t_dbl_1) {
    error_msgs << "4th argument to "
               << ode_fun.integration_function_name_
               << " must be type real[]"
               << " for time"
               << "; found type="
               << ode_fun.time_.bare_type()
               << ". ";
    pass = false;
  }
  if (ode_fun.amt_.bare_type() != torsten_types::t_dbl_1) {
    error_msgs << "5th argument to "
               << ode_fun.integration_function_name_
               << " must be type real[]"
               << " for amount"
               << "; found type="
               << ode_fun.amt_.bare_type()
               << ". ";
    pass = false;
  }
  if (ode_fun.rate_.bare_type() != torsten_types::t_dbl_1) {
    error_msgs << "6th argument to "
               << ode_fun.integration_function_name_
               << " must be type real[]"
               << " for rate"
               << "; found type="
               << ode_fun.rate_.bare_type()
               << ". ";
    pass = false;
  }
  if (ode_fun.ii_.bare_type() != torsten_types::t_dbl_1) {
    error_msgs << "7th argument to "
               << ode_fun.integration_function_name_
               << " must be type real[]"
               << " for inter-dose interval"
               << "; found type="
               << ode_fun.ii_.bare_type()
               << ". ";
    pass = false;
  }
  if (ode_fun.evid_.bare_type() != torsten_types::t_int_1) {
    error_msgs << "8th argument to "
               << ode_fun.integration_function_name_
               << " must be type int[]"
               << " for event ID"
               << "; found type="
               << ode_fun.evid_.bare_type()
               << ". ";
    pass = false;
  }
  if (ode_fun.cmt_.bare_type() != torsten_types::t_int_1) {
    error_msgs << "9th argument to "
               << ode_fun.integration_function_name_
               << " must be type int[]"
               << " for compartment ID"
               << "; found type="
               << ode_fun.cmt_.bare_type()
               << ". ";
    pass = false;
  }
  if (ode_fun.addl_.bare_type() != torsten_types::t_int_1) {
    error_msgs << "10th argument to "
               << ode_fun.integration_function_name_
               << " must be type int[]"
               << " for number of additional doses"
               << "; found type="
               << ode_fun.addl_.bare_type()
               << ". ";
    pass = false;
  }
  if (ode_fun.ss_.bare_type() != torsten_types::t_int_1) {
    error_msgs << "11th argument to "
               << ode_fun.integration_function_name_
               << " must be type int[]"
               << " for steady state flags"
               << "; found type="
               << ode_fun.ss_.bare_type()
               << ". ";
    pass = false;
  }
  if ((ode_fun.pMatrix_.bare_type() != torsten_types::t_dbl_2)
      && (ode_fun.pMatrix_.bare_type() != torsten_types::t_dbl_1)) {
    error_msgs << "12th argument to "
               << ode_fun.integration_function_name_
               << " must be type real[] or real[ , ]"
               << " for ODE parameters"
               << "; found type="
               << ode_fun.pMatrix_.bare_type()
               << ". ";
    pass = false;
  }
  if ((ode_fun.biovar_.bare_type() != torsten_types::t_dbl_2)
      && (ode_fun.biovar_.bare_type() != torsten_types::t_dbl_1)) {
    error_msgs << "13th argument to "
               << ode_fun.integration_function_name_
               << " must be type real[] or real[ , ]"
               << " for bioavailability"
               << "; found type="
               << ode_fun.biovar_.bare_type()
               << ". ";
    pass = false;
  }
  if ((ode_fun.tlag_.bare_type() != torsten_types::t_dbl_2)
      && (ode_fun.tlag_.bare_type() != torsten_types::t_dbl_1)) {
    error_msgs << "14th argument to "
               << ode_fun.integration_function_name_
               << " must be type real[] or real[ , ]"
               << " for lag times"
               << "; found type="
               << ode_fun.tlag_.bare_type()
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
               << " for event ID"
               << " must be data only and not reference parameters";
    pass = false;
  }
  if (has_var(ode_fun.cmt_, var_map)) {
    error_msgs << "eighth argument to "
               << ode_fun.integration_function_name_
               << " for compartment ID"
               << " must be data only and not reference parameters";
    pass = false;
  }
  if (has_var(ode_fun.addl_, var_map)) {
    error_msgs << "ninth argument to "
               << ode_fun.integration_function_name_
               << " for number of additional doses"
               << " must be data only and not reference parameters";
    pass = false;
  }
  if (has_var(ode_fun.ss_, var_map)) {
    error_msgs << "tenth argument to "
               << ode_fun.integration_function_name_
               << " for steady state flags"
               << " must be data only and not reference parameters";
    pass = false;
  }
}

template <class T>
void validate_pmx_solve_group_control_args(const T& ode_fun,
                                           const variable_map& var_map,
                                           bool& pass,
                                           std::ostream& error_msgs) {
  if (!ode_fun.rel_tol_.bare_type().is_primitive()) {
    error_msgs << "15th argument to "
               << ode_fun.integration_function_name_
               << " must be type real or int"
               << " for relative tolerance"
               << "; found type="
               << ode_fun.rel_tol_.bare_type()
               << ". ";
    pass = false;
  }
  if (!ode_fun.abs_tol_.bare_type().is_primitive()) {
    error_msgs << "16th argument to "
               << ode_fun.integration_function_name_
               << " must be type real or int"
               << " for absolute tolerance"
               << "; found type="
               << ode_fun.abs_tol_.bare_type()
               << ". ";
    pass = false;
  }
  if (!ode_fun.max_num_steps_.bare_type().is_primitive()) {
    error_msgs << "17th argument to "
               << ode_fun.integration_function_name_
               << " must be type real or int"
               << " for maximum number of steps"
               << "; found type="
               << ode_fun.max_num_steps_.bare_type()
               << ". ";
    pass = false;
  }

  // test data-only variables do not have parameters (int locals OK)
  if (has_var(ode_fun.rel_tol_, var_map)) {
    error_msgs << "15th argument to "
               << ode_fun.integration_function_name_
               << " for relative tolerance"
               << " must be data only and not reference parameters";
    pass = false;
  }
  if (has_var(ode_fun.abs_tol_, var_map)) {
    error_msgs << "16th argument to "
               << ode_fun.integration_function_name_
               << " for absolute tolerance"
               << " must be data only and not reference parameters";
    pass = false;
  }
  if (has_var(ode_fun.max_num_steps_, var_map)) {
    error_msgs << "17th argument to "
               << ode_fun.integration_function_name_
               << " for maximum number of steps"
               << " must be data only and not reference parameters";
    pass = false;
  }
}

/*********************************
  pmx_integrate_ode
 *********************************/
template <class T>
void validate_pmx_integrate_ode_non_control_args(const T& ode_fun,
                                                 const variable_map& var_map,
                                                 bool& pass,
                                                 std::ostream& error_msgs) {
  pass = true;
  // test function argument type
  torsten_types pmx_t(torsten_types::t_dbl_1);
  std::vector<bare_expr_type> sys_arg_types;
  sys_arg_types.push_back(torsten_types::t_dbl);
  sys_arg_types.push_back(torsten_types::t_dbl_1);
  sys_arg_types.push_back(torsten_types::t_dbl_1);
  sys_arg_types.push_back(torsten_types::t_dbl_1);
  sys_arg_types.push_back(torsten_types::t_int_1);
  function_signature_t system_signature(pmx_t.sys_result_type, sys_arg_types);
  if (!function_signatures::instance()
      .is_defined(ode_fun.system_function_name_, system_signature)) {
    error_msgs << "1st argument to "
               << ode_fun.integration_function_name_
               << " must be a function with signature"
               << " (real, real[], real[], real[], int[]) : real[] ";
    pass = false;
  }

  // test regular argument types
  if (ode_fun.y0_.bare_type() != torsten_types::t_dbl_1) {
    error_msgs << "2nd argument to "
               << ode_fun.integration_function_name_
               << " must have type real[] for intial system state;"
               << " found type="
               << ode_fun.y0_.bare_type()
               << ". ";
    pass = false;
  }
  if (!ode_fun.t0_.bare_type().is_primitive()) {
    error_msgs << "3rd argument to "
               << ode_fun.integration_function_name_
               << " must have type real or int for initial time;"
               << " found type="
               << ode_fun.t0_.bare_type()
               << ". ";
    pass = false;
  }
  if (ode_fun.ts_.bare_type() != torsten_types::t_dbl_1) {
    error_msgs << "4th argument to "
               << ode_fun.integration_function_name_
               << " must have type real[]"
               << " for requested solution times; found type="
               << ode_fun.ts_.bare_type()
               << ". ";
    pass = false;
  }
  if (ode_fun.theta_.bare_type() != torsten_types::t_dbl_1) {
    error_msgs << "5th argument to "
               << ode_fun.integration_function_name_
               << " must have type real[] for parameters; found type="
               << ode_fun.theta_.bare_type()
               << ". ";
    pass = false;
  }
  if (ode_fun.x_.bare_type() != torsten_types::t_dbl_1) {
    error_msgs << "6th argument to "
               << ode_fun.integration_function_name_
               << " must have type real[] for real data; found type="
               << ode_fun.x_.bare_type()
               << ". ";
    pass = false;
  }
  if (ode_fun.x_int_.bare_type() != torsten_types::t_int_1) {
    error_msgs << "7th argument to "
               << ode_fun.integration_function_name_
               << " must have type int[] for integer data; found type="
               << ode_fun.x_int_.bare_type()
               << ". ";
    pass = false;
  }

  // test data-only variables do not have parameters (int locals OK)
  if (has_var(ode_fun.t0_, var_map)) {
    error_msgs << "3rd argument to "
               << ode_fun.integration_function_name_
               << " (initial times)"
               << " must be data only and not reference parameters";
    pass = false;
  }
  // if (has_var(ode_fun.ts_, var_map)) {
  //   error_msgs << "4th argument to "
  //              << ode_fun.integration_function_name_
  //              << " (solution times)"
  //              << " must be data only and not reference parameters";
  //   pass = false;
  // }
  if (has_var(ode_fun.x_, var_map)) {
    error_msgs << "6th argument to "
               << ode_fun.integration_function_name_
               << " for real data"
               << " must be data only and not reference parameters";
    pass = false;
  }
}

void validate_pmx_integrate_ode::operator()(const pmx_integrate_ode& ode_fun,
                                            const variable_map& var_map,
                                            bool& pass,
                                            std::ostream& error_msgs) const {
  validate_pmx_integrate_ode_non_control_args(ode_fun, var_map, pass, error_msgs);
}
boost::phoenix::function<validate_pmx_integrate_ode>
validate_pmx_integrate_ode_f;

bool data_only_expression::operator()(const pmx_integrate_ode& x) const {
  return boost::apply_visitor(*this, x.y0_.expr_)
    && boost::apply_visitor(*this, x.theta_.expr_);
}

template void assign_lhs::operator()(expression&, const pmx_integrate_ode&) const;

/*********************************
  pmx_integrate_ode_control
 *********************************/
void validate_pmx_integrate_ode_control:: operator()(const pmx_integrate_ode_control &ode_fun,
                                                     const variable_map &var_map,
                                                     bool &pass, std::ostream &error_msgs) const {
  validate_pmx_integrate_ode_non_control_args(ode_fun, var_map, pass, error_msgs);
  if (!ode_fun.rel_tol_.bare_type().is_primitive()) {
    error_msgs << "Eighth argument to " << ode_fun.integration_function_name_
               << " (relative tolerance) must have type real or int;"
               << " found type=" << ode_fun.rel_tol_.bare_type() << ". ";
    pass = false;
  }
  if (!ode_fun.abs_tol_.bare_type().is_primitive()) {
    error_msgs << "Ninth argument to " << ode_fun.integration_function_name_
               << " (absolute tolerance) must have type real or int;"
               << " found type=" << ode_fun.abs_tol_.bare_type() << ". ";
    pass = false;
  }
  if (!ode_fun.max_num_steps_.bare_type().is_primitive()) {
    error_msgs << "Tenth argument to " << ode_fun.integration_function_name_
               << " (max steps) must have type real or int;"
               << " found type=" << ode_fun.max_num_steps_.bare_type() << ". ";
    pass = false;
  }

  // test data-only variables do not have parameters (int locals OK)
  if (has_var(ode_fun.rel_tol_, var_map)) {
    error_msgs << "Eighth argument to " << ode_fun.integration_function_name_
               << " (relative tolerance) must be data only"
               << " and not depend on parameters.";
    pass = false;
  }
  if (has_var(ode_fun.abs_tol_, var_map)) {
    error_msgs << "Ninth argument to " << ode_fun.integration_function_name_
               << " (absolute tolerance ) must be data only"
               << " and not depend parameters.";
    pass = false;
  }
  if (has_var(ode_fun.max_num_steps_, var_map)) {
    error_msgs << "Tenth argument to " << ode_fun.integration_function_name_
               << " (max steps) must be data only"
               << " and not depend on parameters.";
    pass = false;
  }
}
boost::phoenix::function<validate_pmx_integrate_ode_control> validate_pmx_integrate_ode_control_f;

bool data_only_expression::operator()(const pmx_integrate_ode_control& x) const {
  return boost::apply_visitor(*this, x.y0_.expr_)
    && boost::apply_visitor(*this, x.theta_.expr_);
}

template void assign_lhs::operator()(expression &, const pmx_integrate_ode_control &) const;

/*********************************
  pmx_integrate_ode_group
 *********************************/
    template <class T>
    void validate_pmx_integrate_ode_group_non_control_args(const T& ode_fun,
                                                           const variable_map& var_map,
                                                           bool& pass,
                                                           std::ostream& error_msgs) {
      pass = true;
      // test function argument type
      torsten_types pmx_t(torsten_types::t_dbl_1);
      std::vector<bare_expr_type> sys_arg_types;
      sys_arg_types.push_back(torsten_types::t_dbl);
      sys_arg_types.push_back(torsten_types::t_dbl_1);
      sys_arg_types.push_back(torsten_types::t_dbl_1);
      sys_arg_types.push_back(torsten_types::t_dbl_1);
      sys_arg_types.push_back(torsten_types::t_int_1);
      function_signature_t system_signature(pmx_t.sys_result_type, sys_arg_types);
      if (!function_signatures::instance()
          .is_defined(ode_fun.system_function_name_, system_signature)) {
        error_msgs << "1st argument to "
                   << ode_fun.integration_function_name_
                   << " must be a function with signature"
                   << " (real, real[], real[], real[], int[]) : real[] ";
        pass = false;
      }

      // test regular argument types
      if (ode_fun.y0_.bare_type() != torsten_types::t_dbl_2) {
        error_msgs << "2nd argument to "
                   << ode_fun.integration_function_name_
                   << " must have type real[ , ] for intial system state;"
                   << " found type="
                   << ode_fun.y0_.bare_type()
                   << ". ";
        pass = false;
      }
      if (!ode_fun.t0_.bare_type().is_primitive()) {
        error_msgs << "3rd argument to "
                   << ode_fun.integration_function_name_
                   << " must have type real or int for initial time;"
                   << " found type="
                   << ode_fun.t0_.bare_type()
                   << ". ";
        pass = false;
      }
      if (ode_fun.len_.bare_type() != torsten_types::t_int_1) {
        error_msgs << "4th argument to "
                   << ode_fun.integration_function_name_
                   << " must have type int[]"
                   << " for length of each ODE's times within ragged array; found type="
                   << ode_fun.len_.bare_type()
                   << ". ";
        pass = false;
      }
      if (ode_fun.ts_.bare_type() != torsten_types::t_dbl_1) {
        error_msgs << "5th argument to "
                   << ode_fun.integration_function_name_
                   << " must have type real[]"
                   << " for requested solution times; found type="
                   << ode_fun.ts_.bare_type()
                   << ". ";
        pass = false;
      }
      if (ode_fun.theta_.bare_type() != torsten_types::t_dbl_2) {
        error_msgs << "6th argument to "
                   << ode_fun.integration_function_name_
                   << " must have type real[ , ] for parameters; found type="
                   << ode_fun.theta_.bare_type()
                   << ". ";
        pass = false;
      }
      if (ode_fun.x_.bare_type() != torsten_types::t_dbl_2) {
        error_msgs << "7th argument to "
                   << ode_fun.integration_function_name_
                   << " must have type real[ , ] for real data; found type="
                   << ode_fun.x_.bare_type()
                   << ". ";
        pass = false;
      }
      if (ode_fun.x_int_.bare_type() != torsten_types::t_int_2) {
        error_msgs << "8th argument to "
                   << ode_fun.integration_function_name_
                   << " must have type int[ , ] for integer data; found type="
                   << ode_fun.x_int_.bare_type()
                   << ". ";
        pass = false;
      }

      // test data-only variables do not have parameters (int locals OK)
      if (has_var(ode_fun.t0_, var_map)) {
        error_msgs << "3rd argument to "
                   << ode_fun.integration_function_name_
                   << " (initial times)"
                   << " must be data only and not reference parameters";
        pass = false;
      }
      if (has_var(ode_fun.ts_, var_map)) {
        error_msgs << "5th argument to "
                   << ode_fun.integration_function_name_
                   << " (solution times)"
                   << " must be data only and not reference parameters";
        pass = false;
      }
      if (has_var(ode_fun.x_, var_map)) {
        error_msgs << "7th argument to "
                   << ode_fun.integration_function_name_
                   << " for real data"
                   << " must be data only and not reference parameters";
        pass = false;
      }
    }

/*****************
 generalOdeModel_control
*****************/

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

/*****************
 generalOdeModel
*****************/

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

/*****************
 pmx_solve_group
*****************/

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

/***********************
 pmx_solve_group_control
***********************/

void validate_pmx_solve_group_control::operator()(
                      const pmx_solve_group_control& ode_fun,
                      const variable_map& var_map,
                      bool& pass,
                      std::ostream& error_msgs) const {
  validate_pmx_solve_group_non_control_args(ode_fun, var_map, pass, error_msgs);
  validate_pmx_solve_group_control_args(ode_fun, var_map, pass, error_msgs);
}
boost::phoenix::function<validate_pmx_solve_group_control>
validate_pmx_solve_group_control_f;

bool data_only_expression::operator()(const pmx_solve_group_control& x)
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
                                     const pmx_solve_group_control&)
  const;

/*************************
 pmx_integrate_ode_group
*************************/

void validate_pmx_integrate_ode_group::operator()(
                      const pmx_integrate_ode_group& ode_fun,
                      const variable_map& var_map,
                      bool& pass,
                      std::ostream& error_msgs) const {
  validate_pmx_integrate_ode_group_non_control_args(ode_fun, var_map, pass, error_msgs);

  // collect ODE functor names to be used in MPI master-slave control
  pmx_integrate_ode_group::CALLED_FUNCTORS.push_back(ode_fun.system_function_name_);
}
boost::phoenix::function<validate_pmx_integrate_ode_group>
validate_pmx_integrate_ode_group_f;

bool data_only_expression::operator()(const pmx_integrate_ode_group& x) const {
  return boost::apply_visitor(*this, x.y0_.expr_)
    && boost::apply_visitor(*this, x.theta_.expr_);
}

template void assign_lhs::operator()(expression&, const pmx_integrate_ode_group&) const;

/*********************************
  pmx_integrate_ode_group_control
 *********************************/
void validate_pmx_integrate_ode_group_control:: operator()(const pmx_integrate_ode_group_control &ode_fun,
                                                     const variable_map &var_map,
                                                     bool &pass, std::ostream &error_msgs) const {
  validate_pmx_integrate_ode_group_non_control_args(ode_fun, var_map, pass, error_msgs);
  if (!ode_fun.rel_tol_.bare_type().is_primitive()) {
    error_msgs << "Eighth argument to " << ode_fun.integration_function_name_
               << " (relative tolerance) must have type real or int;"
               << " found type=" << ode_fun.rel_tol_.bare_type() << ". ";
    pass = false;
  }
  if (!ode_fun.abs_tol_.bare_type().is_primitive()) {
    error_msgs << "Ninth argument to " << ode_fun.integration_function_name_
               << " (absolute tolerance) must have type real or int;"
               << " found type=" << ode_fun.abs_tol_.bare_type() << ". ";
    pass = false;
  }
  if (!ode_fun.max_num_steps_.bare_type().is_primitive()) {
    error_msgs << "Tenth argument to " << ode_fun.integration_function_name_
               << " (max steps) must have type real or int;"
               << " found type=" << ode_fun.max_num_steps_.bare_type() << ". ";
    pass = false;
  }

  // test data-only variables do not have parameters (int locals OK)
  if (has_var(ode_fun.rel_tol_, var_map)) {
    error_msgs << "Eighth argument to " << ode_fun.integration_function_name_
               << " (relative tolerance) must be data only"
               << " and not depend on parameters.";
    pass = false;
  }
  if (has_var(ode_fun.abs_tol_, var_map)) {
    error_msgs << "Ninth argument to " << ode_fun.integration_function_name_
               << " (absolute tolerance ) must be data only"
               << " and not depend parameters.";
    pass = false;
  }
  if (has_var(ode_fun.max_num_steps_, var_map)) {
    error_msgs << "Tenth argument to " << ode_fun.integration_function_name_
               << " (max steps) must be data only"
               << " and not depend on parameters.";
    pass = false;
  }
}
boost::phoenix::function<validate_pmx_integrate_ode_group_control> validate_pmx_integrate_ode_group_control_f;

bool data_only_expression::operator()(const pmx_integrate_ode_group_control& x) const {
  return boost::apply_visitor(*this, x.y0_.expr_)
    && boost::apply_visitor(*this, x.theta_.expr_);
}

template void assign_lhs::operator()(expression&, const pmx_integrate_ode_group_control&) const;

#endif
