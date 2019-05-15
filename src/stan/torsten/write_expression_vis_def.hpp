#ifndef STAN_LANG_TORSTEN_AST_FUN_WRITE_EXPRESSION_VIS_DEF_HPP
#define STAN_LANG_TORSTEN_AST_FUN_WRITE_EXPRESSION_VIS_DEF_HPP

#include <stan/lang/ast.hpp>
#include <boost/lexical_cast.hpp>
#include <sstream>
#include <string>
#include <vector>

namespace stan {
namespace lang {

std::string write_expression_vis::operator()(const univariate_integral_control& e) const {
  std::stringstream ss;
  ss << e.integration_function_name_ << "(" << e.system_function_name_ << ", "
     << e.t0_.to_string() << ", " << e.t1_.to_string() << ", "
     << e.x_r_.to_string() << ", " << e.x_i_.to_string() << ")";
  return ss.str();
}

std::string write_expression_vis::operator()(const pmx_integrate_ode& e) const {
  std::stringstream ss;
  ss << e.integration_function_name_ << "(" << e.system_function_name_ << ", "
     << e.y0_.to_string() << ", " << e.t0_.to_string() << ", "
     << e.ts_.to_string() << ", " << e.x_.to_string() << ", "
     << e.x_int_.to_string() << ")";
  return ss.str();
}

std::string write_expression_vis::operator()(const pmx_integrate_ode_control& e) const {
  std::stringstream ss;
  ss << e.integration_function_name_ << "(" << e.system_function_name_ << ", "
     << e.y0_.to_string() << ", " << e.t0_.to_string() << ", "
     << e.ts_.to_string() << ", " << e.x_.to_string() << ", "
     << e.x_int_.to_string() << ", " << e.rel_tol_.to_string() << ", "
     << e.abs_tol_.to_string() << ", " << e.max_num_steps_.to_string() << ")";
  return ss.str();
}

std::string write_expression_vis::operator()(const pmx_integrate_ode_group& e) const {
  std::stringstream ss;
  ss << e.integration_function_name_ << "(" << e.system_function_name_ << ", "
     << e.y0_.to_string() << ", " << e.t0_.to_string() << ", "
     << e.len_.to_string() << ", "
     << e.ts_.to_string() << ", " << e.x_.to_string() << ", "
     << e.x_int_.to_string() << ")";
  return ss.str();
}

std::string write_expression_vis::operator()(const pmx_integrate_ode_group_control& e) const {
  std::stringstream ss;
  ss << e.integration_function_name_ << "(" << e.system_function_name_ << ", "
     << e.y0_.to_string() << ", " << e.t0_.to_string() << ", "
     << e.len_.to_string() << ", "
     << e.ts_.to_string() << ", " << e.x_.to_string() << ", "
     << e.x_int_.to_string() << ", " << e.rel_tol_.to_string() << ", "
     << e.abs_tol_.to_string() << ", " << e.max_num_steps_.to_string() << ")";
  return ss.str();
}

std::string write_expression_vis::operator()(const generalOdeModel& e) const {
  std::stringstream ss;
  ss << e.integration_function_name_ << "(" << e.system_function_name_ << ", "
     << e.nCmt_.to_string() << ", "
     << e.time_.to_string() << ", "
     << e.amt_.to_string() << ", "
     << e.rate_.to_string() << ", "
     << e.ii_.to_string() << ", "
     << e.evid_.to_string() << ", "
     << e.cmt_.to_string() << ", "
     << e.addl_.to_string() << ", "
     << e.ss_.to_string() << ", "
     << e.pMatrix_.to_string() << ", "
     << e.biovar_.to_string() << ", "
     << e.tlag_.to_string() << ")";
  return ss.str();
}

std::string write_expression_vis::operator()(const generalOdeModel_control& e) const {
  std::stringstream ss;
  ss << e.integration_function_name_ << "(" << e.system_function_name_ << ", "
     << e.nCmt_.to_string() << ", "
     << e.time_.to_string() << ", "
     << e.amt_.to_string() << ", "
     << e.rate_.to_string() << ", "
     << e.ii_.to_string() << ", "
     << e.evid_.to_string() << ", "
     << e.cmt_.to_string() << ", "
     << e.addl_.to_string() << ", "
     << e.ss_.to_string() << ", "
     << e.pMatrix_.to_string() << ", "
     << e.biovar_.to_string() << ", "
     << e.tlag_.to_string() << ", "
     << e.rel_tol_.to_string() << ", "
     << e.abs_tol_.to_string() << ", "
     << e.max_num_steps_.to_string() << ")";
  return ss.str();
}

std::string write_expression_vis::operator()(const pmx_solve_group& e) const {
  std::stringstream ss;
  ss << e.integration_function_name_ << "(" << e.system_function_name_ << ", "
     << e.nCmt_.to_string() << ", "
     << e.len_.to_string() << ", "
     << e.time_.to_string() << ", "
     << e.amt_.to_string() << ", "
     << e.rate_.to_string() << ", "
     << e.ii_.to_string() << ", "
     << e.evid_.to_string() << ", "
     << e.cmt_.to_string() << ", "
     << e.addl_.to_string() << ", "
     << e.ss_.to_string() << ", "
     << e.pMatrix_.to_string() << ", "
     << e.biovar_.to_string() << ", "
     << e.tlag_.to_string() << ")";
  return ss.str();
}

std::string write_expression_vis::operator()(const pmx_solve_group_control& e) const {
  std::stringstream ss;
  ss << e.integration_function_name_ << "(" << e.system_function_name_ << ", "
     << e.nCmt_.to_string() << ", "
     << e.len_.to_string() << ", "
     << e.time_.to_string() << ", "
     << e.amt_.to_string() << ", "
     << e.rate_.to_string() << ", "
     << e.ii_.to_string() << ", "
     << e.evid_.to_string() << ", "
     << e.cmt_.to_string() << ", "
     << e.addl_.to_string() << ", "
     << e.ss_.to_string() << ", "
     << e.pMatrix_.to_string() << ", "
     << e.biovar_.to_string() << ", "
     << e.tlag_.to_string() << ", "
     << e.rel_tol_.to_string() << ", "
     << e.abs_tol_.to_string() << ", "
     << e.max_num_steps_.to_string() << ")";
  return ss.str();
}

}  // namespace lang
}  // namespace stan
#endif
