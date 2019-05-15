#ifndef STAN_LANG_TORSTEN_GENERATOR_EXPRESSION_VISGEN_HPP
#define STAN_LANG_TORSTEN_GENERATOR_EXPRESSION_VISGEN_HPP

void operator()(const univariate_integral_control& fx) const {
  o_ << fx.integration_function_name_
     << '('
     << fx.system_function_name_
     << "_functor__(), ";
  generate_expression(fx.t0_, user_facing_, o_);
  o_ << ", ";
  generate_expression(fx.t1_, user_facing_, o_);
  o_ << ", ";
  generate_expression(fx.theta_, user_facing_, o_);
  o_ << ", ";
  generate_expression(fx.x_r_, NOT_USER_FACING, o_);
  o_ << ", ";
  generate_expression(fx.x_i_, NOT_USER_FACING, o_);
  o_ << ", pstream__)";
}

void operator()(const generalOdeModel_control& fx) const {
  o_ << fx.integration_function_name_
     << '('
     << fx.system_function_name_
     << "_functor__(), ";

  generate_expression(fx.nCmt_, NOT_USER_FACING, o_);
  o_ << ", ";

  generate_expression(fx.time_, NOT_USER_FACING, o_);
  o_ << ", ";

  generate_expression(fx.amt_, NOT_USER_FACING, o_);
  o_ << ", ";

  generate_expression(fx.rate_, NOT_USER_FACING, o_);
  o_ << ", ";

  generate_expression(fx.ii_, NOT_USER_FACING, o_);
  o_ << ", ";

  generate_expression(fx.evid_, NOT_USER_FACING, o_);
  o_ << ", ";

  generate_expression(fx.cmt_, NOT_USER_FACING, o_);
  o_ << ", ";

  generate_expression(fx.addl_, NOT_USER_FACING, o_);
  o_ << ", ";

  generate_expression(fx.ss_, NOT_USER_FACING, o_);
  o_ << ", ";

  generate_expression(fx.pMatrix_, NOT_USER_FACING, o_);
  o_ << ", ";

  generate_expression(fx.biovar_, NOT_USER_FACING, o_);
  o_ << ", ";

  generate_expression(fx.tlag_, NOT_USER_FACING, o_);
  o_ << ", pstream__, ";

  generate_expression(fx.rel_tol_, NOT_USER_FACING, o_);
  o_ << ", ";

  generate_expression(fx.abs_tol_, NOT_USER_FACING, o_);
  o_ << ", ";

  generate_expression(fx.max_num_steps_, NOT_USER_FACING, o_);
  o_ << ")";
}

void operator()(const generalOdeModel& fx) const {
  o_ << fx.integration_function_name_
     << '('
     << fx.system_function_name_
     << "_functor__(), ";

  generate_expression(fx.nCmt_, NOT_USER_FACING, o_);
  o_ << ", ";

  generate_expression(fx.time_, NOT_USER_FACING, o_);
  o_ << ", ";

  generate_expression(fx.amt_, NOT_USER_FACING, o_);
  o_ << ", ";

  generate_expression(fx.rate_, NOT_USER_FACING, o_);
  o_ << ", ";

  generate_expression(fx.ii_, NOT_USER_FACING, o_);
  o_ << ", ";

  generate_expression(fx.evid_, NOT_USER_FACING, o_);
  o_ << ", ";

  generate_expression(fx.cmt_, NOT_USER_FACING, o_);
  o_ << ", ";

  generate_expression(fx.addl_, NOT_USER_FACING, o_);
  o_ << ", ";

  generate_expression(fx.ss_, NOT_USER_FACING, o_);
  o_ << ", ";

  generate_expression(fx.pMatrix_, NOT_USER_FACING, o_);
  o_ << ", ";

  generate_expression(fx.biovar_, NOT_USER_FACING, o_);
  o_ << ", ";

  generate_expression(fx.tlag_, NOT_USER_FACING, o_);
  o_ << ", pstream__)";
}

void operator()(const pmx_integrate_ode& fx) const {
  o_ << fx.integration_function_name_
     << '('
     << fx.system_function_name_
     << "_functor__(), ";
  generate_expression(fx.y0_, NOT_USER_FACING, o_);
  o_ << ", ";
  generate_expression(fx.t0_, NOT_USER_FACING, o_);
  o_ << ", ";
  generate_expression(fx.ts_, NOT_USER_FACING, o_);
  o_ << ", ";
  generate_expression(fx.theta_, user_facing_, o_);
  o_ << ", ";
  generate_expression(fx.x_, NOT_USER_FACING, o_);
  o_ << ", ";
  generate_expression(fx.x_int_, NOT_USER_FACING, o_);
  o_ << ", pstream__)";
}

void operator()(const pmx_integrate_ode_control& fx) const {
  o_ << fx.integration_function_name_
     << '('
     << fx.system_function_name_
     << "_functor__(), ";
  generate_expression(fx.y0_, NOT_USER_FACING, o_);
  o_ << ", ";
  generate_expression(fx.t0_, NOT_USER_FACING, o_);
  o_ << ", ";
  generate_expression(fx.ts_, NOT_USER_FACING, o_);
  o_ << ", ";
  generate_expression(fx.theta_, user_facing_, o_);
  o_ << ", ";
  generate_expression(fx.x_, NOT_USER_FACING, o_);
  o_ << ", ";
  generate_expression(fx.x_int_, NOT_USER_FACING, o_);
  o_ << ", pstream__, ";
  generate_expression(fx.rel_tol_, NOT_USER_FACING, o_);
  o_ << ", ";
  generate_expression(fx.abs_tol_, NOT_USER_FACING, o_);
  o_ << ", ";
  generate_expression(fx.max_num_steps_, NOT_USER_FACING, o_);
  o_ << ")";
}

void operator()(const pmx_integrate_ode_group& fx) const {
  o_ << fx.integration_function_name_
     << '('
     << fx.system_function_name_
     << "_functor__(), ";
  generate_expression(fx.y0_, NOT_USER_FACING, o_);
  o_ << ", ";
  generate_expression(fx.t0_, NOT_USER_FACING, o_);
  o_ << ", ";
  generate_expression(fx.len_, NOT_USER_FACING, o_);
  o_ << ", ";
  generate_expression(fx.ts_, NOT_USER_FACING, o_);
  o_ << ", ";
  generate_expression(fx.theta_, user_facing_, o_);
  o_ << ", ";
  generate_expression(fx.x_, NOT_USER_FACING, o_);
  o_ << ", ";
  generate_expression(fx.x_int_, NOT_USER_FACING, o_);
  o_ << ", pstream__)";
}

void operator()(const pmx_integrate_ode_group_control& fx) const {
  o_ << fx.integration_function_name_
     << '('
     << fx.system_function_name_
     << "_functor__(), ";
  generate_expression(fx.y0_, NOT_USER_FACING, o_);
  o_ << ", ";
  generate_expression(fx.t0_, NOT_USER_FACING, o_);
  o_ << ", ";
  generate_expression(fx.len_, NOT_USER_FACING, o_);
  o_ << ", ";
  generate_expression(fx.ts_, NOT_USER_FACING, o_);
  o_ << ", ";
  generate_expression(fx.theta_, user_facing_, o_);
  o_ << ", ";
  generate_expression(fx.x_, NOT_USER_FACING, o_);
  o_ << ", ";
  generate_expression(fx.x_int_, NOT_USER_FACING, o_);
  o_ << ", pstream__, ";
  generate_expression(fx.rel_tol_, NOT_USER_FACING, o_);
  o_ << ", ";
  generate_expression(fx.abs_tol_, NOT_USER_FACING, o_);
  o_ << ", ";
  generate_expression(fx.max_num_steps_, NOT_USER_FACING, o_);
  o_ << ")";
}

void operator()(const pmx_solve_group& fx) const {
  o_ << fx.integration_function_name_
     << '('
     << fx.system_function_name_
     << "_functor__(), ";

  generate_expression(fx.nCmt_, NOT_USER_FACING, o_);
  o_ << ", ";

  generate_expression(fx.len_, NOT_USER_FACING, o_);
  o_ << ", ";

  generate_expression(fx.time_, NOT_USER_FACING, o_);
  o_ << ", ";

  generate_expression(fx.amt_, NOT_USER_FACING, o_);
  o_ << ", ";

  generate_expression(fx.rate_, NOT_USER_FACING, o_);
  o_ << ", ";

  generate_expression(fx.ii_, NOT_USER_FACING, o_);
  o_ << ", ";

  generate_expression(fx.evid_, NOT_USER_FACING, o_);
  o_ << ", ";

  generate_expression(fx.cmt_, NOT_USER_FACING, o_);
  o_ << ", ";

  generate_expression(fx.addl_, NOT_USER_FACING, o_);
  o_ << ", ";

  generate_expression(fx.ss_, NOT_USER_FACING, o_);
  o_ << ", ";

  generate_expression(fx.pMatrix_, NOT_USER_FACING, o_);
  o_ << ", ";

  generate_expression(fx.biovar_, NOT_USER_FACING, o_);
  o_ << ", ";

  generate_expression(fx.tlag_, NOT_USER_FACING, o_);
  o_ << ", pstream__)";
}

void operator()(const pmx_solve_group_control& fx) const {
  o_ << fx.integration_function_name_
     << '('
     << fx.system_function_name_
     << "_functor__(), ";

  generate_expression(fx.nCmt_, NOT_USER_FACING, o_);
  o_ << ", ";

  generate_expression(fx.len_, NOT_USER_FACING, o_);
  o_ << ", ";

  generate_expression(fx.time_, NOT_USER_FACING, o_);
  o_ << ", ";

  generate_expression(fx.amt_, NOT_USER_FACING, o_);
  o_ << ", ";

  generate_expression(fx.rate_, NOT_USER_FACING, o_);
  o_ << ", ";

  generate_expression(fx.ii_, NOT_USER_FACING, o_);
  o_ << ", ";

  generate_expression(fx.evid_, NOT_USER_FACING, o_);
  o_ << ", ";

  generate_expression(fx.cmt_, NOT_USER_FACING, o_);
  o_ << ", ";

  generate_expression(fx.addl_, NOT_USER_FACING, o_);
  o_ << ", ";

  generate_expression(fx.ss_, NOT_USER_FACING, o_);
  o_ << ", ";

  generate_expression(fx.pMatrix_, NOT_USER_FACING, o_);
  o_ << ", ";

  generate_expression(fx.biovar_, NOT_USER_FACING, o_);
  o_ << ", ";

  generate_expression(fx.tlag_, NOT_USER_FACING, o_);
  o_ << ", pstream__, ";

  generate_expression(fx.rel_tol_, NOT_USER_FACING, o_);
  o_ << ", ";

  generate_expression(fx.abs_tol_, NOT_USER_FACING, o_);
  o_ << ", ";

  generate_expression(fx.max_num_steps_, NOT_USER_FACING, o_);
  o_ << ")";
}
#endif
