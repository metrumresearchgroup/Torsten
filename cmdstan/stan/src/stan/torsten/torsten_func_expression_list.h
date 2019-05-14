#ifndef STAN_LANG_TORSTEN_FUNC_EXPRESSION_LIST_H
#define STAN_LANG_TORSTEN_FUNC_EXPRESSION_LIST_H

/* list of Torsten's functions that accept functor argument */
/* (function name, return type) */
#define TORSTEN_LANG_FUNCTORS_EXPRESSION_LIST                                            \
  TORSTEN_FUNC_EXPR(univariate_integral_control     , bare_expr_type(double_type()))     \
  TORSTEN_FUNC_EXPR(generalOdeModel_control         , bare_expr_type(matrix_type()))     \
  TORSTEN_FUNC_EXPR(generalOdeModel                 , bare_expr_type(matrix_type()))     \
  TORSTEN_FUNC_EXPR(pmx_solve_group_control         , bare_expr_type(matrix_type()))     \
  TORSTEN_FUNC_EXPR(pmx_solve_group                 , bare_expr_type(matrix_type()))     \
  TORSTEN_FUNC_EXPR(pmx_integrate_ode_control       , bare_array_type(double_type(), 2)) \
  TORSTEN_FUNC_EXPR(pmx_integrate_ode               , bare_array_type(double_type(), 2)) \
  TORSTEN_FUNC_EXPR(pmx_integrate_ode_group_control , bare_expr_type(matrix_type()))     \
  TORSTEN_FUNC_EXPR(pmx_integrate_ode_group         , bare_expr_type(matrix_type()))

#endif
