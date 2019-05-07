#ifndef STAN_LANG_TORSTEN_FUNC_EXPRESSION_LIST_H
#define STAN_LANG_TORSTEN_FUNC_EXPRESSION_LIST_H

/* list of Torsten's functions that accept functor argument */
/* (function name, return type) */
#define TORSTEN_LANG_FUNCTORS_EXPRESSION_LIST                                 \
  TORSTEN_FUNC_EXPR(univariate_integral_control, expr_type(double_type())   ) \
  TORSTEN_FUNC_EXPR(generalOdeModel_control    , expr_type(matrix_type(), 0)) \
  TORSTEN_FUNC_EXPR(generalOdeModel            , expr_type(matrix_type(), 0)) \
  TORSTEN_FUNC_EXPR(pmx_solve_group            , expr_type(matrix_type(), 0)) \
  TORSTEN_FUNC_EXPR(pmx_solve_group_control    , expr_type(matrix_type(), 0)) \
  TORSTEN_FUNC_EXPR(pmx_integrate_ode          , expr_type(double_type(), 2)) \
  TORSTEN_FUNC_EXPR(pmx_integrate_ode_group    , expr_type(matrix_type(), 0))

#endif
