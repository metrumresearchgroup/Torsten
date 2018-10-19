#ifndef STAN_LANG_TORSTEN_AST_NODE_EXPRESSION_TYPE_VIS_HPP
#define STAN_LANG_TORSTEN_AST_NODE_EXPRESSION_TYPE_VIS_HPP

expr_type operator()(const univariate_integral_control& e) const;
expr_type operator()(const generalOdeModel_control& e) const;
expr_type operator()(const generalOdeModel& e) const;

#endif
