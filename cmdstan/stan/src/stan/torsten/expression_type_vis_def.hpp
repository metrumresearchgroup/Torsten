#ifndef STAN_LANG_TORSTEN_AST_NODE_EXPRESSION_TYPE_VIS_DEF_HPP
#define STAN_LANG_TORSTEN_AST_NODE_EXPRESSION_TYPE_VIS_DEF_HPP

#include <stan/lang/ast.hpp>

namespace stan {
  namespace lang {

    expr_type
    expression_type_vis::operator()(const univariate_integral_control& e) const {
      return expr_type(double_type());
    }
    expr_type
    expression_type_vis::operator()(const generalOdeModel_control& e) const {
      return expr_type(matrix_type(), 0);
    }
    expr_type
    expression_type_vis::operator()(const generalOdeModel& e) const {
      return expr_type(matrix_type(), 0);
    }

  }
}

#endif
