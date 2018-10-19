#ifndef STAN_LANG_TORSTEN_AST_NODE_EXPRESSION_DEF_HPP
#define STAN_LANG_TORSTEN_AST_NODE_EXPRESSION_DEF_HPP

#include <stan/lang/ast.hpp>

namespace stan {
  namespace lang {

    expression::expression(const univariate_integral_control& expr) :
      expr_(expr) { }
    expression::expression(const generalOdeModel_control& expr) :
      expr_(expr) { }
    expression::expression(const generalOdeModel& expr) :
      expr_(expr) { }
  }
}
#endif
