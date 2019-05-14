#ifndef STAN_LANG_TORSTEN_AST_NODE_EXPRESSION_DEF_HPP
#define STAN_LANG_TORSTEN_AST_NODE_EXPRESSION_DEF_HPP

// #include <stan/lang/ast.hpp>
#include <stan/torsten/torsten_func_expression_list.h>

namespace stan {
  namespace lang {

#define TORSTEN_FUNC_EXPR(F, R) expression::expression(const F& expr) : expr_(expr) { }
    TORSTEN_LANG_FUNCTORS_EXPRESSION_LIST
#undef TORSTEN_FUNC_EXPR

  }
}
#endif
