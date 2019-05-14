#ifndef STAN_LANG_TORSTEN_AST_NODE_EXPRESSION_TYPE_VIS_DEF_HPP
#define STAN_LANG_TORSTEN_AST_NODE_EXPRESSION_TYPE_VIS_DEF_HPP

#include <stan/lang/ast.hpp>
#include <stan/torsten/torsten_func_expression_list.h>

namespace stan {
  namespace lang {

#define TORSTEN_FUNC_EXPR(F, R) expr_type expression_type_vis::operator()(const F& e) const {return R;}
    TORSTEN_LANG_FUNCTORS_EXPRESSION_LIST
#undef TORSTEN_FUNC_EXPR

  }
}

#endif
