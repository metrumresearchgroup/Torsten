#ifndef STAN_LANG_TORSTEN_AST_FUN_VAR_OCCURS_VIS_DEF_HPP
#define STAN_LANG_TORSTEN_AST_FUN_VAR_OCCURS_VIS_DEF_HPP

#include <stan/lang/ast.hpp>
#include <stan/torsten/torsten_func_expression_list.h>

namespace stan {
  namespace lang {

// no refs persist out of these calls
#define TORSTEN_FUNC_EXPR(F, R) bool var_occurs_vis::operator()(const F& e) const {return false;}
    TORSTEN_LANG_FUNCTORS_EXPRESSION_LIST
#undef TORSTEN_FUNC_EXPR

  }
}

#endif
