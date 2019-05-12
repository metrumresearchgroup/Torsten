#ifndef STAN_LANG_TORSTEN_AST_FUN_IS_NIL_VIS_DEF_HPP
#define STAN_LANG_TORSTEN_AST_FUN_IS_NIL_VIS_DEF_HPP

#include <stan/torsten/torsten_func_expression_list.h>

namespace stan {
  namespace lang {

#define TORSTEN_FUNC_EXPR(F, R) bool is_nil_vis::operator()(const F& /* x */) const {return false;}
    TORSTEN_LANG_FUNCTORS_EXPRESSION_LIST
#undef TORSTEN_FUNC_EXPR

  }
}

#endif
