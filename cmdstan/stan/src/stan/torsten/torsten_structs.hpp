#ifndef STAN_LANG_TORSTEN_STRUCTS_HPP
#define STAN_LANG_TORSTEN_STRUCTS_HPP

#include <stan/torsten/torsten_func_expression_list.h>

namespace stan {
  namespace lang {

#define TORSTEN_FUNC_EXPR(F, R) struct F; 
    TORSTEN_LANG_FUNCTORS_EXPRESSION_LIST
#undef TORSTEN_FUNC_EXPR

 }
}

#endif
