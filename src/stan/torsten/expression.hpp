#ifndef STAN_LANG_TORSTEN_AST_NODE_EXPRESSION_HPP
#define STAN_LANG_TORSTEN_AST_NODE_EXPRESSION_HPP

#include <stan/torsten/torsten_func_expression_list.h>

#define TORSTEN_FUNC_EXPR(F, R) expression(const F& expr); 
    TORSTEN_LANG_FUNCTORS_EXPRESSION_LIST
#undef TORSTEN_FUNC_EXPR

#endif
