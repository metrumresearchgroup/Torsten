#ifndef STAN_LANG_TORSTEN_AST_NODE_EXPRESSION_BARE_TYPE_VIS_HPP
#define STAN_LANG_TORSTEN_AST_NODE_EXPRESSION_BARE_TYPE_VIS_HPP

#include <stan/torsten/torsten_func_expression_list.h>

#define TORSTEN_FUNC_EXPR(F, R) bare_expr_type operator()(const F& e) const;
    TORSTEN_LANG_FUNCTORS_EXPRESSION_LIST
#undef TORSTEN_FUNC_EXPR

#endif
