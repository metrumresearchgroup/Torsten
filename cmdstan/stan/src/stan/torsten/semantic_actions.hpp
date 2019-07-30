#ifndef STAN_LANG_TORSTEN_GRAMMARS_SEMANTIC_ACTIONS_HPP
#define STAN_LANG_TORSTEN_GRAMMARS_SEMANTIC_ACTIONS_HPP

#include <stan/torsten/torsten_func_expression_list.h>

// use macro to declare the following given torsten expresion F
// struct validate_F : public phoenix_functor_quaternary {
//   void operator()(const F& func,
//                   const variable_map& var_map, bool& pass,
//                   std::ostream& error_msgs) const;
// };
// extern boost::phoenix::function<validate_F>
// validate_F_f;

// called from: term_grammar
#define TORSTEN_FUNC_EXPR(F, R) struct validate_##F : public phoenix_functor_quaternary {void operator()(const F& func, const variable_map& var_map, bool& pass, std::ostream& error_msgs) const;}; 
  TORSTEN_LANG_FUNCTORS_EXPRESSION_LIST
#undef TORSTEN_FUNC_EXPR

#define TORSTEN_FUNC_EXPR(F, R) extern boost::phoenix::function<validate_##F> validate_##F##_f; 
  TORSTEN_LANG_FUNCTORS_EXPRESSION_LIST
#undef TORSTEN_FUNC_EXPR

#endif
