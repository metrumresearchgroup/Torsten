#ifndef STAN_LANG_TORSTEN_GRAMMARS_TERM_GRAMMAR_HPP
#define STAN_LANG_TORSTEN_GRAMMARS_TERM_GRAMMAR_HPP

#include <stan/torsten/torsten_func_expression_list.h>

/* use macro to declare the following given torsten expresion F */
// boost::spirit::qi::rule<Iterator,
//                         F(scope),
//                         whitespace_grammar<Iterator> >
// F_r;

#define TORSTEN_FUNC_EXPR(F, R) boost::spirit::qi::rule<Iterator, F(scope), whitespace_grammar<Iterator> > F##_r; 
  TORSTEN_LANG_FUNCTORS_EXPRESSION_LIST
#undef TORSTEN_FUNC_EXPR

#endif
