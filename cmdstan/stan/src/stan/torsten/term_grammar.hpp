#ifndef STAN_LANG_TORSTEN_GRAMMARS_TERM_GRAMMAR_HPP
#define STAN_LANG_TORSTEN_GRAMMARS_TERM_GRAMMAR_HPP

boost::spirit::qi::rule<Iterator,
                        univariate_integral_control(scope),
                        whitespace_grammar<Iterator> >
univariate_integral_control_r;

boost::spirit::qi::rule<Iterator,
                        generalOdeModel_control(scope),
                        whitespace_grammar<Iterator> >
generalOdeModel_control_r;

boost::spirit::qi::rule<Iterator,
                        generalOdeModel(scope),
                        whitespace_grammar<Iterator> >
generalOdeModel_r;


#endif
