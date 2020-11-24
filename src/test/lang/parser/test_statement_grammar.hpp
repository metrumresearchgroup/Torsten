#ifndef STAN_LANG_GRAMMARS_TEST_STATEMENT_GRAMMAR_HPP
#define STAN_LANG_GRAMMARS_TEST_STATEMENT_GRAMMAR_HPP

#include <stan/io/program_reader.hpp>
#include <stan/lang/ast.hpp>
#include <stan/lang/grammars/whitespace_grammar.hpp>
#include <stan/lang/grammars/expression_grammar.hpp>
#include <stan/lang/grammars/statement_grammar.hpp>
#include <stan/lang/grammars/semantic_actions.hpp>
#include <boost/spirit/include/qi.hpp>
#include <string>
#include <sstream>
#include <vector>
#include <utility>

namespace stan {

namespace lang {

template <typename Iterator>
struct test_statement_grammar
    : boost::spirit::qi::grammar<Iterator, boost::spirit::qi::locals<scope>,
                                 statement, whitespace_grammar<Iterator> > {
  const io::program_reader& reader_;
  variable_map var_map_;
  std::stringstream& error_msgs_;
  statement_grammar<Iterator> statement_g;

  test_statement_grammar(const io::program_reader& reader,
                         variable_map& var_map, std::stringstream& error_msgs);

  boost::spirit::qi::rule<Iterator, boost::spirit::qi::locals<scope>, statement,
                          whitespace_grammar<Iterator> >
      test_statement_r;
};

}  // namespace lang
}  // namespace stan
#endif
