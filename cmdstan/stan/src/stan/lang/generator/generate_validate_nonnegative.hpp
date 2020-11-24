#ifndef STAN_LANG_GENERATOR_GENERATE_VALIDATE_NONNEGATIVE_HPP
#define STAN_LANG_GENERATOR_GENERATE_VALIDATE_NONNEGATIVE_HPP

#include <stan/lang/ast.hpp>
#include <stan/lang/generator/constants.hpp>
#include <stan/lang/generator/generate_indent.hpp>
#include <stan/lang/generator/generate_quoted_expression.hpp>
#include <ostream>
#include <string>

namespace stan {
namespace lang {

/**
 * Generate call to stan_math lib function validate_non_negative_index
 * which will throw an informative error if dim size is < 0
 *
 * This check should precede the variable declaration in order to
 * avoid bad alloc runtime error.
 * Called by
 * <br />generate_validate_var_dims
 * <br/> generate_set_param_ranges
 *

 * @param[in] name variable name
 * @param[in] expr dim size expression
 * @param[in] indent indentation level
 * @param[in,out] o output stream for generated code
 */
void generate_validate_nonnegative(const std::string& name,
                                   const expression& expr, int indent,
                                   std::ostream& o) {
  generate_indent(indent, o);
  o << "validate_non_negative_index(\"" << name << "\", ";
  generate_quoted_expression(expr, o);
  o << ", ";
  generate_expression(expr, NOT_USER_FACING, o);
  o << ");" << EOL;
}
}  // namespace lang
}  // namespace stan
#endif
