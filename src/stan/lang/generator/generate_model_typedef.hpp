#ifndef STAN_LANG_GENERATOR_GENERATE_MODEL_TYPEDEF_HPP
#define STAN_LANG_GENERATOR_GENERATE_MODEL_TYPEDEF_HPP

#include <stan/lang/ast.hpp>
#include <stan/lang/generator/constants.hpp>
#include <ostream>
#include <string>

namespace stan {
namespace lang {

/**
 * Generate reusable typedef of <code>stan_model</code> for
 * specified model name writing to the specified stream along with
 * a factory method `new_model()` to create a reference to a base model.
 *
 * @param model_name name of model
 * @param o stream for generating
 */
void generate_model_typedef(const std::string& model_name, std::ostream& o) {
  o << "typedef " << model_name << "_namespace::" << model_name
    << " stan_model;" << EOL2;

  o << "#ifndef USING_R" << EOL2;

  o << "stan::model::model_base& new_model(" << EOL
    << "        stan::io::var_context& data_context," << EOL
    << "        unsigned int seed," << EOL
    << "        std::ostream* msg_stream) {" << EOL
    << "  stan_model* m = new stan_model(data_context, seed, msg_stream);"
    << EOL << "  return *m;" << EOL << "}" << EOL2;

  o << "#endif" << EOL2;
}
}  // namespace lang
}  // namespace stan
#endif
