#ifndef STAN_LANG_GENERATOR_GENERATE_TORSTEN_MPI_HPP
#define STAN_LANG_GENERATOR_GENERATE_TORSTEN_MPI_HPP

#include <stan/lang/ast.hpp>
#include <string>
#include <iostream>

namespace stan {
namespace lang {

  /*
   * In MPI calls that employ master-slave pattern, slaves
   * don't have functor infomation as they are locked in a
   * waiting-loop. So we delcare @c pmx_ode_group_mpi_functor 
   * in @c stan::math but define its operation here. We also
   * use perfect forwarding so the generator can be applied
   * to various functors.
   */
void generate_torsten_mpi(const std::string& model_name,
                           std::ostream& o) {
  // if (pmx_integrate_ode_group::CALLED_FUNCTORS.empty()) return;

  const std::set<std::string> funcs(pmx_integrate_ode_group::CALLED_FUNCTORS.begin(), // NOLINT
                                    pmx_integrate_ode_group::CALLED_FUNCTORS.end()); // NOLINT

  // generate @c pmx_ode_group_mpi_functor::operator()
  o << "namespace torsten {" << EOL
    << "namespace dsolve {" << EOL
    << INDENT << "template<typename... Args>" << EOL
    << INDENT << "inline auto pmx_ode_group_mpi_functor::operator()(Args&&... args) const {" << EOL; // NOLINT

  int i = 0;
  for (auto fs : funcs) {
    o << INDENT2 << "if (id == " << i << ") { "
      << model_name << "_namespace::" << fs << "_functor__ f; return f(std::forward<Args>(args)...); }" << EOL; // NOLINT
    i++;
  }
  // default to dummy
  o << INDENT2 << "dummy_functor f; return f(std::forward<Args>(args)...);" << EOL // NOLINT
    << INDENT << "}" << EOL2;

  // generate specializations of @c pmx_ode_group_mpi_functor_id
  i = 0;
  for (auto fs : funcs) {
    o << INDENT << "template<>" << EOL
      << INDENT << "struct pmx_ode_group_mpi_functor_id<"
      << model_name << "_namespace::" << fs << "_functor__> { static constexpr int value = " << i << "; };" // NOLINT
      << EOL;
    i++;
  }

  o << "}" << EOL << "}" << EOL;
}

}
}
#endif
