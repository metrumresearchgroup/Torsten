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

  const std::set<std::string> funcs(pmx_integrate_ode_group::CALLED_FUNCTORS.begin(),
                                    pmx_integrate_ode_group::CALLED_FUNCTORS.end());

  // generate @c pmx_ode_group_mpi_functor::operator()
  o << "namespace torsten {" << EOL
    << "namespace dsolve {" << EOL
    << INDENT << "template<typename... Args>" << EOL
    << INDENT << "inline auto pmx_ode_group_mpi_functor::operator()(Args&&... args) const {" << EOL;
    
  int i = 0;
  for (auto fs : funcs) {
    o << INDENT2 << "if (id == " << i << ") { " 
      << model_name << "_namespace::" << fs << "_functor__ f; return f(std::forward<Args>(args)...); }" << EOL;
    i++;
  }
  // default to dummy
  o << INDENT2 << "dummy_functor f; return f(std::forward<Args>(args)...);" << EOL
    << INDENT << "}" << EOL2;

  // generate specializations of @c pmx_ode_group_mpi_functor_id
  i = 0;
  for (auto fs : funcs) {
    o << INDENT << "template<>" << EOL
      << INDENT << "struct pmx_ode_group_mpi_functor_id<"
      << model_name << "_namespace::" << fs << "_functor__> { static constexpr int value = " << i << "; };"
      << EOL;
    i++;
  }

  o << "}" << EOL << "}" << EOL;
}

// void generate_torsten_mpi(const std::string& model_name,
//                            std::ostream& o) {
//   o << "struct pmx_ode_group_mpi_functor {" << std::endl
//     << "int id;" << std::endl << std::endl
//     << "pmx_ode_group_mpi_functor(int i) : id(i) {}" << std::endl << std::endl
//     << "template<typename T0, typename T1, typename T2>" << std::endl
//     << "inline std::vector<typename stan::return_type<T1, T2>::type>" << std::endl
//     << "operator()("
//     << "const T0& t, "
//     << "const std::vector<T1>& y, "
//     << "const std::vector<T2>& theta, "
//     << "const std::vector<double>& x_r, "
//     << "const std::vector<int>& x_i, "
//     << "std::ostream* msgs) const {" << std::endl;
    
//   size_t n = pmx_integrate_ode_group::CALLED_FUNCTORS.size();
//   for (size_t i = 0; i < n; ++i) {
//     std::string fs = pmx_integrate_ode_group::CALLED_FUNCTORS[i];
//     o << "if (id == " << i << ") {" << std::endl
//       << "return " << model_name << "_namespace::" << fs << "_functor__(t, y, theta, x_r, x_i, msgs);" << std::endl
//       << "}" << std::endl;
//   }

//   o << "}" << std::endl;
//   o << "};" << std::endl;
// }

// void generate_torsten_mpi(const std::string& model_name,
//                            std::ostream& o) {
//   o << "#undef PMX_ODE_GROUP_FUNCTOR_CALL" << std::endl
//     << "#define PMX_ODE_GROUP_FUNCTOR_CALL(CALLID, RESULT, __VA_ARGS__) \\" << std::endl << "{ \\" << std::endl;

//   size_t n = pmx_integrate_ode_group::CALLED_FUNCTORS.size();
//   for (size_t i = 0; i < n; ++i) {
//     std::string fs = pmx_integrate_ode_group::CALLED_FUNCTORS[i];
//     o << "if (CALLID == " << i << ") RESULT = " << model_name << "_namespace::" << fs << "_functor__(__VA_ARGS__);\\" << std::endl;
//   }

//   o << "}" << std::endl;
// }

}
}
#endif
