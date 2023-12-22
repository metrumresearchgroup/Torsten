#ifndef STAN_LANG_AST_NODE_TORSTEN_PMX_INTEGRATE_ODE_GROUP_CONTROL_DEF_HPP
#define STAN_LANG_AST_NODE_TORSTEN_PMX_INTEGRATE_ODE_GROUP_CONTROL_DEF_HPP

#include <stan/lang/ast.hpp>
#include <string>

namespace stan {
  namespace lang {

    pmx_integrate_ode_group_control::pmx_integrate_ode_group_control() { }

    pmx_integrate_ode_group_control::pmx_integrate_ode_group_control(const std::string& integration_function_name, // NOLINT
                                                                     const std::string& system_function_name, // NOLINT
                                                                     const expression& y0, const expression& t0, // NOLINT
                                                                     const expression& len, // NOLINT
                                                                     const expression& ts, const expression& theta, // NOLINT
                                                                     const expression& x, const expression& x_int, // NOLINT
                                                                     const expression& rel_tol, // NOLINT
                                                                     const expression& abs_tol, // NOLINT
                                                                     const expression& max_steps) // NOLINT
    : integration_function_name_(integration_function_name),
      system_function_name_(system_function_name),
      y0_(y0), t0_(t0), len_(len), ts_(ts), theta_(theta), x_(x), x_int_(x_int),
      rel_tol_(rel_tol), abs_tol_(abs_tol), max_num_steps_(max_steps) {}
  }
}
#endif
