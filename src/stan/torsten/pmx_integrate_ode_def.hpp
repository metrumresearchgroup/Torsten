#ifndef STAN_LANG_AST_NODE_TORSTEN_INTEGRATE_ODE_DEF_HPP
#define STAN_LANG_AST_NODE_TORSTEN_INTEGRATE_ODE_DEF_HPP

#include <stan/lang/ast.hpp>
#include <string>

namespace stan {
  namespace lang {

    pmx_integrate_ode::pmx_integrate_ode() { }

    pmx_integrate_ode::pmx_integrate_ode(const std::string& integration_function_name, // NOLINT
                                         const std::string& system_function_name, // NOLINT
                                         const expression& y0, const expression& t0, // NOLINT
                                         const expression& ts, const expression& theta, // NOLINT
                                         const expression& x, const expression& x_int) // NOLINT
      : integration_function_name_(integration_function_name),
        system_function_name_(system_function_name),
        y0_(y0), t0_(t0), ts_(ts), theta_(theta), x_(x), x_int_(x_int) {  }

  }
}
#endif
