#ifndef STAN_LANG_TORSTEN_AST_NODE_UNIVARIATE_INTEGRAL_CONTROL_DEF_HPP
#define STAN_LANG_TORSTEN_AST_NODE_UNIVARIATE_INTEGRAL_CONTROL_DEF_HPP

#include <stan/lang/ast.hpp>
#include <string>

namespace stan {
  namespace lang {

    univariate_integral_control::univariate_integral_control() { }

    univariate_integral_control::univariate_integral_control(
                                const std::string& integration_function_name,
                                const std::string& system_function_name,
                                const expression& t0,
                                const expression& t1,
                                const expression& theta,
                                const expression& x_r,
                                const expression& x_i)
      : integration_function_name_(integration_function_name),
        system_function_name_(system_function_name),
        t0_(t0),
        t1_(t1),
        theta_(theta),
        x_r_(x_r),
        x_i_(x_i)
    { }
  }
}
#endif
