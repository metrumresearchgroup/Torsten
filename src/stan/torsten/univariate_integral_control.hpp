#ifndef STAN_LANG_TORSTEN_AST_NODE_UNIVARIATE_INTEGRAL_CONTROL_HPP
#define STAN_LANG_TORSTEN_AST_NODE_UNIVARIATE_INTEGRAL_CONTROL_HPP

#include <stan/lang/ast/node/expression.hpp>
#include <string>

namespace stan {
  namespace lang {

    struct expression;

    /**
     * Structure for a Torsten univariate_integral function with control
     * parameters for the integrator.
     */
    struct univariate_integral_control {
      /* The name of the function (* tells us which integrator is being
       * called).
       */
      std::string integration_function_name_;
      /**
       * Name of the functor that's be integrated
       */
      std::string system_function_name_;

      /**
       * left limit
       */
      expression t0_;

      /**
       * right limit
       */
      expression t1_;

      /**
       * parameters
       */
      expression theta_;

      /**
       * Real-valued data (array of real).
       */
      expression x_r_;

      /**
       * Integer-valued data (array of int).
       */
      expression x_i_;

      /**
       * default constructor.
       */
      univariate_integral_control();

      /**
       * Construct an univarate integral
       *
       * @param f functor for that's being integrated
       * @param t0 left end of interval
       * @param t1 right end of interval
       */
      univariate_integral_control(
                                  const std::string& integration_function_name,
                                  const std::string& system_function_name,
                                  const expression& t0,
                                  const expression& t1,
                                  const expression& theta,
                                  const expression& x_r,
                                  const expression& x_i);
    };

  }
}
#endif

