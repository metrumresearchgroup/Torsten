#ifndef STAN_LANG_AST_NODE_TORSTEN_PMX_INTEGRATE_ODE_GROUP_CONTROL_HPP
#define STAN_LANG_AST_NODE_TORSTEN_PMX_INTEGRATE_ODE_GROUP_CONTROL_HPP

#include <stan/lang/ast/node/expression.hpp>
#include <string>

namespace stan {
  namespace lang {

    struct expression;

    /**
     * Structure for integrate diff eq statement.
     */
    struct pmx_integrate_ode_group_control {
      /**
       * The name of the integrator.
       */
      std::string integration_function_name_;

      /**
       * Name of the ODE system.
       */
      std::string system_function_name_;

      /**
       * Initial state.
       */
      expression y0_;

      /**
       * Initial time.
       */
      expression t0_;

      /**
       * length of each part of ragged array @c ts
       */
      expression len_;

      /**
       * Solution times.
       */
      expression ts_;

      /**
       * Parameters.
       */
      expression theta_;  // params

      /**
       * Real-valued data.
       */
      expression x_;

      /**
       * Integer-valued data.
       */
      expression x_int_;

      /**
       * Relative tolerance (real).
       */
      expression rel_tol_;

      /**
       * Absolute tolerance (real).
       */
      expression abs_tol_;

      /**
       * Maximum number of steps (integer).
       */
      expression max_num_steps_;

      /**
       * Construct a default integrate ODE node.
       */
      pmx_integrate_ode_group_control();

      /**
       * Construct an integrate ODE node with the specified
       * components. 
       *
       * @param integration_function_name name of integrator
       * @param system_function_name name of ODE system
       * @param y0 initial value
       * @param t0 initial time
       * @param len length of each component of ragged array @c ts
       * @param ts solution times
       * @param theta parameters
       * @param x real-valued data
       * @param x_int integer-valued data
       */
      pmx_integrate_ode_group_control(const std::string& integration_function_name,
                                      const std::string& system_function_name,
                                      const expression& y0,
                                      const expression& t0,
                                      const expression& len,
                                      const expression& ts,
                                      const expression& theta,
                                      const expression& x,
                                      const expression& x_int,
                                      const expression& rel_tol,
                                      const expression& abs_tol,
                                      const expression& max_steps);
    };

  }
}
#endif
