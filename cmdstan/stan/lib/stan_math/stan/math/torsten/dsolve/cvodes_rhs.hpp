#ifndef STAN_MATH_TORSTEN_DSOLVE_CVODES_RHS_HPP
#define STAN_MATH_TORSTEN_DSOLVE_CVODES_RHS_HPP

#include <cvodes/cvodes.h>
#include <cvodes/cvodes_direct.h>
#include <nvector/nvector_serial.h>

namespace torsten {
  namespace dsolve {
    /**
     * return a closure for CVODES residual callback using a
     * non-capture lambda.
     *
     * @tparam Ode type of Ode
     * @return RHS function for Cvodes
     */
    template <typename Ode>
    inline CVRhsFn cvodes_rhs() {
      return [](double t, N_Vector y, N_Vector ydot, void* user_data) -> int {
        Ode* ode = static_cast<Ode*>(user_data);
        ode -> eval_rhs(t, y, ydot);
        return 0;
      };
    }

  }
}

#endif
