#ifndef STAN_MATH_TORSTEN_DSOLVE_CVODES_JAC_HPP
#define STAN_MATH_TORSTEN_DSOLVE_CVODES_JAC_HPP

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
    inline CVDlsJacFn cvodes_jac() {
      return [](realtype t, N_Vector y, N_Vector fy, SUNMatrix J, void* user_data, // NOLINT
                N_Vector tmp1, N_Vector tmp2, N_Vector tmp3) -> int {
        Ode* ode = static_cast<Ode*>(user_data);
        ode -> eval_jac(t, y, fy, J);
        return 0;
      };
    }

  }
}

#endif
