#ifndef TORSTEN_DSOLVE_CVODES_RHS_HPP
#define TORSTEN_DSOLVE_CVODES_RHS_HPP

#include <stan/math/rev/scal/meta/is_var.hpp>
#include <cvodes/cvodes.h>
#include <cvodes/cvodes_direct.h>
#include <sunmatrix/sunmatrix_dense.h>
#include <sunlinsol/sunlinsol_dense.h>
#include <nvector/nvector_serial.h>
#include <ostream>
#include <vector>
#include <algorithm>

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
