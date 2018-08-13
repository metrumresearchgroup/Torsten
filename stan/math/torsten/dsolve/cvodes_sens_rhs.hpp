#ifndef TORSTEN_DSOLVE_CVODES_SENS_RHS_HPP
#define TORSTEN_DSOLVE_CVODES_SENS_RHS_HPP

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
    static CVSensRhsFn cvodes_sens_rhs() {
    return [](int ns, double t, N_Vector y, N_Vector ydot,
              N_Vector* ys, N_Vector* ysdot, void* user_data,
              N_Vector temp1, N_Vector temp2) -> int {
      Ode* ode = static_cast<Ode*>(user_data);
      ode -> eval_sens_rhs(ns, t, y, ydot, ys, ysdot, temp1, temp2);
      return 0;
    };
  }

  }
}

#endif
