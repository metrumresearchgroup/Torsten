#ifndef STAN_MATH_TORSTEN_DSOLVE_CVODES_SENS_RHS_HPP
#define STAN_MATH_TORSTEN_DSOLVE_CVODES_SENS_RHS_HPP

#include <cvodes/cvodes.h>
#include <cvodes/cvodes_direct.h>
#include <nvector/nvector_serial.h>

namespace torsten {
  namespace dsolve {

    template <typename Ode, bool needs_sens>
    struct cvodes_sens_rhs_impl {
      static CVSensRhsFn f() {
        return [](int ns, double t, N_Vector y, N_Vector ydot,
                  N_Vector* ys, N_Vector* ysdot, void* user_data,
                  N_Vector temp1, N_Vector temp2) -> int {
          Ode* ode = static_cast<Ode*>(user_data);
          ode -> eval_sens_rhs(ns, t, y, ydot, ys, ysdot, temp1, temp2);
          return 0;
        };
      }
    };

    template <typename Ode>
    struct cvodes_sens_rhs_impl<Ode, false> {
      static CVSensRhsFn f() {
        return [](int ns, double t, N_Vector y, N_Vector ydot,
                  N_Vector* ys, N_Vector* ysdot, void* user_data,
                  N_Vector temp1, N_Vector temp2) -> int {return 0;};
      }
    };

    /**
     * return a closure for CVODES sensitivity RHS callback using a
     * non-capture lambda.
     *
     * @tparam Ode type of Ode
     * @return RHS function for Cvodes
     */
    template <typename Ode>
    static CVSensRhsFn cvodes_sens_rhs() {
      return cvodes_sens_rhs_impl<Ode, Ode::need_fwd_sens>::f();
    }
  }
}

#endif
