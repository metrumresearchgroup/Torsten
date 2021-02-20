#ifndef STAN_MATH_TORSTEN_DSOLVE_CVODES_INTEGRATOR_HPP
#define STAN_MATH_TORSTEN_DSOLVE_CVODES_INTEGRATOR_HPP

#include <stan/math/rev/fun/value_of.hpp>
#include <stan/math/prim/fun/value_of.hpp>
#include <stan/math/torsten/dsolve/sundials_check.hpp>
#include <stan/math/torsten/dsolve/cvodes_service.hpp>
#include <type_traits>

namespace torsten {
namespace dsolve {

/**
 * CVODES ODE integrator.
 */
  template<int lmm_type, int ism_type>
  struct PMXCvodesIntegrator {
    static constexpr int lmm_t = lmm_type;
    static constexpr int CVODES_MAX_STEPS = 500;

    const double rtol_;
    const double atol_;
    const int64_t max_num_steps_;

    /**
     * constructor
     * @param[in] rtol relative tolerance
     * @param[in] atol absolute tolerance
     * @param[in] max_num_steps max nb. of times steps
     */
    PMXCvodesIntegrator(const double rtol, const double atol,
                        const int64_t max_num_steps = CVODES_MAX_STEPS)
      : rtol_(rtol), atol_(atol), max_num_steps_(max_num_steps)
    {}
      
    /**
     * Return the solutions for the specified ODE
     * given the specified initial state,
     * initial times, times of desired solution, and parameters and
     * data, writing error and warning messages to the specified
     * stream contained in the ODE system.
     *
     * @tparam ODE type of ODE system
     * @param[in] ode ODE system
     * increasing order, all greater than the initial time.
     * @return a vector of states, each state being a vector of the
     * same size as the state variable, corresponding to a time in ts.
     */
    template <typename Ode, typename Observer>
    inline void integrate(Ode& ode, Observer& observer) {
      PMXOdeService<Ode, lmm_type> serv(ode.N, ode.M, ode.ns, ode);

      void* mem = serv.mem;
      N_Vector& y = serv.nv_y;
      N_Vector* ys = serv.nv_ys;
      const size_t n = ode.N;
      const size_t ns = ode.ns;

      // Initial condition is from nv_y, which has changed
      // from previous solution, we need to reset it. Likewise
      // we also reset ys.
      for (size_t i = 0; i < n; ++i) {
        NV_Ith_S(y, i) = stan::math::value_of(ode.y0_[i]);
      }
      for (size_t is = 0; is < ns; ++is) {
        N_VConst(RCONST(0.0), ys[is]);
      }

      CHECK_SUNDIALS_CALL(CVodeReInit(mem, ode.t0_, y));
      CHECK_SUNDIALS_CALL(CVodeSStolerances(mem, rtol_, atol_));
      CHECK_SUNDIALS_CALL(CVodeSetMaxNumSteps(mem, max_num_steps_));
      CHECK_SUNDIALS_CALL(CVodeSetMaxErrTestFails(mem, 20));
      CHECK_SUNDIALS_CALL(CVodeSetMaxConvFails(mem, 20));
#ifndef TORSTEN_CVS_JAC_DQ
      CHECK_SUNDIALS_CALL(CVodeSetJacFn(mem, ode.cvodes_jac()));
#endif

      /** if y0 is parameter, the first n sensitivity vector
       * are regarding y0, thus they form a unit matrix.
       **/
      if (Ode::use_fwd_sens) {
        if (Ode::is_var_y0) {
          for (size_t i = 0; i < n; ++i) NV_Ith_S(ys[i], i) = 1.0;          
        }
        CHECK_SUNDIALS_CALL(CVodeSensReInit(mem, ism_type, ys));
        CHECK_SUNDIALS_CALL(CVodeSensEEtolerances(mem));
      }

      double t1 = ode.t0_;
      
      for (size_t i = 0; i < ode.ts_.size(); ++i) {
        CHECK_SUNDIALS_CALL(CVode(mem, stan::math::value_of(ode.ts_[i]), y, &t1, CV_NORMAL));
        if (ode.use_fwd_sens) {
          CHECK_SUNDIALS_CALL(CVodeGetSens(mem, &t1, ys));
        }
        observer(y, ys, t1);
      }
    }
  };  // cvodes integrator
}  // namespace dsolve
}  // namespace torsten

#endif
