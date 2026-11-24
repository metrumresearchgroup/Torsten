#ifndef STAN_MATH_TORSTEN_DSOLVE_ARKODE_INTEGRATOR_HPP
#define STAN_MATH_TORSTEN_DSOLVE_ARKODE_INTEGRATOR_HPP

#include <stan/math/rev/fun/value_of.hpp>
#include <stan/math/prim/fun/value_of.hpp>
#include <stan/math/torsten/dsolve/sundials_check.hpp>
#include <stan/math/torsten/dsolve/cvodes_service.hpp>
// #include <arkode/arkode_butcher_erk.h>
#include <type_traits>

namespace torsten {
namespace dsolve {

/**
 * ARKODE ODE integrator.
 */
  template<int butcher_tab>
  struct PMXArkodeIntegrator {
    static constexpr int ARKODE_MAX_STEPS = 500;

    const double rtol_;
    const double atol_;
    const int64_t max_num_steps_;

    /**
     * constructor
     * @param[in] rtol relative tolerance
     * @param[in] atol absolute tolerance
     * @param[in] max_num_steps max nb. of times steps
     */
    PMXArkodeIntegrator(const double rtol, const double atol,
                        const int64_t max_num_steps = ARKODE_MAX_STEPS)
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
      PMXOdeService<Ode, 0, butcher_tab> serv(ode.N, ode.M, ode.ns, ode);
      N_Vector& y = serv.nv_y;
      const size_t n = ode.N;
      const size_t ns = ode.ns;

      // Initial condition is from nv_y, which has changed
      // from previous solution, we need to reset it. Likewise
      for (size_t i = 0; i < ode.system_size; ++i) {
        NV_Ith_S(y, i) = ode.y0_fwd_system[i];
      }

      void* mem = serv.mem;
      
      CHECK_SUNDIALS_CALL(ERKStepReInit(mem, Ode::arkode_combined_rhs, ode.t0_, y));
      CHECK_SUNDIALS_CALL(ERKStepSStolerances(mem, rtol_, atol_));
      CHECK_SUNDIALS_CALL(ERKStepSetMaxNumSteps(mem, max_num_steps_));
      CHECK_SUNDIALS_CALL(ERKStepSetUserData(mem, static_cast<void*>(&ode)));      

      CHECK_SUNDIALS_CALL(ERKStepSetTableNum(mem, butcher_tab));
      CHECK_SUNDIALS_CALL(ERKStepSetAdaptivityMethod(mem, 3, SUNTRUE, SUNFALSE, NULL));

      double t1 = ode.t0_;
      
      for (size_t i = 0; i < ode.ts_.size(); ++i) {
        CHECK_SUNDIALS_CALL(ERKStepEvolve(mem, stan::math::value_of(ode.ts_[i]), y, &t1, ARK_NORMAL));
        for (size_t j = 0; j < ode.system_size; ++j) {
          ode.y0_fwd_system[j] = NV_Ith_S(y, j);
        }
        observer(ode.y0_fwd_system, t1);
      }
    }

  };  // arkode integrator
}  // namespace dsolve
}  // namespace torsten

#endif
