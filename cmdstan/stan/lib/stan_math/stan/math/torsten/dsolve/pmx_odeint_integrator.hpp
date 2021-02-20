#ifndef STAN_MATH_TORSTEN_DSOLVE_ODEINT_INTEGRATOR_HPP
#define STAN_MATH_TORSTEN_DSOLVE_ODEINT_INTEGRATOR_HPP

#include <stan/math/torsten/dsolve/pmx_ode_system.hpp>
#include <stan/math/torsten/mpi/precomputed_gradients.hpp>
#include <stan/math/torsten/dsolve/ode_observer.hpp>
#include <boost/numeric/odeint.hpp>
#include <type_traits>

namespace torsten {
namespace dsolve {

/**
 * @c boost::odeint ODE integrator.
 */
  template<typename scheme_t>
  struct PMXOdeintIntegrator {
    static constexpr int lmm_t = -1;

    const double rtol_;
    const double atol_;
    const int64_t max_num_steps_;

  public:
    /**
     * constructor
     * @param[in] rtol relative tolerance
     * @param[in] atol absolute tolerance
     * @param[in] max_num_steps max nb. of times steps
     */
    PMXOdeintIntegrator(const double rtol, const double atol,
                        const int64_t max_num_steps)
      : rtol_(rtol), atol_(atol), max_num_steps_(max_num_steps)
    {}

    template <typename Ode, typename Observer>
    inline void integrate(Ode& ode, Observer& observer) {
      std::vector<double> ts_vec(ode.ts_.size() + 1);
      ts_vec[0] = ode.t0_;
      for (size_t i = 0; i < ode.ts_.size(); ++i) {
        ts_vec[i + 1] = stan::math::value_of(ode.ts_[i]);
      }

      const double init_dt = 0.1;
      integrate_times(make_dense_output(atol_, rtol_, scheme_t()),
                      boost::ref(ode), ode.y0_fwd_system,
                      ts_vec.begin(), ts_vec.end(),
                      init_dt, boost::ref(observer),
                      boost::numeric::odeint::max_step_checker(max_num_steps_));
    }
  };

}
}
#endif
