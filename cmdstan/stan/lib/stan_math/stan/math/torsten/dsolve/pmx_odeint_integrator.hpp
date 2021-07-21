#ifndef STAN_MATH_TORSTEN_DSOLVE_ODEINT_INTEGRATOR_HPP
#define STAN_MATH_TORSTEN_DSOLVE_ODEINT_INTEGRATOR_HPP

#include <stan/math/torsten/dsolve/pmx_ode_system.hpp>
#include <stan/math/torsten/mpi/precomputed_gradients.hpp>
#include <stan/math/torsten/dsolve/ode_observer.hpp>
#include <boost/numeric/odeint.hpp>
#include <boost/numeric/odeint/external/eigen/eigen.hpp>
#include <type_traits>

namespace torsten {
namespace dsolve {

  using odeint_scheme_rk45 = boost::numeric::odeint::runge_kutta_dopri5<std::vector<double>, double, std::vector<double>, double>;
  using odeint_scheme_ckrk = boost::numeric::odeint::runge_kutta_cash_karp54<std::vector<double>, double, std::vector<double>, double>;

  // FIXME: vector_space_algebra for Eigen::VectorXd is not efficient
  // compared to std::vector
  // using odeint_scheme_rk45 = boost::numeric::odeint::runge_kutta_dopri5<Eigen::VectorXd, double, Eigen::VectorXd, double, boost::numeric::odeint::vector_space_algebra>;
  // using odeint_scheme_ckrk = boost::numeric::odeint::runge_kutta_cash_karp54<Eigen::VectorXd, double, Eigen::VectorXd, double, boost::numeric::odeint::vector_space_algebra>;

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
      // integrate_times(make_dense_output(atol_, rtol_, scheme_t()),
      integrate_times(make_stepper(scheme_t()),
                      boost::ref(ode), ode.y0_fwd_system,
                      ts_vec.begin(), ts_vec.end(),
                      init_dt, boost::ref(observer),
                      boost::numeric::odeint::max_step_checker(max_num_steps_));

      // integrate_times(make_controlled(absolute_tolerance, relative_tolerance,
      //                                 runge_kutta_cash_karp54<std::vector<double>, double,
      //                                 std::vector<double>, double>()),
      //                 std::ref(coupled_system), initial_coupled_state, std::begin(ts_vec),
      //                 std::end(ts_vec), step_size, filtered_observer,
      //                 max_step_checker(max_num_steps));
    }

    template<typename... Ts>
    auto make_stepper(boost::numeric::odeint::runge_kutta_dopri5<Ts...> const& stepper)
    {
      return make_dense_output(atol_, rtol_, stepper);
    }

    template<typename... Ts>
    auto make_stepper(boost::numeric::odeint::runge_kutta_cash_karp54<Ts...> const& stepper)
    {
      return make_controlled(atol_, rtol_, stepper);
    }
  };
}
}
#endif
