#ifndef STAN_MATH_TORSTEN_DSOLVE_ODEINT_INTEGRATOR_HPP
#define STAN_MATH_TORSTEN_DSOLVE_ODEINT_INTEGRATOR_HPP

#include <stan/math/torsten/dsolve/pmx_odeint_system.hpp>
#include <stan/math/torsten/mpi/precomputed_gradients.hpp>
#include <type_traits>

namespace torsten {
namespace dsolve {

/**
 * @c boost::odeint ODE integrator.
 */
  template<typename scheme_t>
  class PMXOdeintIntegrator {
    const double rtol_;
    const double atol_;
    const int64_t max_num_steps_;

    template<typename Ode, bool GenVar>
    struct SolObserver {
      const Ode& ode_;
      const size_t n;
      const size_t m;
      std::vector<std::vector<typename Ode::scalar_t>> y;
      int step_counter_;

      SolObserver(const Ode& ode) :
        ode_(ode), n(ode.N_), m(ode.M_),
        y(ode.ts_.size(), std::vector<typename Ode::scalar_t>(ode.N_, 0.0)),
        step_counter_(0)
      {}

      inline void operator()(const std::vector<double>& curr_result, double t) {
        if(t > ode_.t0_) {
          observer_impl(y[step_counter_], curr_result, ode_.y0_, ode_.theta_);
          step_counter_++;
        }
      }

    private:
      inline void observer_impl(std::vector<double>& y_res,
                         const std::vector<double>& y,
                         const std::vector<double>& y0,
                         const std::vector<double>& theta) const {
        std::copy(y.begin(), y.end(), y_res.begin());
      }

      inline void observer_impl(std::vector<stan::math::var>& y_res,
                         const std::vector<double>& y,
                         const std::vector<double>& y0,
                         const std::vector<stan::math::var>& theta) const {
        y_res = torsten::precomputed_gradients(y, theta);
      }

      inline void observer_impl(std::vector<stan::math::var>& y_res,
                         const std::vector<double>& y,
                         const std::vector<stan::math::var>& y0,
                         const std::vector<double>& theta) const {
        y_res = torsten::precomputed_gradients(y, y0);
      }

      inline void observer_impl(std::vector<stan::math::var>& y_res,
                         const std::vector<double>& y,
                         const std::vector<stan::math::var>& y0,
                         const std::vector<stan::math::var>& theta) const {
        std::vector<stan::math::var> vars = y0;
        vars.insert(vars.end(), theta.begin(), theta.end());
        y_res = torsten::precomputed_gradients(y, vars);
      }
    };

    template<typename Ode>
    struct SolObserver<Ode, false> {
      const Ode& ode_;
      Eigen::MatrixXd y;
      int step_counter_;

      SolObserver(const Ode& ode) :
        ode_(ode),
        y(Eigen::MatrixXd::Zero(ode_.size_, ode_.ts_.size())),
        step_counter_(0)
      {}

      inline void operator()(const std::vector<double>& curr_result, double t) {
        if(t > ode_.t0_) {
          y.col(step_counter_) = Eigen::VectorXd::Map(curr_result.data(), ode_.size_);
          step_counter_++;
        }
      }
    };

  public:
    /**
     * constructor
     * @param[in] rtol relative tolerance
     * @param[in] atol absolute tolerance
     * @param[in] max_num_steps max nb. of times steps
     */
    PMXOdeintIntegrator(const double rtol, const double atol,
                        const int64_t max_num_steps)
      : rtol_(rtol), atol_(atol), max_num_steps_(max_num_steps) {
      using stan::math::invalid_argument;
      if (rtol_ <= 0)
        invalid_argument("cvodes_integrator", "relative tolerance,", rtol_, "",
                         ", must be greater than 0");
      if (rtol_ > 1.0E-3)
        invalid_argument("cvodes_integrator", "relative tolerance,", rtol_, "",
                         ", must be less than 1.0E-3");
      if (atol_ <= 0)
        invalid_argument("cvodes_integrator", "absolute tolerance,", atol_, "",
                         ", must be greater than 0");
      if (max_num_steps_ <= 0)
        invalid_argument("cvodes_integrator", "max_num_steps,",
                         max_num_steps_, "",
                         ", must be greater than 0");
    }

    template <typename Ode, bool GenVar = true>
    auto integrate(Ode& ode) {
      std::vector<double> ts_vec(ode.ts_.size() + 1);
      ts_vec[0] = ode.t0_;
      std::copy(ode.ts_.begin(), ode.ts_.end(), ts_vec.begin() + 1);

      SolObserver<Ode, GenVar> observer(ode);

      const double init_dt = 0.1;
      integrate_times(make_dense_output(atol_, rtol_, scheme_t()),
                      boost::ref(ode), ode.y0_fwd_system,
                      ts_vec.begin(), ts_vec.end(),
                      init_dt, boost::ref(observer),
                      boost::numeric::odeint::max_step_checker(max_num_steps_));
      return observer.y;
    }
  };

}
}
#endif
