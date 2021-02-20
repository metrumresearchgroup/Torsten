#ifndef STAN_MATH_TORSTEN_DSOLVE_ODEINT_INTEGRATOR_HPP
#define STAN_MATH_TORSTEN_DSOLVE_ODEINT_INTEGRATOR_HPP

#include <stan/math/torsten/dsolve/pmx_odeint_system.hpp>
#include <stan/math/torsten/mpi/precomputed_gradients.hpp>
#include <boost/numeric/odeint.hpp>
#include <type_traits>

namespace torsten {
namespace dsolve {

/**
 * @c boost::odeint ODE integrator.
 */
  template<typename scheme_t>
  struct PMXOdeintIntegrator {
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

      /*
       * use observer to convert y value and gradient to var
       * results, if necessary.
       */
      inline void operator()(const std::vector<double>& curr_result, double t) {
        if(t > ode_.t0_) {
          observer_impl(y[step_counter_], curr_result, ode_.ts_, ode_.y0_, ode_.theta_, step_counter_);
          step_counter_++;
        }
      }

    private:
      /*
       * All data, return data
       */
      inline void observer_impl(std::vector<double>& y_res,
                                const std::vector<double>& y,
                                const std::vector<double>& ts,
                                const std::vector<double>& y0,
                                const std::vector<double>& theta,
                                int i) const {
        std::copy(y.begin(), y.end(), y_res.begin());
      }

      /*
       * When only @c ts is @c var, we don't solve
       * sensitivity ODE since the sensitivity is simply the RHS.
       */
      inline void observer_impl(std::vector<stan::math::var>& y_res,
                                const std::vector<double>& y,
                                const std::vector<stan::math::var>& ts,
                                const std::vector<double>& y0,
                                const std::vector<double>& theta,
                                int i) const {
        int n = y0.size();
        std::vector<double> g(n * (1 + ts.size()), 0.0);
        std::copy(y.begin(), y.end(), g.begin());        
        std::vector<double> dy_dt(n);
        ode_.dbl_rhs_impl(y, dy_dt, ts[i].val());
        std::copy(dy_dt.begin(), dy_dt.end(), g.begin() + n + i * n);
        y_res = torsten::precomputed_gradients(g, ts);
      }

      /*
       * Only @c theta is @c var
       */
      inline void observer_impl(std::vector<stan::math::var>& y_res,
                                const std::vector<double>& y,
                                const std::vector<double>& ts,
                                const std::vector<double>& y0,
                                const std::vector<stan::math::var>& theta,
                                int i) const {
        y_res = torsten::precomputed_gradients(y, theta);
      }

      /*
       * only @c y0 is @c var
       */
      inline void observer_impl(std::vector<stan::math::var>& y_res,
                                const std::vector<double>& y,
                                const std::vector<double>& ts,
                                const std::vector<stan::math::var>& y0,
                                const std::vector<double>& theta,
                                int i) const {
        y_res = torsten::precomputed_gradients(y, y0);
      }

      /*
       * @c y0 and @c theta are @c var
       */
      inline void observer_impl(std::vector<stan::math::var>& y_res,
                                const std::vector<double>& y,
                                const std::vector<double>& ts,
                                const std::vector<stan::math::var>& y0,
                                const std::vector<stan::math::var>& theta,
                                int i) const {
        y_res = torsten::precomputed_gradients(y, ode_.vars());
      }

      /*
       * @c theta and/or &c y0 are @c var, together with @c ts.
       */
      template<typename T1, typename T2>
      inline void observer_impl(std::vector<stan::math::var>& y_res,
                                const std::vector<double>& y,
                                const std::vector<stan::math::var>& ts,
                                const std::vector<T1>& y0,
                                const std::vector<T2>& theta,
                                int i) const {
        int ns = ode_.ns;
        int n = y0.size();
        std::vector<double> g(n * (1 + ns + ts.size()), 0.0);
        std::copy(y.begin(), y.end(), g.begin());        
        std::vector<double> dy_dt(n), y_dbl(y.begin(), y.begin() + n);
        ode_.dbl_rhs_impl(y_dbl, dy_dt, ts[i].val());
        std::copy(dy_dt.begin(), dy_dt.end(), g.begin() + n + ns * n + i * n);
        y_res = torsten::precomputed_gradients(g, ode_.vars());
      }
    };

    template<typename Ode>
    struct SolObserver<Ode, false> {
      const Ode& ode_;
      Eigen::MatrixXd y;
      int step_counter_;

      SolObserver(const Ode& ode) :
        ode_(ode),
        y(Eigen::MatrixXd::Zero(ode_.size_ + ode_.N_*(Ode::is_var_ts ? ode_.ts_.size() : 0), ode_.ts_.size())),
        step_counter_(0)
      {}

      inline void operator()(const std::vector<double>& curr_result, double t) {
        if(t > ode_.t0_) {
          observer_impl(y, curr_result, ode_.ts_, step_counter_);
          step_counter_++;
        }
      }

    private:
      /*
       * @@c ts is data
       */
      inline void observer_impl(Eigen::MatrixXd& y_res,
                                const std::vector<double>& y,
                                const std::vector<double>& ts,
                                int i) const {
        using Eigen::VectorXd;
        y_res.col(i) = VectorXd::Map(y.data(), ode_.size_);
      }

      /*
       * @@c ts is @c var
       */
      inline void observer_impl(Eigen::MatrixXd& y_res,
                                const std::vector<double>& y,
                                const std::vector<stan::math::var>& ts,
                                int i) const {
        using Eigen::VectorXd;
        for (size_t j = 0; j < ode_.size_; ++j) y_res(j, i) = y[j];
        int n = ode_.N_;
        std::vector<double> dy_dt(n);
        std::vector<double> y_tmp(y.begin(), y.begin() + n);
        ode_.dbl_rhs_impl(y_tmp, dy_dt, ts[i].val());
        for (size_t j = 0; j < n; ++j) y_res(ode_.size_ + i * n + j, i) = dy_dt[j];
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
      for (size_t i = 0; i < ode.ts_.size(); ++i) {
        ts_vec[i + 1] = stan::math::value_of(ode.ts_[i]);
      }

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
