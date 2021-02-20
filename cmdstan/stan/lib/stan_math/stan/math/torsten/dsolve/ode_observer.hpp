#ifndef STAN_MATH_TORSTEN_DSOLVE_ODE_OBSERVER_HPP
#define STAN_MATH_TORSTEN_DSOLVE_ODE_OBSERVER_HPP

#include <stan/math/torsten/mpi/precomputed_gradients.hpp>
#include <stan/math/torsten/is_var.hpp>
#include <type_traits>

namespace torsten {
namespace dsolve {

  template<typename Ode>
  struct OdeObserver {
    const Ode& ode;
    const size_t n;
    const size_t m;
    const size_t ns;
    std::vector<std::vector<typename Ode::scalar_t>> y;
    int step_counter_;
    double curr_t_;

    OdeObserver(const Ode& ode0) :
      ode(ode0), n(ode.N), m(ode.M), ns(ode.ns),
      y(ode.ts_.size(), std::vector<typename Ode::scalar_t>(ode.N, 0.0)),
      step_counter_(0)
    {}

    inline void reset() {
      step_counter_ = 0;
      for (size_t i = 0; i < ode.ts_.size(); ++i) {
        std::fill(y[i].begin(), y[i].end(), 0);        
      }
    }

    /*
     * use observer to convert y value and gradient to var
     * results, if necessary.
     */
    inline void operator()(const std::vector<double>& curr_result, double t) {
      curr_t_ = t;
      if(curr_t_ > ode.t0_) {
        observer_impl(y[step_counter_], curr_result, ode);
        step_counter_++;
      }
    }

    inline void operator()(const N_Vector& curr_y, const N_Vector* curr_ys, double t) {
      curr_t_ = t;
      if(curr_t_ > ode.t0_) {
        observer_impl(y[step_counter_], curr_y, curr_ys, ode);
        step_counter_++;
      }
    }

  private:

    /*
     * All data, return data
     */
    template<typename F>
    inline void observer_impl(std::vector<double>& y_res,
                              const std::vector<double> & y,
                              const PMXOdeSystem<F, double, double, double>& system) const {
      std::copy(y.begin(), y.end(), y_res.begin());
    }

    template<typename F>
    inline void observer_impl(std::vector<double>& y_res,
                              const N_Vector& y,
                              const N_Vector* ys,
                              const PMXOdeSystem<F, double, double, double>& system) const {
      std::copy(NV_DATA_S(y), NV_DATA_S(y) + n, y_res.begin());
    }

    /*
     * When only @c ts is @c var, we don't solve
     * sensitivity ODE since the sensitivity is simply the RHS.
     */
    template<typename F>
    inline void observer_impl(std::vector<stan::math::var>& y_res,
                              const std::vector<double> & y,
                              const PMXOdeSystem<F, stan::math::var, double, double>& system) const {
      std::vector<double> g(n * (1 + system.ts_.size()), 0.0);
      std::copy(y.begin(), y.end(), g.begin());        
      std::vector<double> dydt(n);
      system.dbl_rhs_impl(curr_t_, y, dydt);
      std::copy(dydt.begin(), dydt.end(), g.begin() + n + step_counter_ * n);
      y_res = torsten::precomputed_gradients(g, system.ts_);
    }

    /*
     * When only @c ts is @c var, we don't solve
     * sensitivity ODE since the sensitivity is simply the RHS.
     */
    template<typename F>
    inline void observer_impl(std::vector<stan::math::var>& y_res,
                              const N_Vector & y,
                              const N_Vector* ys,
                              const PMXOdeSystem<F, stan::math::var, double, double>& system) const {
      std::vector<double> g(system.ts_.size(), 0.0);
      std::vector<double> dydt(n);
      system.dbl_rhs_impl(curr_t_, y, dydt);
      for (size_t j = 0; j < n; ++j) {
        g[step_counter_] = dydt[j];
        // FIXME: use ts[i] instead of ts
        y_res[j] = stan::math::precomputed_gradients(NV_Ith_S(y, j), system.ts_, g);
        g[step_counter_] = 0.0;
      }
    }

    /*
     * Only @c theta is @c var
     */
    template<typename F, typename T_init, typename T_par,
             std::enable_if_t<torsten::is_var<T_init, T_par>::value>* = nullptr>
    inline void observer_impl(std::vector<stan::math::var>& y_res,
                              const std::vector<double>& y,
                              const PMXOdeSystem<F, double, T_init, T_par>& system) const {
      y_res = torsten::precomputed_gradients(y, system.vars());
    }

    /*
     * Only @c theta is @c var
     */
    template<typename F, typename T_init, typename T_par,
             std::enable_if_t<torsten::is_var<T_init, T_par>::value>* = nullptr>
    inline void observer_impl(std::vector<stan::math::var>& y_res,
                              const N_Vector& y,
                              const N_Vector* ys,
                              const PMXOdeSystem<F, double, T_init, T_par>& system) const {
      y_res = torsten::precomputed_gradients(y, ys, system.vars());
    }

    /*
     * @c theta and/or &c y0 are @c var, together with @c ts.
     */
    template<typename F, typename T_init, typename T_par,
             std::enable_if_t<torsten::is_var<T_init, T_par>::value>* = nullptr>
    inline void observer_impl(std::vector<stan::math::var>& y_res,
                              const std::vector<double>& y,
                              const PMXOdeSystem<F, stan::math::var, T_init, T_par>& system) const {
      std::vector<double> g(n * (1 + ns + system.ts_.size()), 0.0);
      std::copy(y.begin(), y.end(), g.begin());        
      std::vector<double> dydt(n), y_dbl(y.begin(), y.begin() + n);
      system.dbl_rhs_impl(curr_t_, y_dbl, dydt);
      std::copy(dydt.begin(), dydt.end(), g.begin() + n + ns * n + step_counter_ * n);
      y_res = torsten::precomputed_gradients(g, system.vars());
    }

    // /*
    //  * @c theta and/or &c y0 are @c var, together with @c ts.
    //  */
    // template<typename F, typename T_init, typename T_par,
    //          std::enable_if_t<torsten::is_var<T_init, T_par>::value>* = nullptr>
    // inline void observer_impl(std::vector<stan::math::var>& y_res,
    //                           const std::vector<double>& y,
    //                           const PMXOdeSystem<F, stan::math::var, T_init, T_par>& system) const {
    //   std::vector<double> g(n * (1 + ns + system.ts_.size()), 0.0);
    //   std::copy(y.begin(), y.end(), g.begin());        
    //   std::vector<double> dydt(n), y_dbl(y.begin(), y.begin() + n);
    //   system.dbl_rhs_impl(curr_t_, y, dydt);
    //   std::copy(dydt.begin(), dydt.end(), g.begin() + n + ns * n + i * n);
    //   y_res = torsten::precomputed_gradients(g, ode_.vars());
    // }

    /*
     * @c theta and/or &c y0 are @c var, together with @c ts.
     */
    template<typename F, typename T_init, typename T_par,
             std::enable_if_t<torsten::is_var<T_init, T_par>::value>* = nullptr>
    inline void observer_impl(std::vector<stan::math::var>& y_res,
                              const N_Vector& y,
                              const N_Vector* ys,
                              const PMXOdeSystem<F, stan::math::var, T_init, T_par>& system) const {
      std::vector<double> g(ns + system.ts_.size(), 0.0);
      std::vector<double> dydt(n);
      system.dbl_rhs_impl(curr_t_, y, dydt);
      for (size_t k = 0; k < n; ++k) {
        for (size_t j = 0; j < ns; ++j) {
          g[j] = NV_Ith_S(ys[j], k);
        }
        g[step_counter_ + ns] = dydt[k];
        y_res[k] = precomputed_gradients(NV_Ith_S(y, k), system.vars(), g);
        g[step_counter_ + ns] = 0.0;
      }
    }
  };

  template<typename Ode>
  struct OdeDataObserver {
    const Ode& ode;
    const size_t n;
    const size_t m;
    const size_t ns;
    const size_t nt;
    Eigen::MatrixXd y;
    int step_counter_;
    double curr_t_;

    OdeDataObserver(const Ode& ode0) :
      ode(ode0), n(ode.N), m(ode.M), ns(ode.ns), nt(ode.ts_.size()),
      y(Eigen::MatrixXd::Zero(ode.system_size + n*(Ode::is_var_ts ? nt : 0), nt)),
      step_counter_(0)
    {}

    inline void reset(const Ode& ode_) {
      step_counter_ = 0;
      y.setZero();
    }

    inline void operator()(const std::vector<double>& curr_result, double t) {
      curr_t_ = t;
      if(curr_t_ > ode.t0_) {
        observer_impl(y, curr_result, ode);
        step_counter_++;
      }
    }

    inline void operator()(const N_Vector& curr_y, const N_Vector* curr_ys, double t) {
      curr_t_ = t;
      if(curr_t_ > ode.t0_) {
        observer_impl(y, curr_y, curr_ys, ode);
        step_counter_++;
      }
    }

  private:
    /*
     * @c ts is data
     */
    template<typename F, typename T_init, typename T_par>
    inline void observer_impl(Eigen::MatrixXd& y_res,
                              const std::vector<double>& y,
                              const PMXOdeSystem<F, double, T_init, T_par>& system) const {
      y_res.col(step_counter_) = Eigen::VectorXd::Map(y.data(), system.system_size);
    }

    /*
     * @c ts is @c var
     */
    template<typename F, typename T_init, typename T_par>
    inline void observer_impl(Eigen::MatrixXd& y_res,
                              const std::vector<double>& y,
                              const PMXOdeSystem<F, stan::math::var, T_init, T_par>& system) const {
      for (size_t j = 0; j < system.system_size; ++j) y_res(j, step_counter_) = y[j];
      std::vector<double> dydt(n);
      std::vector<double> y_tmp(y.begin(), y.begin() + n);
      system.dbl_rhs_impl(curr_t_, y_tmp, dydt);
      for (size_t j = 0; j < n; ++j) {
        y_res(system.system_size + step_counter_ * n + j, step_counter_) = dydt[j];        
      }
    }

    /*
     * @c ts is @c var
     */
    template<typename F, typename Tt, typename T_init, typename T_par>
    inline void observer_impl(Eigen::MatrixXd& y_res,
                              const N_Vector& y,
                              const N_Vector* ys,
                              const PMXOdeSystem<F, Tt, T_init, T_par>& system) const {
      y_res.block(0, step_counter_, n, 1) = Eigen::VectorXd::Map(NV_DATA_S(y), n);

      if (system.use_fwd_sens) {
        for (size_t j = 0; j < ns; ++j) {
          y_res.block(n + j * n, step_counter_, n, 1) = Eigen::VectorXd::Map(NV_DATA_S(ys[j]), n);
        }
      }

      if (system.is_var_ts) {
        std::vector<double> dydt(n);
        system.dbl_rhs_impl(curr_t_, y, dydt);
        y_res.block(n + (step_counter_ + ns) * n, step_counter_, n, 1) =
          Eigen::VectorXd::Map(dydt.data(), n);
      }
    }
  };

}
}
#endif
