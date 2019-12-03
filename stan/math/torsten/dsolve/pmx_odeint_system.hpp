#ifndef STAN_MATH_TORSTEN_DSOLVE_PMX_ODEINT_SYSTEM_HPP
#define STAN_MATH_TORSTEN_DSOLVE_PMX_ODEINT_SYSTEM_HPP

#include <stan/math/torsten/dsolve/cvodes_service.hpp>
#include <stan/math/torsten/dsolve/ode_forms.hpp>
#include <stan/math/torsten/return_type.hpp>
#include <stan/math/torsten/dsolve/pmx_ode_vars.hpp>
#include <stan/math/prim/arr/meta/get.hpp>
#include <stan/math/prim/arr/meta/length.hpp>
#include <stan/math/prim/arr/fun/value_of.hpp>
#include <stan/math/prim/scal/err/check_size_match.hpp>
#include <stan/math/rev/scal/fun/value_of_rec.hpp>
#include <stan/math/rev/core.hpp>
#include <ostream>
#include <stdexcept>
#include <vector>

namespace torsten {
namespace dsolve {

  /**
   * Boost Odeint ODE system that contains informtion on residual
   * equation functor, sensitivity residual equation functor,
   * as well as initial conditions. This is a base type that
   * is intended to contain common values used by forward
   * sensitivity system.
   *
   * @tparam F type of functor for ODE residual
   * @tparam Tt scalar type of time steps
   * @tparam T_init scalar type of initial unknown values
   * @tparam T_par scalar type of parameters
   */
  template <typename F, typename Tt, typename T_init, typename T_par>
  struct PMXOdeintSystem {
    using Ode = PMXOdeintSystem<F, Tt, T_init, T_par>;
    using scalar_t = typename torsten::return_t<Tt, T_init, T_par>::type;
    static constexpr bool is_var_ts  = stan::is_var<Tt>::value;
    static constexpr bool is_var_y0  = stan::is_var<T_init>::value;
    static constexpr bool is_var_par = stan::is_var<T_par>::value;

    const F& f_;
    const double t0_;
    const std::vector<Tt>& ts_;
    const std::vector<T_init>& y0_;
    const std::vector<T_par>& theta_;
    const std::vector<double> theta_dbl_;
    const std::vector<double>& x_r_;
    const std::vector<int>& x_i_;
    const size_t N_;
    const size_t M_;
    const size_t ns;
    const size_t size_;
    std::ostream* msgs_;
    std::vector<double>& y0_fwd_system;
  private:
    int step_counter_;  

  public:
    template<typename ode_t>
    PMXOdeintSystem(dsolve::PMXOdeService<ode_t>& serv,
                    const F& f,
                    double t0,
                    const std::vector<Tt>& ts,
                    const std::vector<T_init>& y0,
                    const std::vector<T_par>& theta,
                    const std::vector<double>& x_r,
                    const std::vector<int>& x_i,
                    std::ostream* msgs)
      : f_(f),
        t0_(t0),
        ts_(ts),
        y0_(y0),
        theta_(theta),
        theta_dbl_(stan::math::value_of(theta)),
        x_r_(x_r),
        x_i_(x_i),
        N_(y0.size()),
        M_(theta.size()),      
        ns(serv.ns),
        size_(serv.size),
        msgs_(msgs),
        y0_fwd_system(serv.y),
        step_counter_(0)
    {
      // initial state
      std::fill(y0_fwd_system.begin(), y0_fwd_system.end(), 0.0);
      if (is_var_y0)  {
      std::transform(y0.begin(), y0.end(), y0_fwd_system.begin(),
                     [](const T_init& v){ return stan::math::value_of(v); });        
      for (size_t i = 0; i < N_; i++) y0_fwd_system[N_ + i * N_ + i] = 1.0;
      } else {
      std::transform(y0.begin(), y0.end(), y0_fwd_system.begin(),
                     [](const T_init& v){ return stan::math::value_of(v); });
      }
    }

    inline const size_t fwd_system_size() const { return size_; }

    inline const std::vector<Tt> & ts() const { return ts_; }

    inline const std::vector<T_init>& y0() const { return y0_; }

    inline const std::vector<T_par>& theta() const { return theta_; }

    /*
     * retrieving a vector of vars that will be used as parameters
     */
    inline auto vars() const {
      return pmx_ode_vars(y0_, theta_, ts_);
    }

    /*
     * Evaluate RHS of the ODE
     * @param y current dependent value, arranged as {y, dy_dp1, dy_dp2...}
     * @param dy_dt ODE RHS to be filled.
     * @param t current indepedent value
     */
    void operator()(const std::vector<double>& y, std::vector<double>& dy_dt,
                    double t) const {
      stan::math::check_size_match("PMXOdeintSystem", "y", y.size(), "dy_dt", dy_dt.size());
      rhs_impl(y, dy_dt, t);
    }

    /*
     * evaluate RHS with data only inputs.
     */
    void dbl_rhs_impl(const std::vector<double>& y, std::vector<double>& dy_dt, double t) const
    {
      stan::math::check_size_match("PMXOdeintSystem", "y", y.size(), "dy_dt", dy_dt.size());
      dy_dt = f_(t, y, theta_dbl_, x_r_, x_i_, msgs_);
      return;      
    }

    /*
     * evalute RHS of the entire system, possibly including
     * the forward sensitivity equation components in @c y and @c dy_dt.
     */
    void rhs_impl(const std::vector<double>& y, std::vector<double>& dy_dt, double t) const
    {
      using std::vector;
      using stan::math::var;

      if (!(is_var_y0 || is_var_par)) {
        dy_dt = f_(t, y, theta_dbl_, x_r_, x_i_, msgs_);
        return;
      }

      std::fill(dy_dt.begin(), dy_dt.end(), 0.0);
      try {
        stan::math::start_nested();

        std::vector<var> yv(y.begin(), y.begin() + N_);
        std::vector<var> theta_v(theta_dbl_.begin(), theta_dbl_.end());
        std::vector<var> fyv(is_var_par ?
                             f_(t, yv, theta_v, x_r_, x_i_, msgs_) :
                             f_(t, yv, theta_, x_r_, x_i_, msgs_));

        stan::math::check_size_match("PMXOdeintSystem", "dz_dt", fyv.size(), "states", N_);

        for (size_t i = 0; i < N_; ++i) {
          stan::math::set_zero_all_adjoints_nested();
          dy_dt[i] = fyv[i].val();
          fyv[i].grad();

          // df/dy*s_i term, for i = 1...ns
          for (size_t j = 0; j < ns; ++j) {
            for (size_t k = 0; k < N_; ++k) {
              dy_dt.at(N_ + N_ * j + i) += y[N_ + N_ * j + k] * yv[k].adj();
            }
          }

          // df/dp_i term, for i = n...n+m-1
          if (is_var_par) {
            for (size_t j = 0; j < M_; ++j) {
              dy_dt.at(N_ + N_ * (ns - M_ + j) + i) += theta_v[j].adj();
            }
          }
        }
      } catch (const std::exception& e) {
        stan::math::recover_memory_nested();
        throw;
      }
      stan::math::recover_memory_nested();
    }
  };

}  // namespace dsolve
}  // namespace torsten
#endif
