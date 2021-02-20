#ifndef STAN_MATH_TORSTEN_DSOLVE_PMX_ODE_SYSTEM_HPP
#define STAN_MATH_TORSTEN_DSOLVE_PMX_ODE_SYSTEM_HPP

#include <stan/math/rev/core/recover_memory.hpp>
#include <stan/math/prim/meta/return_type.hpp>
#include <stan/math/torsten/dsolve/ode_check.hpp>
#include <stan/math/torsten/dsolve/pmx_ode_vars.hpp>
#include <stan/math/prim/fun/get.hpp>
#include <stan/math/prim/fun/value_of.hpp>
#include <stan/math/prim/err/check_size_match.hpp>
#include <stan/math/rev/fun/value_of_rec.hpp>
#include <stan/math/rev/core.hpp>
#include <cvodes/cvodes.h>
#include <nvector/nvector_serial.h>
#include <sunmatrix/sunmatrix_dense.h>
#include <ostream>
#include <stdexcept>
#include <vector>

namespace torsten {
namespace dsolve {

  /**
   * ODE system that contains informtion on residual
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
  struct  PMXOdeSystem {
    using Ode = PMXOdeSystem<F, Tt, T_init, T_par>;
    using scalar_t = typename stan::return_type_t<Tt, T_init, T_par>;
    static constexpr bool is_var_ts  = stan::is_var<Tt>::value;
    static constexpr bool is_var_y0  = stan::is_var<T_init>::value;
    static constexpr bool is_var_par = stan::is_var<T_par>::value;
    static constexpr bool use_fwd_sens = is_var_y0 || is_var_par;

    const F& f_;
    const double t0_;
    const std::vector<Tt>& ts_;
    const std::vector<T_init>& y0_;
    const std::vector<T_par>& theta_;
    const std::vector<double> theta_dbl_;
    const std::vector<double>& x_r_;
    const std::vector<int>& x_i_;
    const size_t N;
    const size_t M;
    const size_t ns;
    const size_t system_size;
    std::ostream* msgs_;
    std::vector<double> y0_fwd_system;

  public:
    PMXOdeSystem(const F& f,
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
        N(y0.size()),
        M(theta.size()),      
        ns((is_var_y0 ? N : 0) + (is_var_par ? M : 0)),
        system_size(N + N * ns),
        msgs_(msgs),
        y0_fwd_system(system_size, 0.0)
    {
      const char* caller = "PMX ODE System";
      torsten::dsolve::ode_check(y0_, t0_, ts_, theta_, x_r_, x_i_, caller);

      // initial state
      std::transform(y0.begin(), y0.end(), y0_fwd_system.begin(),
                     [](const T_init& v){ return stan::math::value_of(v); });
      if (is_var_y0)  {
        for (size_t i = 0; i < N; ++i) {
          y0_fwd_system[N + i * N + i] = 1.0;        
        }
      }
    }

    /*
     * retrieving a vector of vars that will be used as parameters
     */
    inline auto vars() const {
      return pmx_ode_vars(y0_, theta_, ts_);
    }

    /*
     * Evaluate RHS of the ODE(the combined system)
     * @param y current dependent value, arranged as {y, dy_dp1, dy_dp2...}
     * @param dy_dt ODE RHS to be filled.
     * @param t current indepedent value
     */
    inline void operator()(const std::vector<double>& y, std::vector<double>& dy_dt,
                           double t) const {
      stan::math::check_size_match("PMXOdeSystem", "y", y.size(), "dy_dt", dy_dt.size());
      rhs_impl(y, dy_dt, t);
    }

    /*
     * evaluate RHS with data only inputs.
     */
    inline void dbl_rhs_impl(double t, const std::vector<double>& y, std::vector<double>& dy_dt) const
    {
      stan::math::check_size_match("PMXOdeSystem", "y", y.size(), "dy_dt", dy_dt.size());
      dy_dt = f_(t, y, theta_dbl_, x_r_, x_i_, msgs_);
      return;      
    }

    /*
     * evaluate RHS with data only inputs.
     */
    inline void dbl_rhs_impl(double t, const N_Vector& nv_y, std::vector<double>& dy_dt) const
    {
      stan::math::check_size_match("PMXOdeSystem", "y", NV_LENGTH_S(nv_y), "dy_dt", dy_dt.size());
      std::vector<double> y(NV_DATA_S(nv_y), NV_DATA_S(nv_y) + N);
      dy_dt = f_(t, y, theta_dbl_, x_r_, x_i_, msgs_);
      return;      
    }

    /*
     * evaluate RHS with data only inputs for N_Vector data
     */    
    inline void operator()(double t, N_Vector& nv_y, N_Vector& ydot) const {
      stan::math::check_size_match("PMXOdeSystem", "y", NV_LENGTH_S(nv_y), "dy_dt", NV_LENGTH_S(ydot));
      std::vector<double> y(NV_DATA_S(nv_y), NV_DATA_S(nv_y) + N);
      std::vector<double> dydt(f_(t, y, theta_dbl_, x_r_, x_i_, msgs_));
      for (size_t i = 0; i < N; ++i) {
        NV_Ith_S(ydot, i) = dydt[i];
      }
    }

    static int cvodes_rhs(double t, N_Vector y, N_Vector ydot, void* user_data) {
      Ode* ode = static_cast<Ode*>(user_data);
      (*ode)(t, y, ydot);
      return 0;
    }

    /*
     * evalute RHS of the entire system, possibly including
     * the forward sensitivity equation components in @c y and @c dy_dt.
     */
    void rhs_impl(const std::vector<double>& y,
                  std::vector<double>& dy_dt, double t) const {
      using std::vector;
      using stan::math::var;

      if (!(is_var_y0 || is_var_par)) {
        dy_dt = f_(t, y, theta_dbl_, x_r_, x_i_, msgs_);
        return;
      }

      std::fill(dy_dt.begin(), dy_dt.end(), 0.0);
      stan::math::nested_rev_autodiff nested;

      std::vector<var> yv(y.begin(), y.begin() + N);
      std::vector<var> theta_v(theta_dbl_.begin(), theta_dbl_.end());
      std::vector<var> fyv(is_var_par ?
                           f_(t, yv, theta_v, x_r_, x_i_, msgs_) :
                           f_(t, yv, theta_, x_r_, x_i_, msgs_));

      stan::math::check_size_match("PMXOdeSystem", "dz_dt", fyv.size(), "states", N);

      for (size_t i = 0; i < N; ++i) {
        if (i > 0) {
          nested.set_zero_all_adjoints();            
        }
        dy_dt[i] = fyv[i].val();
        fyv[i].grad();

        // df/dy*s_i term, for i = 1...ns
        for (size_t j = 0; j < ns; ++j) {
          for (size_t k = 0; k < N; ++k) {
            dy_dt.at(N + N * j + i) += y[N + N * j + k] * yv[k].adj();
          }
        }

        // df/dp_i term, for i = n...n+m-1
        if (is_var_par) {
          for (size_t j = 0; j < M; ++j) {
            dy_dt.at(N + N * (ns - M + j) + i) += theta_v[j].adj();
          }
        }
      }
    }

    /**
     * Calculate sensitivity rhs using CVODES vectors. The
     * internal workspace is allocated by @c PMXOdeService.
     */
    inline void operator()(int ns, double t, N_Vector nv_y, N_Vector ydot,
                           N_Vector* ys, N_Vector* ysdot,
                           N_Vector temp1, N_Vector temp2) {
      // for (int i = 0; i < N; ++i) y[i] = NV_Ith_S(nv_y, i);

      // initialize ysdot
      for (int i = 0; i < ns; ++i) N_VConst(0.0, ysdot[i]);

      stan::math::nested_rev_autodiff nested;

      std::vector<stan::math::var> yv_work(NV_DATA_S(nv_y), NV_DATA_S(nv_y) + N);
      std::vector<stan::math::var> theta_work(theta_dbl_.begin(), theta_dbl_.end());
      std::vector<stan::math::var> dy_dt(is_var_par ?
                                         f_(t, yv_work, theta_work, x_r_, x_i_, msgs_) :
                                         f_(t, yv_work, theta_dbl_, x_r_, x_i_, msgs_));

      stan::math::check_size_match("PMXOdeSystem", "dy_dt", dy_dt.size(), "states", N);

      for (int j = 0; j < N; ++j) {
        if (j > 0) {
          nested.set_zero_all_adjoints();
        }
        dy_dt[j].grad();

        // df/dy*s_i term, for i = 1...ns
        for (int i = 0; i < ns; ++i) {
          auto ysp = N_VGetArrayPointer(ys[i]);
          auto nvp = N_VGetArrayPointer(ysdot[i]);
          for (int k = 0; k < N; ++k) {
            nvp[j] += yv_work[k].adj() * ysp[k];              
          }
        }

        // df/dp_i term, for i = n...n+m-1
        if (is_var_par) {
          for (int i = 0; i < M; ++i) {
            auto nvp = N_VGetArrayPointer(ysdot[ns - M + i]);
            nvp[j] += theta_work[i].adj();
          }
        }
      }
    }

    static int cvodes_sens_rhs(int ns, double t, N_Vector y, N_Vector ydot,
                               N_Vector* ys, N_Vector* ysdot, void* user_data,
                               N_Vector temp1, N_Vector temp2) {
      if (use_fwd_sens) {
        Ode* ode = static_cast<Ode*>(user_data);
        (*ode)(ns, t, y, ydot, ys, ysdot, temp1, temp2);            
      }
      return 0;
    }

    /**
     * return a closure for CVODES residual callback using a
     * non-capture lambda.
     *
     * @tparam Ode type of Ode
     * @return RHS function for Cvodes
     */
    inline CVLsJacFn cvodes_jac() {
      return [](realtype t, N_Vector y, N_Vector fy, SUNMatrix J, void* user_data, // NOLINT
                N_Vector tmp1, N_Vector tmp2, N_Vector tmp3) -> int {
        Ode* ode = static_cast<Ode*>(user_data);
        ode -> jac(t, y, fy, J);
        return 0;
      };
    }

    /**
     * evaluate Jacobian matrix using current state, store
     * the result in @c SUNMatrix J.
     *
     * @param t current time
     * @param y current y
     * @param fy current f(y)
     * @param J Jacobian matrix J(i,j) = df_i/dy_j
     */
    inline void jac(double t, N_Vector& nv_y, N_Vector& fy, SUNMatrix& J) {
      stan::math::nested_rev_autodiff nested;

      std::vector<stan::math::var> yv_work(NV_DATA_S(nv_y), NV_DATA_S(nv_y) + N);
      std::vector<stan::math::var> fyv_work(f_(t, yv_work, theta_dbl_, x_r_, x_i_, msgs_));

      for (int i = 0; i < N; ++i) {
        nested.set_zero_all_adjoints();
        fyv_work[i].grad();
        for (int j = 0; j < N; ++j) {
          SM_ELEMENT_D(J, i, j) = yv_work[j].adj();
        }
      }
    }
  };


}  // namespace dsolve
}  // namespace torsten
#endif
