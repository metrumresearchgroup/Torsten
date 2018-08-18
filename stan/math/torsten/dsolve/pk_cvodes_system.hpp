#ifndef STAN_MATH_TORSTEN_DSOLVE_CVODES_SYSTEM_HPP
#define STAN_MATH_TORSTEN_DSOLVE_CVODES_SYSTEM_HPP

#include <stan/math/prim/arr/fun/value_of.hpp>
#include <stan/math/prim/arr/fun/dot_self.hpp>
#include <stan/math/prim/scal/err/check_greater_or_equal.hpp>
#include <stan/math/prim/scal/err/check_less_or_equal.hpp>
#include <stan/math/prim/scal/err/check_finite.hpp>
#include <stan/math/prim/scal/err/check_less.hpp>
#include <stan/math/prim/scal/err/system_error.hpp>
#include <stan/math/prim/arr/err/check_nonzero_size.hpp>
#include <stan/math/prim/arr/err/check_ordered.hpp>
#include <stan/math/torsten/dsolve/sundials_check.hpp>
#include <stan/math/torsten/dsolve/cvodes_rhs.hpp>
#include <stan/math/torsten/dsolve/cvodes_service.hpp>
#include <stan/math/rev/scal/meta/is_var.hpp>
#include <stan/math/prim/scal/meta/return_type.hpp>
#include <stan/math/prim/mat/fun/typedefs.hpp>
#include <stan/math/rev/mat/fun/typedefs.hpp>
#include <cvodes/cvodes.h>
#include <nvector/nvector_serial.h>
#include <ostream>
#include <vector>

namespace torsten {
  namespace dsolve {

    /**
     * CVODES ODE system that contains informtion on residual
     * equation functor, sensitivity residual equation functor,
     * as well as initial conditions. This is a base type that
     * is intended to contain common values used by forward
     * sensitivity system.
     *
     * @tparam F type of functor for ODE residual
     * @tparam Ty0 scalar type of initial unknown values
     * @tparam Tpar scalar type of parameters
     */
    template <typename F, typename Tts,
              typename Ty0, typename Tpar, int Lmm>
    class pk_cvodes_system {
    public:
      using Ode = pk_cvodes_system<F, Tts, Ty0, Tpar, Lmm>;

    protected:
      cvodes_service<Ode>& serv_;
      const F& f_;
      const double t0_;
      const std::vector<Tts>& ts_;
      const std::vector<Ty0>& y0_;
      const std::vector<Tpar>& theta_;
      const std::vector<double> y0_dbl_;
      const std::vector<double> theta_dbl_;
      const std::vector<double>& x_r_;
      const std::vector<int>& x_i_;
      const size_t N_;
      const size_t M_;
      const size_t ns_;  // nb. of sensi params
      N_Vector& nv_y_;
      std::vector<double>& y_vec_;
      std::vector<double>& fval_;
      void* mem_;
      SUNMatrix& A_;
      SUNLinearSolver& LS_;
      std::ostream* msgs_;

    public:
      static constexpr bool is_var_y0 = stan::is_var<Ty0>::value;
      static constexpr bool is_var_par = stan::is_var<Tpar>::value;
      static constexpr bool need_fwd_sens = is_var_y0 || is_var_par;
      static constexpr int lmm_type = Lmm;

      // when ts is param, we don't have to do fwd
      // sensitivity by solving extra ODEs, because in this
      // case the sensitivity regarding ts is just the RHS
      // of ODE. We can append this type of sensitivity
      // results after CVODES solutions.
      static constexpr bool is_var_ts = stan::is_var<Tts>::value;

      using scalar_type = typename stan::return_type<Tts, Ty0, Tpar>::type;
      using return_type = std::vector<std::vector<scalar_type> >;

      /**
       * Construct CVODES ODE system from initial condition and parameters
       *
       * @param[in] f ODE residual functor
       * @param[in] y0 initial condiiton
       * @param[in] theta parameters of the base ODE.
       * @param[in] x_r continuous data vector for the ODE.
       * @param[in] x_i integer data vector for the ODE.
       * @param[in] msgs stream to which messages are printed.
       */
      pk_cvodes_system(cvodes_service<Ode>& serv,
                       const F& f,
                       double t0,
                       const std::vector<Tts>& ts,
                       const std::vector<Ty0>& y0,
                       const std::vector<Tpar>& theta,
                       const std::vector<double>& x_r,
                       const std::vector<int>& x_i,
                       std::ostream* msgs)
        : serv_(serv),
          f_(f),
          t0_(t0),
          ts_(ts),
          y0_(y0),
          theta_(theta),
          y0_dbl_(stan::math::value_of(y0)),
          theta_dbl_(stan::math::value_of(theta)),
          x_r_(x_r),
          x_i_(x_i),
          N_(y0.size()),
          M_(theta.size()),
          ns_((is_var_y0 ? N_ : 0) + (is_var_par ? M_ : 0)),
          nv_y_(serv_.nv_y),
          y_vec_(serv_.y),
          fval_(serv_.fval),
          mem_(serv_.mem),
          A_(serv_.A),
          LS_(serv.LS),
          msgs_(msgs) {
        using stan::math::system_error;

        if (nv_y_ == NULL)
          throw std::runtime_error("N_VMake_Serial failed to allocate memory");

        if (mem_ == NULL)
          throw std::runtime_error("CVodeCreate failed to allocate memory");

        auto t0_data = stan::math::value_of(t0);
        auto ts_data = stan::math::value_of(ts);
        
        static const char* caller = "pk_cvodes_system";
        stan::math::check_finite(caller, "initial time", t0_data);
        stan::math::check_finite(caller, "times", ts_data);
        stan::math::check_ordered(caller, "times", ts_data);
        stan::math::check_nonzero_size(caller, "times", ts_data);
        stan::math::check_less(caller, "initial time", t0_data, ts_data.front());
        stan::math::check_finite(caller, "initial state", y0);
        stan::math::check_finite(caller, "parameter vector", theta);
        stan::math::check_finite(caller, "continuous data", x_r);
        stan::math::check_nonzero_size(caller, "initial state", y0);

        int err = 1;
        if (N_ != serv_.N)
          system_error(caller, "N_", err, "inconsistent allocated memory");
        if (N_ != size_t(N_VGetLength_Serial(nv_y_)))
          system_error(caller, "nv_y", err, "inconsistent allocated memory");
        if (M_ != serv_.M)
          system_error(caller, "M_", err, "inconsistent allocated memory");
        if (ns_ != serv_.ns)
          system_error(caller, "ns_", err, "inconsistent allocated memory");

        // initial condition
        for (size_t i = 0; i < N_; ++i) {
          NV_Ith_S(nv_y_, i) = stan::math::value_of(y0_[i]);
        }
      }

      /**
       * destructor is empty as all CVODES resources are
       * handled by @c cvodes_service
       */
      ~pk_cvodes_system() {
      }

      /**
       * return reference to initial time
       *
       * @return reference to initial time
       */
      const double& t0() { return t0_; }

      /**
       * return reference to time steps
       *
       * @return reference to time steps
       */
      const std::vector<Tts> & ts() const { return ts_; }

      /**
       * return reference to current N_Vector of unknown variable
       *
       * @return reference to current N_Vector of unknown variable
       */
      N_Vector& nv_y() { return nv_y_; }

      /**
       * evaluate RHS function using current state, store
       * the result in internal @c fval_
       */
      void eval_rhs(double t, N_Vector& y) {
        for (size_t i = 0; i < N_; ++i) y_vec_[i] = NV_Ith_S(y, i);
        fval_ = f_(t, y_vec_, theta_dbl_, x_r_, x_i_, msgs_);
      }

      /**
       * evaluate RHS function using current state, store
       * the result in @c N_Vector.
       */
      void eval_rhs(double t, N_Vector& y, N_Vector& ydot) {
        for (size_t i = 0; i < N_; ++i) y_vec_[i] = NV_Ith_S(y, i);
        fval_ = f_(t, y_vec_, theta_dbl_, x_r_, x_i_, msgs_);
        for (size_t i = 0; i < N_; ++i) NV_Ith_S(ydot, i) = fval_[i];
      }

      /**
       * return reference to initial condition
       *
       * @return reference to initial condition
       */
      const std::vector<Ty0>& y0() const { return y0_; }

      /**
       * return reference to parameter
       *
       * @return reference to parameter
       */
      const std::vector<Tpar>& theta() const { return theta_; }

      /**
       * return a vector of vars for that contains the initial
       * condition and parameters in case they are vars. The
       * sensitivity with respect to this vector will be
       * calculated by CVODES.
       *
       * @return vector of vars
       */
      std::vector<scalar_type> vars() const {
        std::vector<scalar_type> res;
        if (is_var_y0) {
          res.insert(res.end(), y0().begin(), y0().end());
        }
        if (is_var_par) {
          res.insert(res.end(), theta().begin(), theta().end());
        }
        if (is_var_ts) {
          res.insert(res.end(), ts().begin(), ts().end());
        }
        return res;
      }

      /**
       * return current @c y_vec(). We also use it for workspace.
       */
      std::vector<double>& y_vec() { return y_vec_; }

      /**
       * return current RHS evaluation.
       */
      const std::vector<double>& fval() { return fval_; }

      /**
       * return number of unknown variables
       */
      const size_t n() { return N_; }

      /**
       * return number of sensitivity parameters
       */
      const size_t ns() { return ns_; }

      /**
       * return size of ODE system for primary and sensitivity unknowns
       */
      const size_t n_sys() { return N_ * (ns_ + 1); }

      /**
       * return theta size
       */
      const size_t n_par() { return theta_.size(); }

      /**
       * return CVODES memory handle
       */
      void* mem() { return mem_; }

      /**
       * return reference to ODE functor
       */
      const F& f() { return f_; }
      /**
       * return a closure for CVODES residual callback using a non-capture lambda
       */
      static CVRhsFn rhs() {
        return cvodes_rhs<pk_cvodes_system<F, Tts, Ty0, Tpar, Lmm>>();
      }
    };

    // // TODO(yizhang): adjoint system construction

  }  // namespace dsolve
}  // namespace torsten

#endif
