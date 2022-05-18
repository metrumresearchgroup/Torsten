#ifndef STAN_MATH_TORSTEN_ODE_MODEL_HPP
#define STAN_MATH_TORSTEN_ODE_MODEL_HPP

#include <stan/math/prim/fun/to_vector.hpp>
#include <stan/math/prim/fun/to_array_1d.hpp>
#include <stan/math/rev/fun/to_var.hpp>
#include <stan/math/rev/functor/algebra_solver_powell.hpp>
#include <stan/math/rev/functor/algebra_solver_newton.hpp>
#include <stan/math/rev/functor/algebra_solver_fp.hpp>
#include <stan/math/prim/err/check_less_or_equal.hpp>
#include <stan/math/torsten/pk_nvars.hpp>
#include <stan/math/torsten/ode_rhs_ostream_adatpor.hpp>
#include <stan/math/torsten/dsolve/pmx_algebra_solver_newton.hpp>
#include <stan/math/torsten/dsolve/pmx_ode_integrator.hpp>

namespace torsten {

  using Eigen::Matrix;
  using Eigen::Dynamic;
  using torsten::dsolve::PMXOdeIntegrator;

  /** 
   * Functor class that solves ODE & apply infusion dosing.
   */
  template<typename F>
  struct PMXOdeFunctorRateAdaptor {
    F const& f_;

    PMXOdeFunctorRateAdaptor(F const& f) : f_(f) {}

    /**
     * Evaluate ODE functor, and add rate.
     */
    template <typename T0, typename T1, typename T2, typename T3>
    Eigen::Matrix<typename stan::return_type<T0, T1, T2, T3>::type, -1, 1>
    operator()(const T0& t,
               const Eigen::Matrix<T1, -1, 1>& y,
               std::ostream* msgs,
               const std::vector<T2>& theta,
               const std::vector<T3>& rate,
               const std::vector<double>& x_r,
               const std::vector<int>& x_i) const {
      Eigen::Matrix<stan::return_type_t<T0, T1, T2, T3>, -1, 1> fy = f_(t, y, msgs, theta, x_r, x_i);
      for (auto i = 0; i < fy.size(); ++i) {
        fy(i) += rate[i];
      }
      return fy;
    }
  };

  template <typename F, typename integrator_type>
  struct PMXOdeFunctorSSAdaptor {
    F const& f_;
    integrator_type const& integrator_;
    
    PMXOdeFunctorSSAdaptor(F const& f, integrator_type const& integrator)
      : f_(f), integrator_(integrator)
    {}

    /**
     * When rate is RV, it's attached to parameters vector.
     * IN this case parameter @c y consists of {theta, rate}
     * This is used AFTER newton iteration to calculate <code>Jxy</code>.
     */
    template <typename T0, typename T1, typename T_amt, typename T_r, typename T_ii>
    Eigen::Matrix<typename stan::return_type_t<T0, T1, T_amt, T_r, T_ii>, -1, 1>
    operator()(const Eigen::Matrix<T0, -1, 1>& x,
               std::ostream* msgs,
               const std::vector<T1>& y,
               const T_amt& amt, const T_r& rate, const T_ii& ii, const int& cmt,
               const std::vector<double>& x_r,
               const std::vector<int>& x_i) const {
      using stan::math::to_vector;
      using stan::math::value_of;

      static const char* func("Steady State Event");
      using scalar_t = typename stan::return_type_t<T0, T1, T_amt, T_r, T_ii>;

      double t0 = 0;            // FIXME: ODE explicitly depdends on time

      const int ncmt = x.size();
      Eigen::Matrix<scalar_t, -1, 1> x0(ncmt);
      Eigen::Matrix<scalar_t, -1, 1> result(ncmt);

      for (auto i = 0; i < ncmt; ++i) {
        x0(i) = x(i);
      }

      if (rate == 0) {  // bolus dose
        x0[cmt - 1] += amt;
        auto pred = integrator_(f_, x0, t0, ii, y, x_r, x_i);
        for (int i = 0; i < ncmt; ++i) {
#ifdef TORSTEN_AS_FP
          result(i) = pred(i);
#else
          result(i) = x(i) - pred(i);
#endif
        }
      } else if (ii > 0) {  // multiple truncated infusions
        auto dt = amt / rate;
        // consider dosing internal > ii
        int n = int(std::floor(value_of(dt) / value_of(ii)) + 0.1);
        auto dt1 = dt - n * ii;
        std::vector<T_r> rate_vec(ncmt, 0.0);
        rate_vec[cmt - 1] = (n + 1) * rate;
        const PMXOdeFunctorRateAdaptor<F> f1(f_);
        x0 = integrator_(f1, x, t0, dt1, y, rate_vec, x_r, x_i);

        auto dt2 = (n + 1) * ii - dt;
        rate_vec[cmt - 1] = n * rate;
        const PMXOdeFunctorRateAdaptor<F> f2(f_);
        auto pred = integrator_(f2, x0, t0, dt2, y, rate_vec, x_r, x_i);
        for (int i = 0; i < ncmt; ++i) {
#ifdef TORSTEN_AS_FP
          result(i) = pred(i);
#else
          result(i) = x(i) - pred(i);
#endif
        }
      } else {  // constant infusion
        stan::math::check_less_or_equal(func, "AMT", amt, 0); 
        std::vector<T_r> rate_vec(ncmt, 0.0);
        rate_vec[cmt - 1] = rate;
#ifdef TORSTEN_AS_FP
        stan::math::throw_domain_error(func, "algebra_solver_fp used for ", 1, "constant infusion");
#else
        const PMXOdeFunctorRateAdaptor<F> f(f_);
        result = f(0, x, nullptr, y, rate_vec, x_r, x_i);
#endif
      }

      return result;
    }
  };

  /**
   * ODE-based PKPD models.
   *
   * @tparam T_par PK parameters type
   * @tparam F ODE functor
   */
  template<typename T_par, typename F0>
  class PKODEModel {
    using F = ode_rhs_ostream_adaptor<F0>;
    static const double dt_min; /**< min step to move ode solution */
    const std::vector<double> x_r_dummy; /**< dummy data to point to*/
    const std::vector<int> x_i_dummy; /**< dummy data to point to */
    const std::vector<T_par> &par_; /**< parameters */
    const std::vector<double>& x_r_; /**< real data */
    const std::vector<int>& x_i_; /**< integer data */
    const F f_;                /**< ODE functor */
    const int ncmt_;            /**< dim of ODE system */
  public:
    using par_type    = T_par;
    using f_type      = F;

    /**
     * Constructor
     *
     * @param par model parameters
     * @param f ODE functor
     * @param ncmt the ODE size.
     */
    PKODEModel(const std::vector<T_par> &par, int ncmt, const F& f) :
      x_r_dummy(), x_i_dummy(),
      par_(par), x_r_(x_r_dummy), x_i_(x_i_dummy),f_(f), ncmt_(ncmt)
    {}

    /**
     * Constructor
     *
     * @param par model parameters
     * @param f ODE functor
     * @param ncmt the ODE size.
     */
    PKODEModel(const std::vector<T_par> &par,
               const std::vector<double>& x_r,
               const std::vector<int>& x_i,
               int ncmt, const F& f) :
      par_(par), x_r_(x_r), x_i_(x_i),f_(f), ncmt_(ncmt)
    {}

    /**
     * Constructor
     *
     * @param par model parameters
     * @param f ODE functor
     * @param ncmt the ODE size.
     */
    PKODEModel(const std::vector<T_par> &par,
               const std::vector<double>& x_r,
               int ncmt, const F& f) :
      x_r_dummy(), x_i_dummy(),
      par_(par), x_r_(x_r), x_i_(x_i_dummy),f_(f), ncmt_(ncmt)
    {}

    /**
     * Constructor from any other model type as long as it
     * provides enough information to build an ODE model.
     *
     * @tparam T_model model type that can be used as ODE model.
     * @tparam Ts type parameters for @c T_model
     * @param m model that provides ODE information.
     */
    template<template<typename...> class T_model, typename... Ts>    
    PKODEModel(const T_model<Ts...>& m) :
      par_(m.par()),
      x_r_(m.x_r_), x_i_(m.x_i_),
      f_(m.f()),
      ncmt_(m.ncmt())
    {}
    
    template<typename T0, typename T1, typename T2, typename T3>
    static size_t nvars(const T0& t0,
                        const PKRec<T1>& y0,
                        const std::vector<T2> &rate,
                        const std::vector<T3> &par) {
      using stan::is_var;
      size_t res = 0;
      if (is_var<T0>::value) res += 1;
      if (is_var<T1>::value) res += y0.size();
      if (is_var<T2>::value) res += rate.size();
      if (is_var<T3>::value) res += par.size();
      return res;
    }

    /*
     * generate @c var vector that consisting of all the
     * parameters used passed to an ODE integrator. The output
     * will be arranged in order (y0, par, rate, dt)
     */
    template<typename T0, typename T1, typename T2, typename T3>
    static std::vector<stan::math::var> vars(const T0& t1,
                                             const PKRec<T1>& y0,
                                             const std::vector<T2> &rate,
                                             const std::vector<T3> &par) {
      using stan::math::to_var;      
      using stan::is_var;      
      using stan::math::var;
      std::vector<stan::math::var> res(nvars(t1, y0, rate, par));
      int n = nvars(t1, y0, rate, par);

      int i = 0;
      if (is_var<T1>::value) {
        for (int j = 0; j < y0.size(); ++j) {
          res[i] = y0(j);
          i++;
        }
      }
      if (is_var<T2>::value) {
        // FIXME: don't prepend par if it's not var
        if (is_var<T3>::value) {
          for (size_t j = 0; j < par.size(); ++j) {
            res[i] = par[j];
            i++;
          }
        }
        for (size_t j = 0; j < rate.size(); ++j) {
          res[i] = rate[j];
          i++;
        }
      } else if (is_var<T3>::value) {
        for (size_t j = 0; j < par.size(); ++j) {
          res[i] = par[j];
          i++;
        }
      }

      if (is_var<T0>::value) {
        res[i] = t1;
      }
      return res;
    }

    /**
     * @return model parameters
     */
    const std::vector<T_par>  & par() const { return par_; }
    /**
     * @return RHS functor
     */
    const F                   & f () const { return f_; }
    /**
     * @return ODE size
     */
    const int                 & ncmt () const { return ncmt_; }

  private:
    /**
     * When time is param, the time step needs to
     * incorporate information of initial time.
     */
    template<typename Tt0>
    std::vector<stan::math::var> time_step(const Tt0& t0, const stan::math::var& t1) const {
      return {stan::math::value_of(t0) + t1 - t0};
    }

    /**
     * When time is data, the time step is trivial
     */
    template<typename Tt0>
    std::vector<double> time_step(const Tt0& t0, const double t1) const {
      return {t1};
    }

    /**
     * For steady-state with data rate. 
     * Same as the non-steady version but with given @c init.
     */
    template<typename integrator_type>
    Eigen::Matrix<double, Eigen::Dynamic, 1>
    integrate(double t1,
              const std::vector<double> &rate,
              const Eigen::Matrix<double, 1, Eigen::Dynamic>& y0,
              const double& dt,
              const integrator_type& integrator) const {
      using stan::math::value_of;

      const double t0 = t1 - dt;
      std::vector<double> ts{t1};
      Eigen::Matrix<double, Eigen::Dynamic, 1> res;
      Eigen::Matrix<double, -1, 1> y0_(y0);
      if (ts[0] == t0) {
        res = y0;
      } else {
        const std::vector<double> pars{value_of(par_)};
        PMXOdeFunctorRateAdaptor<F> f(f_);
        res = integrator(f, y0_, t0, ts[0], pars, rate, x_r_, x_i_);
      }
      return res;
    }

  public:
    /**
     * solve ODE system. The different cases when @c rate is
     * @c var or data are handled by private methods @c integrate.
     *
     * @parm[in] t_next next time point when result is to be
     *           solved. The actual time point will be @c t0+dt
     * @parm[in] rtol relative tolerance for ODE solver
     * @parm[in] atol absolute tolerance for ODE solver
     * @parm[in] max_num_steps max number of steps between @c t0 and @c t0+dt.
     * @parm[in] msgs output stream.
     * @return an Eigen matrix with each row being the
     *         solution at certain time step. Hence the returned
     *         matrix is of dim (numer of time steps) x (siez of ODE system).
     */
    template<typename Tt0, typename Tt1, typename T, typename T1, typename integrator_type>
    void solve(Eigen::Matrix<T, -1, 1>& y,
               const Tt0& t0, const Tt1& t1,
               const std::vector<T1>& rate,
               const integrator_type& integrator) const {
      const double t0_d = stan::math::value_of(t0);
      std::vector<Tt1> ts(time_step(t0, t1));
      PMXOdeFunctorRateAdaptor<F> f_rate(f_);
      y = integrator(f_rate, y, t0_d, ts[0], par_, rate, x_r_, x_i_);
    }

    /** 
     * Torsten's integrators can return 
     * results in form of data directly,
     * thanks to @c pmx_cvodes/arkode/odeint_integrator's observer implmenetation.
     * 
     * @param yd output data
     * @param y initial condition
     * @param t0 starting time
     * @param t1 end tme
     * @param rate infusion rate
     * @param integrator ODE integrator
     */
    template<typename T0, typename T, typename T1, typename integrator_type>
    void solve_d(Eigen::VectorXd& yd,
                 const PKRec<T>& y,
                 const T0& t0, const T0& t1,
                 const std::vector<T1>& rate,
                 const integrator_type& integrator) const {
      using stan::math::var;
      using stan::math::value_of;
      using stan::math::to_var;

      const double t0_d = value_of(t0);
      std::vector<T0> ts(time_step(t0, t1));
      PMXOdeFunctorRateAdaptor<F> f_rate(f_);

      if (t1 - t0 > dt_min) {
        yd = integrator.solve_d(f_rate, y, t0_d, ts, par_, rate, x_r_, x_i_).col(0);
      }
    }

    /**
     * Solve steady state ODE-based PKPD model with @c amt
     * and @c rate both being data. Specifically for
     * constant rate we solve ODEs of the form
     *
     * x′(t)=F(x,t,θ)+ Rate,
     *
     * and solution x satisfies
     *
     * F(x,t,θ) + Rate = 0.
     *
     * For the periodic steady-state resulting from a
     * periodic input D into one or more cpts at a constant
     * interval τ, the root finder is combined with the ODE
     * solver to solve
     *
     * x(t0 + τ) − x(t0) = 0.
     *
     * Denote x(t0) by x0, it is equivalent to solve x0 in
     *
     * x(t0 + τ) = x0, where in the interval (t0, t0 + τ)
     * x(t) satisfies the ODE with x0 as init condition.
     *
     * @tparam It ODE numerical integrator ID.
     * @tparam T_ii dosing interval type
     * @param amt dosing amount
     * @param rate dosing rate
     * @param ii dosing interval
     * @param cmt compartment where the dosing occurs
     * @param rtol relative tolerance of ODE integrator.
     * @param atol absolute tolerance of ODE integrator.
     * @param max_num_steps max number of steps of ODE integrator.
     * @param msgs @c ostream output of ODE integrator.
     * @return col vector of the steady state ODE solution
     * Solve steady state ODE-based PKPD model with @c amt
     * being param and @c rate both being data.
     *
     */
    template<typename integrator_type, typename T_amt, typename T_r, typename T_ii>
    Eigen::Matrix<typename stan::return_type_t<T_amt, T_r, T_par, T_ii>, Eigen::Dynamic, 1> // NOLINT
    solve(double t0, const T_amt& amt, const T_r& rate, const T_ii& ii, const int& cmt,
          const integrator_type& integrator) const {
      using stan::math::value_of;
      using stan::math::algebra_solver_powell;
      using stan::math::algebra_solver_newton;
      using stan::math::algebra_solver_fp;

      typedef typename stan::return_type_t<T_amt, T_r, T_par, T_ii> scalar;

      double ii_dbl = value_of(ii);
      Eigen::Matrix<double, 1, -1> init_dbl(Eigen::Matrix<double, 1, -1>::Zero(ncmt_));
      std::vector<double> rate_vec(ncmt_, 0);
      // std::vector<int> x_i{cmt, ncmt_, int(par_.size()), int(integrator.max_num_step)};

      if (rate == 0) {                     // bolus dose
        init_dbl(cmt - 1) = value_of(amt); // bolus as initial condition
      } else {                             // infusion
        rate_vec[cmt - 1] = value_of(rate);
      }

      const double init_dt = (rate == 0.0 || ii > 0) ? ii_dbl : 24.0;
      using fss_func_t = PMXOdeFunctorSSAdaptor<F, integrator_type>;
      // fss_func_t fss(f_, par_, amt, rate, ii, cmt, ncmt_, integrator);
      fss_func_t fss(f_, integrator);
      try {
#ifdef TORSTEN_AS_POWELL
      return algebra_solver_powell(fss, integrate(t0, rate_vec, init_dbl, init_dt, integrator),
                                   fss.adapted_param(),
                                   x_r_, x_i_, 0,
                                   integrator.as_rtol, integrator.as_atol, integrator.as_max_num_step);
#elif defined(TORSTEN_AS_FP)
      std::vector<double> u_scale(ncmt_, 1.0);
      std::vector<double> f_scale(ncmt_, 1.0);
      return algebra_solver_fp(fss, integrate(t0, rate_vec, init_dbl, init_dt, integrator),
                               fss.adapted_param(),
                               x_r_, x_i_,
                               u_scale, f_scale, 0,
                               integrator.as_atol, integrator.as_max_num_step);
#else
      std::vector<double> scaling(ncmt_, 1.0);
      return pmx_algebra_solver_newton_tol(fss, integrate(t0, rate_vec, init_dbl, init_dt, integrator),
                                           scaling, scaling,
                                           integrator.as_rtol, integrator.as_atol, integrator.as_max_num_step,
                                           nullptr, par_, amt, rate, ii, cmt, x_r_, x_i_);
#endif
      } catch (const std::exception& e) {
        const char *text =
          "Failed to find steady state, due to "
          "either system's intrinsic lack of such a state "
          "or improper algebra solver controls. \nDetails: ";
        throw std::runtime_error(std::string(text) + e.what());
      }
    }
  };

  template<typename T_par, typename F>
  const double PKODEModel<T_par, F>::dt_min = 1.e-12;
}

#endif
