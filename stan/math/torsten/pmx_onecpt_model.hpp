#ifndef STAN_MATH_TORSTEN_ONECPT_MODEL_HPP
#define STAN_MATH_TORSTEN_ONECPT_MODEL_HPP

#include <stan/math/torsten/PKModel/functors/check_mti.hpp>
#include <stan/math/torsten/PKModel/Pred/unpromote.hpp>
#include <stan/math/torsten/PKModel/Pred/PolyExp.hpp>
#include <stan/math/torsten/model_solve_d.hpp>
#include <stan/math/torsten/pmx_ode_integrator.hpp>
#include <stan/math/torsten/dsolve/pk_vars.hpp>
#include <stan/math/torsten/pk_nvars.hpp>
#include <stan/math/prim/scal/err/check_positive_finite.hpp>
#include <stan/math/prim/scal/err/check_finite.hpp>

namespace refactor {

  using boost::math::tools::promote_args;

  /**
   * standard one compartment PK ODE functor
   */
  struct PMXOneCptODE {
  /**
   * standard one compartment PK ODE RHS function
   * @tparam T0 t type
   * @tparam T1 initial condition type
   * @tparam T2 parameter type
   * @tparam T3 real data/rate type
   * @param t type
   * @param x initial condition type
   * @param parms parameters
   * @param rate dosing rate
   * @param dummy dummy
   */
    template <typename T0, typename T1, typename T2, typename T3>
    inline
    std::vector<typename stan::return_type<T0, T1, T2, T3>::type>
    operator()(const T0& t,
               const std::vector<T1>& x,
               const std::vector<T2>& parms,
               const std::vector<T3>& rate,
               const std::vector<int>& dummy,
               std::ostream* pstream__) const {
      typedef typename stan::return_type<T0, T1, T2, T3>::type scalar;

      scalar CL = parms.at(0), V1 = parms.at(1), ka = parms.at(2), k10 = CL / V1;
      std::vector<scalar> y(2, 0);

      y.at(0) = -ka * x.at(0);
      y.at(1) = ka * x.at(0) - k10 * x.at(1);

      return y;
    }
  };

  /**
   * One-compartment PK model. The static memebers provide
   * universal information, i.e. nb. of compartments,
   * nb. of parameters, and the RHS functor.
   *
   * @tparam T_time t type
   * @tparam T_init initial condition type
   * @tparam T_rate dosing rate type
   * @tparam T_par PK parameters type
   */
  template<typename T_time, typename T_init, typename T_rate, typename T_par>
  class PMXOneCptModel {
    const T_time &t0_;
    const refactor::PKRec<T_init>& y0_;
    const std::vector<T_rate> &rate_;
    const T_par &CL_;
    const T_par &V2_;
    const T_par &ka_;
    const T_par k10_;
    const std::vector<T_par> alpha_;
    const std::vector<T_par> par_;

  public:
    static constexpr int Ncmt = 2;
    static constexpr int Npar = 3;
    static constexpr PMXOneCptODE f_ = PMXOneCptODE();

    using scalar_type = typename stan::return_type<T_time, T_init, T_rate, T_par>::type;
    using init_type   = T_init;
    using time_type   = T_time;
    using par_type    = T_par;
    using rate_type   = T_rate;

  /**
   * One-compartment PK model constructor
   *
   * @param t0 initial time
   * @param y0 initial condition
   * @param rate dosing rate
   * @param par model parameters
   * @param CL clearance
   * @param V2 central cpt vol
   * @param ka absorption
   */
    PMXOneCptModel(const T_time& t0,
                  const refactor::PKRec <T_init>& y0,
                  const std::vector<T_rate> &rate,
                  const T_par& CL,
                  const T_par& V2,
                  const T_par& ka) :
      t0_(t0),
      y0_(y0),
      rate_(rate),
      CL_(CL),
      V2_(V2),
      ka_(ka),
      k10_(CL / V2),
      alpha_{k10_, ka_},
      par_{CL_, V2_, ka_}
    {
      using stan::math::check_positive_finite;
      using stan::math::check_finite;
      const char* fun = "PMXOneCptModel";
      check_positive_finite(fun, "CL", CL_);
      check_positive_finite(fun, "V2", V2_);
      check_positive_finite(fun, "ka", ka_);
    }

  /**
   * One-compartment PK model constructor
   * FIXME need to remove parameter as this is for linode only.
   *
   * @tparam T_mp parameters class
   * @tparam Ts parameter types
   * @param t0 initial time
   * @param y0 initial condition
   * @param rate dosing rate
   * @param par model parameters
   * @param parameter ModelParameter type
   */
    template<template<typename...> class T_mp, typename... Ts>
    PMXOneCptModel(const T_time& t0,
                  const Eigen::Matrix<T_init, 1, Eigen::Dynamic>& y0,
                  const std::vector<T_rate> &rate,
                  const std::vector<T_par> & par,
                  const T_mp<Ts...> &parameter) :
      PMXOneCptModel(t0, y0, rate, par.at(0), par.at(1), par.at(2))
    {}

  /**
   * One-compartment PK model constructor
   *
   * @param t0 initial time
   * @param y0 initial condition
   * @param rate dosing rate
   * @param par model parameters
   */
    PMXOneCptModel(const T_time& t0,
                  const Eigen::Matrix<T_init, 1, Eigen::Dynamic>& y0,
                  const std::vector<T_rate> &rate,
                  const std::vector<T_par> & par) :
      PMXOneCptModel(t0, y0, rate, par.at(0), par.at(1), par.at(2))
    {}

    /*
     * calculate number of @c vars for transient dosing.
     */
    static int nvars(int ncmt, int npar) {
      using stan::is_var;
      int n = 0;
      if (is_var<T_time>::value) n++; // t0
      if (is_var<T_init>::value) n += Ncmt; // y0 is fixed for onecpt model
      if (is_var<T_rate>::value) n += Ncmt; // rate is fixed for onecpt model
      if (is_var<T_par>::value) n += Npar; // par is fixed for onecpt model
      return n;
    }

    /*
     * calculate number of @c vars for steady-state dosing.
     */
    template<typename T_a, typename T_r, typename T_ii>
    static int nvars(int npar) {
      using stan::is_var;
      int n = 0;
      if (is_var<T_a>::value) n++; // amt
      if (is_var<T_r>::value) n++; // rate
      if (is_var<T_ii>::value) n++; // ii
      if (is_var<T_par>::value) n += Npar; // par is fixed for onecpt model
      return n;
    }

    /*
     * return the number @c var that will be the parameters
     * of the trasient dosing event's solution
     */
    template<typename T0>
    int nvars(const T0& t0) {
      return torsten::pk_nvars(t0, y0_, rate_, par_);
    }

    /*
     * return the number @c var that will be the parameters
     * of the stead-state dosing event's solution
     */
    template<typename T_a, typename T_r, typename T_ii>
    int nvars(const T_a& a, const T_r& r, const T_ii& ii) {
      return torsten::pk_nvars(a, r, ii, par_);
    }

    /*
     * return @c vars that will be solution
     */
    template<typename T0>
    std::vector<stan::math::var> vars(const T0 t1) {
      return torsten::dsolve::pk_vars(t1, y0_, rate_, par_);
    }

    /*
     * return @c vars that will be steady-state
     * solution. For SS solution @c rate_ or @ y0_ will not
     * be in the solution.
     */
    template<typename T_a, typename T_r, typename T_ii>
    std::vector<stan::math::var> vars(const T_a& a, const T_r& r, const T_ii& ii) {
      return torsten::dsolve::pk_vars(a, r, ii, par_);
    }

  /**
   * One-compartment PK model get methods
   */
    const T_time              & t0()      const { return t0_;    }
    const PKRec<T_init>       & y0()      const { return y0_;    }
    const std::vector<T_rate> & rate()    const { return rate_;  }
    const std::vector<T_par>  & alpha()   const { return alpha_; }
    const std::vector<T_par>  & par()     const { return par_;   }
    const PMXOneCptODE         & f()       const { return f_;     }
    const int                 & ncmt ()   const { return Ncmt;   }

    /**
     * Solve one-cpt model using analytical solution.
     *
     * @tparam T_time time type
     * @tparam T_model ODE model type
     */
    Eigen::Matrix<scalar_type, Eigen::Dynamic, 1>
    solve(const T_time& t_next) const {
      using Eigen::Matrix;
      using Eigen::Dynamic;

      T_time dt = t_next - t0_;

      std::vector<scalar_type> a(Ncmt, 0);
      Matrix<scalar_type, 1, Dynamic> pred = PKRec<scalar_type>::Zero(Ncmt);

      if ((y0_[0] != 0) || (rate_[0] != 0)) {
        pred(0, 0) = y0_[0] * exp(-ka_ * dt) + rate_[0] * (1 - exp(-ka_ * dt)) / ka_;
        a[0] = ka_ / (ka_ - alpha_[0]);
        a[1] = -a[0];
        pred(0, 1) += torsten::PolyExp(dt, y0_[0], 0, 0, 0, false, a, alpha_, 2) +
          torsten::PolyExp(dt, 0, rate_[0], dt, 0, false, a, alpha_, 2);
      }

      if ((y0_[1] != 0) || (rate_[1] != 0)) {
        a[0] = 1;
        pred(0, 1) += torsten::PolyExp(dt, y0_[1], 0, 0, 0, false, a, alpha_, 1) +
          torsten::PolyExp(dt, 0, rate_[1], dt, 0, false, a, alpha_, 1);
      }
      return pred;
    }

    /*
     * Solve the transient problem and return the result in
     * form of data, arranged as (solution value, grad1, grad2...)
     */
    Eigen::VectorXd solve_d(const T_time& t_next) const {
      return torsten::model_solve_d(*this, t_next);
    }

  /**
   * Solve one-cpt model: steady state solution
   *
   * @tparam T_time dosing interval time type
   * @tparam T_model ODE model type
   * @tparam T_amt dosing amount type
   */
    template<typename T_amt, typename T_r, typename T_ii>
    Eigen::Matrix<typename stan::return_type<T_par, T_amt, T_r, T_ii>::type, Eigen::Dynamic, 1>
    solve(const T_amt& amt, const T_r& rate, const T_ii& ii, const int& cmt) const { //NOLINT
      using Eigen::Matrix;
      using Eigen::Dynamic;
      using std::vector;

      using ss_scalar_type = typename stan::return_type<T_par, T_amt, T_r, T_ii>::type;

      const double inf = std::numeric_limits<double>::max();  // "infinity"

      stan::math::check_positive_finite("steady state one-cpt solver", "cmt", cmt);
      stan::math::check_less("steady state one-cpt solver", "cmt", cmt, 3);

      std::vector<ss_scalar_type> a(2, 0);
      Matrix<ss_scalar_type, 1, Dynamic> pred = Matrix<ss_scalar_type, 1, Dynamic>::Zero(2);
      if (stan::math::value_of(rate) == 0) {  // bolus dose
        if (cmt == 1) {
          a[0] = 0;
          a[1] = 1;
          pred(0) = torsten::PolyExp(ii, amt, 0, 0, ii, true, a, alpha_, 2);
          a[0] = ka_ / (ka_ - alpha_[0]);
          a[1] = -a[0];
          pred(1) = torsten::PolyExp(ii, amt, 0, 0, ii, true, a, alpha_, 2);
        } else {  // cmt=2
          a[0] = 1;
          pred(1) = torsten::PolyExp(ii, amt, 0, 0, ii, true, a, alpha_, 1);
        }
      } else if (ii > 0) {  // multiple truncated infusions
        double delta = torsten::unpromote(amt / rate);
        static const char* function("Steady State Event");
        torsten::check_mti(amt, delta, ii, function);

        if (cmt == 1) {
          a[0] = 0;
          a[1] = 1;
          pred(0) = torsten::PolyExp(ii, 0, rate, amt / rate, ii, true, a, alpha_, 2);
          a[0] = ka_ / (ka_ - alpha_[0]);
          a[1] = -a[0];
          pred(1) = torsten::PolyExp(ii, 0, rate, amt / rate, ii, true, a, alpha_, 2);
        } else {  // cmt = 2
          a[0] = 1;
          pred(1) = torsten::PolyExp(ii, 0, rate, amt / rate, ii, true, a, alpha_, 1);
        }
      } else {  // constant infusion
        if (cmt == 1) {
          a[0] = 0;
          a[1] = 1;
          pred(0) = torsten::PolyExp(0, 0, rate, inf, 0, true, a, alpha_, 2);
          a[0] = ka_ / (ka_ - alpha_[0]);
          a[1] = -a[0];
          pred(1) = torsten::PolyExp(0, 0, rate, inf, 0, true, a, alpha_, 2);
        } else {  // cmt = 2
          a[0] = 1;
          pred(1) = torsten::PolyExp(0, 0, rate, inf, 0, true, a, alpha_, 1);
        }
      }
      return pred;
    }

    /*
     * Solve the steady-state problem and return the result in
     * form of data, arranged as (solution value, grad1, grad2...)
     */
    template<typename T_amt, typename T_r, typename T_ii>
    Eigen::VectorXd solve_d(const T_amt& amt, const T_r& rate, const T_ii& ii, const int& cmt) const {
      return torsten::model_solve_d(*this, amt, rate, ii, cmt);
    }

    /*
     * wrapper to fit @c PrepWrapper's call signature
     */
    template<torsten::PMXOdeIntegratorId It, typename T_amt, typename T_r, typename T_ii>
    Eigen::Matrix<scalar_type, Eigen::Dynamic, 1>
    solve(const T_amt& amt, const T_r& rate, const T_ii& ii, const int& cmt,
          const torsten::PMXOdeIntegrator<It>& integrator) const {
      return solve(amt, rate, ii, cmt);
    }
  };

  template<typename T_time, typename T_init, typename T_rate, typename T_par>
  constexpr int PMXOneCptModel<T_time, T_init, T_rate, T_par>::Ncmt;

  template<typename T_time, typename T_init, typename T_rate, typename T_par>
  constexpr int PMXOneCptModel<T_time, T_init, T_rate, T_par>::Npar;

  template<typename T_time, typename T_init, typename T_rate, typename T_par>
  constexpr PMXOneCptODE PMXOneCptModel<T_time, T_init, T_rate, T_par>::f_;



}

#endif
