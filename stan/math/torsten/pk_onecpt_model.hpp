#ifndef STAN_MATH_TORSTEN_ONECPT_MODEL_HPP
#define STAN_MATH_TORSTEN_ONECPT_MODEL_HPP

#include <stan/math/torsten/torsten_def.hpp>

namespace refactor {

  using boost::math::tools::promote_args;

  /**
   * standard one compartment PK ODE functor
   */
  struct PKOneCptODE {
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
  class PKOneCptModel {
    const T_time &t0_;
    const Eigen::Matrix<T_init, 1, Eigen::Dynamic>& y0_;
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
    static constexpr PKOneCptODE f_ = PKOneCptODE();

    using scalar_type = typename
      stan::return_type<T_time, T_init, T_rate, T_par>::type;
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
    PKOneCptModel(const T_time& t0,
                  const Eigen::Matrix<T_init, 1, Eigen::Dynamic>& y0,
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
    {}

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
    PKOneCptModel(const T_time& t0,
                  const Eigen::Matrix<T_init, 1, Eigen::Dynamic>& y0,
                  const std::vector<T_rate> &rate,
                  const std::vector<T_par> & par,
                  const T_mp<Ts...> &parameter) :
      PKOneCptModel(t0, y0, rate, par.at(0), par.at(1), par.at(2))
    {}

  /**
   * One-compartment PK model constructor
   *
   * @param t0 initial time
   * @param y0 initial condition
   * @param rate dosing rate
   * @param par model parameters
   */
    PKOneCptModel(const T_time& t0,
                  const Eigen::Matrix<T_init, 1, Eigen::Dynamic>& y0,
                  const std::vector<T_rate> &rate,
                  const std::vector<T_par> & par) :
      PKOneCptModel(t0, y0, rate, par.at(0), par.at(1), par.at(2))
    {}

  /**
   * One-compartment PK model get methods
   */
    const T_time              & t0()      const { return t0_;    }
    const PKRec<T_init>       & y0()      const { return y0_;    }
    const std::vector<T_rate> & rate()    const { return rate_;  }
    const std::vector<T_par>  & alpha()   const { return alpha_; }
    const std::vector<T_par>  & par()     const { return par_;   }
    const PKOneCptODE         & f()       const { return f_;     }
    const int                 & ncmt ()   const { return Ncmt;   }

    /**
     * Solve one-cpt model using analytical solution.
     *
     * @tparam T_time time type
     * @tparam T_model ODE model type
     */
    Eigen::Matrix<scalar_type, Eigen::Dynamic, 1>
    solve(const T_time& dt) {
      using Eigen::Matrix;
      using Eigen::Dynamic;

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

  /**
   * Solve one-cpt model: steady state solution
   *
   * @tparam T_time dosing interval time type
   * @tparam T_model ODE model type
   * @tparam T_amt dosing amount type
   */
    template<typename T_amt>
    Eigen::Matrix<scalar_type, Eigen::Dynamic, 1>
    solve(const T_amt& amt,
          const T_time& ii,
          const int& cmt) {
      using Eigen::Matrix;
      using Eigen::Dynamic;
      using std::vector;

      const double inf = std::numeric_limits<double>::max();  // "infinity"

      stan::math::check_positive("steady state one-cpt solver", "cmt", cmt);
      stan::math::check_less("steady state one-cpt solver", "cmt", cmt, 3);

      const auto rate = rate_.at(cmt - 1);

      std::vector<scalar_type> a(2, 0);
      Matrix<scalar_type, 1, Dynamic> pred = Matrix<scalar_type, 1, Dynamic>::Zero(2);
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

  };

  template<typename T_time, typename T_init, typename T_rate, typename T_par>
  constexpr int PKOneCptModel<T_time, T_init, T_rate, T_par>::Ncmt;

  template<typename T_time, typename T_init, typename T_rate, typename T_par>
  constexpr int PKOneCptModel<T_time, T_init, T_rate, T_par>::Npar;

  template<typename T_time, typename T_init, typename T_rate, typename T_par>
  constexpr PKOneCptODE PKOneCptModel<T_time, T_init, T_rate, T_par>::f_;



}

#endif
