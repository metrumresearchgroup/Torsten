#ifndef STAN_MATH_TORSTEN_PKMODEL_FUNCTORS_SS_SYSTEM2_HPP
#define STAN_MATH_TORSTEN_PKMODEL_FUNCTORS_SS_SYSTEM2_HPP

#include <stan/math/prim/mat/fun/to_vector.hpp>
#include <stan/math/torsten/PKModel/Pred/Pred1_general.hpp>
#include <stan/math/torsten/PKModel/functors/check_mti.hpp>
#include <stan/math/rev/core.hpp>
#include <stan/math/fwd/core.hpp>
#include <vector>
#include <iostream>

namespace torsten {

/**
 * A structure to store the algebraic system
 * which gets solved when computing the steady
 * state solution.
 * 
 * In this structure, both amt and rate are fixed
 * variables.
 */
  template <typename F, typename F2, typename T_integrator>
struct SS_system2_dd {
  F f_;
  F2 f2_;
  double ii_;
  int cmt_;  // dosing compartment
  const T_integrator integrator_;
  int nPK_;

  SS_system2_dd() {}

  SS_system2_dd(const F& f,
               const F2& f2,
               double ii,
               int cmt,
               const T_integrator& integrator)
    : f_(f), f2_(f2), ii_(ii), cmt_(cmt), integrator_(integrator),
      nPK_(0) { }

  SS_system2_dd(const F& f,
               const F2& f2,
               double ii,
               int cmt,
               const T_integrator& integrator,
               int nPK)
    : f_(f), f2_(f2), ii_(ii), cmt_(cmt), integrator_(integrator),
      nPK_(nPK) { }

  /**
   *  dd regime.
   *  dat contains the rates in each compartment followed
   *  by the adjusted amount (biovar * amt).
   */
  template <typename T0, typename T1>
  inline
  Eigen::Matrix<typename boost::math::tools::promote_args<T0, T1>::type,
                Eigen::Dynamic, 1>
  operator()(const Eigen::Matrix<T0, Eigen::Dynamic, 1>& x,
             const Eigen::Matrix<T1, Eigen::Dynamic, 1>& y,
             const std::vector<double>& dat,
             const std::vector<int>& dat_int,
             std::ostream* msgs) const {
    using stan::math::to_array_1d;
    using stan::math::to_vector;
    using stan::math::to_vector;
    using Eigen::Matrix;
    using Eigen::Dynamic;
    using std::vector;

    typedef typename boost::math::tools::promote_args<T0, T1>::type scalar;
    typedef typename stan::return_type<T0, T1>::type T_deriv;

    double t0 = 0;
    vector<double> ts(1);

    vector<scalar> x0(x.size());
    for (size_t i = 0; i < x0.size(); i++) x0[i] = x(i);
    double amt = dat[dat.size() - 1];
    double rate = dat[cmt_ - 1];

    // real data, which gets passed to the integrator, shoud not have
    // amt in it. Important for the mixed solver where the last element
    // is expected to be the absolute time (in this case, 0).
    vector<double> dat_ode = dat;
    dat_ode.pop_back();

    Eigen::Matrix<scalar, Eigen::Dynamic, 1> result(x.size());

    if (rate == 0) {  // bolus dose
      if ((cmt_ - nPK_) >= 0) x0[cmt_ - nPK_ - 1] += amt;
      ts[0] = ii_;

      vector<scalar> pred = integrator_(f_, x0, t0, ts, to_array_1d(y),
                                        dat_ode, dat_int)[0];

      for (int i = 0; i < result.size(); i++)
        result(i) = x(i) - pred[i];

    } else if (ii_ > 0) {  // multiple truncated infusions
      double delta = amt / rate;

      static const char* function("Steady State Event");
      check_mti(amt, delta, ii_, function);

      vector<scalar> pred;
      ts[0] = delta;  // time at which infusion stops
      x0 = integrator_(f_, to_array_1d(x), t0, ts, to_array_1d(y),
                       dat_ode, dat_int)[0];

      Matrix<T1, Dynamic, 1> y2(y.size());
      int nParms = y.size() - nPK_;
      for (int i = 0; i < nParms; i++) y2(i) = y(i);

      if (nPK_ != 0) {
        Matrix<T1, 1, Dynamic> x0_pk(nPK_);
        for (int i = 0; i < nPK_; i++) x0_pk(i) = y(nParms + i);

        Matrix<T1, 1, Dynamic>
          x_pk = f2_(delta,
                ModelParameters<double, T1, double, double>
                  (0, to_array_1d(y), vector<double>(),
                   vector<double>()),
                   x0_pk, dat_ode);

        for (int i = 0; i < nPK_; i++) y2(nParms + i) = x_pk(i);
      }

      ts[0] = ii_ - delta;
      dat_ode[cmt_ - 1] = 0;
      pred = integrator_(f_, x0, t0, ts, to_array_1d(y2), dat_ode, dat_int)[0];

      for (int i = 0; i < result.size(); i++)
        result(i) = x(i) - pred[i];
    } else {  // constant infusion
      vector<T_deriv> derivative = f_(0, to_array_1d(x), to_array_1d(y),
                                      dat_ode, dat_int, 0);
      result = to_vector(derivative);
    }

    return result;
  }
};

/**
 * A structure to store the algebraic system
 * which gets solved when computing the steady
 * state solution.
 *
 * In this structure, amt is a random variable
 * and rate a fixed variable (vd regime).
 */
  template <typename F, typename T_integrator>
struct SS_system2_vd {
  F f_;
  double ii_;
  int cmt_;  // dosing compartment
    const T_integrator integrator_;
  int nPK_;

  SS_system2_vd() { }

  SS_system2_vd(const F& f,
               double ii,
               int cmt,
               const T_integrator& integrator)
    : f_(f), ii_(ii), cmt_(cmt), integrator_(integrator),
      nPK_(0) { }

  SS_system2_vd(const F& f,
               double ii,
               int cmt,
               const T_integrator& integrator,
               int nPK)
    : f_(f), ii_(ii), cmt_(cmt), integrator_(integrator),
      nPK_(nPK) { }

 /**
  *  Case where the modified amt is a random variable. This
  *  will usually happen because biovar is a parameter, making 
  *  amt a transformed parameter.
  *  The last element of y is contains amt.
  *  dat stores the rate.
  */
  template <typename T0, typename T1>
  inline
    Eigen::Matrix<typename boost::math::tools::promote_args<T0, T1>::type,
                  Eigen::Dynamic, 1>
  operator()(const Eigen::Matrix<T0, Eigen::Dynamic, 1>& x,
             const Eigen::Matrix<T1, Eigen::Dynamic, 1>& y,
             const std::vector<double>& dat,
             const std::vector<int>& dat_int,
             std::ostream* msgs) const {
    using stan::math::to_array_1d;
    using stan::math::to_vector;
    using std::vector;
    using stan::math::to_vector;
    using stan::math::invalid_argument;

    typedef typename boost::math::tools::promote_args<T0, T1>::type scalar;
    typedef typename stan::return_type<T0, T1>::type T_deriv;

    double t0 = 0;
    vector<double> ts(1);
    vector<double> rate_v(dat.size(), 0);
    for (size_t i = 0; i < rate_v.size(); i++) rate_v[i] = dat[i];

    vector<scalar> x0(x.size());
    for (size_t i = 0; i < x0.size(); i++) x0[i] = x(i);
    scalar amt = y(y.size() - 1);
    double rate = ((cmt_ - 1) >= 0) ? dat[cmt_ - 1] : 0;

    Eigen::Matrix<scalar, Eigen::Dynamic, 1> result(x.size());
    std::vector<scalar> parms(y.size() - 1);
    for (size_t i = 0; i < parms.size(); i++) parms[i] = y(i);

    if (rate == 0) {  // bolus dose
      if ((cmt_ - nPK_) >= 0) x0[cmt_ - nPK_ - 1] += amt;
      ts[0] = ii_;
      vector<scalar> pred = integrator_(f_, x0, t0, ts, parms,
                                        dat, dat_int)[0];

      for (int i = 0; i < result.size(); i++)
        result(i) = x(i) - pred[i];

    } else if (ii_ > 0) {  // multiple truncated infusions
      invalid_argument("Steady State Event",
                       "Current version does not handle the case of",
                       "", " multiple truncated infusions ",
                       "(i.e ii > 0 and rate > 0) when F * amt is a parameter.");  // NOLINT

    } else {  // constant infusion
      if (amt != 0)
        invalid_argument("Steady State Event",
                         "amt should be 0 when specifying a constant",
                         "", " infusion (i.e. rate > 0 and ii = 0).",
                         "");

      vector<T_deriv> derivative = f_(0, to_array_1d(x), parms,
                                      dat, dat_int, 0);
      result = to_vector(derivative);
    }

    return result;
  }
};

}

#endif
