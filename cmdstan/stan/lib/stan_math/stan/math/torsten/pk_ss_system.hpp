#ifndef STAN_MATH_TORSTEN_REFACTOR_SS_SYSTEM_HPP
#define STAN_MATH_TORSTEN_REFACTOR_SS_SYSTEM_HPP

#include <stan/math/prim/mat/fun/to_vector.hpp>
#include <stan/math/torsten/PKModel/Pred/Pred1_general.hpp>
#include <stan/math/torsten/PKModel/functors/check_mti.hpp>
#include <stan/math/torsten/pmx_ode_integrator.hpp>
#include <vector>
#include <iostream>

namespace torsten {

  template <PMXOdeIntegratorId It, typename T_amt, typename T_rate, typename F, typename F2> // NOLINT
    struct SSFunctor;

/**
 * A structure to store the algebraic system
 * which gets solved when computing the steady
 * state solution.
 * 
 * In this structure, both amt and rate are fixed
 * variables.
 */
  template <PMXOdeIntegratorId It, typename F, typename F2>
    struct SSFunctor<It, double, double, F, F2> {
    F f_;
    double ii_;
    int cmt_;  // dosing compartment
    const PMXOdeIntegrator<It> integrator_;
    int nPK_;

    SSFunctor() {}

    SSFunctor(const F& f,
                      double ii,
                      int cmt,
                      const PMXOdeIntegrator<It>& integrator) :
      f_(f), ii_(ii), cmt_(cmt), integrator_(integrator), nPK_(0)
    {}

    SSFunctor(const F& f,
                      double ii,
                      int cmt,
                      const PMXOdeIntegrator<It>& integrator,
                      int nPK) :
      f_(f), ii_(ii), cmt_(cmt), integrator_(integrator), nPK_(nPK)
    {}

  template<typename Ty, typename T_pksolve>
  struct AttachPKSol {
    void attach(const double& delta,
                const int& npk,
                const std::vector<double> &dat_ode,
                const Eigen::Matrix<Ty, Eigen::Dynamic, 1>& y,
                Eigen::Matrix<Ty, Eigen::Dynamic, 1>& y2) {
      const int nParms = y.size() - npk;
      if (npk != 0) {
        refactor::PKRec<Ty> x0_pk(npk);
        for (int i = 0; i < npk; i++) x0_pk(i) = y(nParms + i);

        using T_pkmodel =
          typename T_pksolve::template default_model<double, Ty, double, Ty>;
        const double t00 = 0.0;
        auto y_par = stan::math::to_array_1d(y);
        T_pkmodel pkmodel {t00, x0_pk, dat_ode, y_par};
        refactor::PKRec<Ty> x_pk = T_pksolve().solve(pkmodel, delta);
        for (int i = 0; i < npk; i++) y2(nParms + i) = x_pk(i);
      }
    }
  };

  template<typename Ty>
  struct AttachPKSol<Ty, void> {
    void attach(const double& delta,
                const int& npk,
                const std::vector<double> &dat_ode,
                const Eigen::Matrix<Ty, Eigen::Dynamic, 1>& y,
                Eigen::Matrix<Ty, Eigen::Dynamic, 1>& y2)
    {}
  };

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

      AttachPKSol<T1, F2>().attach(delta, nPK_, dat_ode, y, y2);

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
  template <PMXOdeIntegratorId It, typename F>
    struct SSFunctor<It, stan::math::var, double, F, void> {
    F f_;
    double ii_;
    int cmt_;  // dosing compartment
    const PMXOdeIntegrator<It> integrator_;
    int nPK_;

    SSFunctor() {}

  SSFunctor(const F& f,
                    double ii,
                    int cmt,
                    const PMXOdeIntegrator<It>& integrator) :
    f_(f), ii_(ii), cmt_(cmt), integrator_(integrator), nPK_(0)
  {}

  SSFunctor(const F& f,
                    double ii,
                    int cmt,
                    const PMXOdeIntegrator<It>& integrator,
                    int nPK) :
    f_(f), ii_(ii), cmt_(cmt), integrator_(integrator), nPK_(nPK)
  {}

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
