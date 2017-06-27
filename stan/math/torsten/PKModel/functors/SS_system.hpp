#ifndef STAN_MATH_TORSTEN_PKMODEL_FUNCTORS_SS_SYSTEM_HPP
#define STAN_MATH_TORSTEN_PKMODEL_FUNCTORS_SS_SYSTEM_HPP

#include <stan/math/rev/core.hpp>
#include <stan/math/fwd/core.hpp>
#include <vector>
#include <iostream>

/**
 * A structure to store the algebraic system
 * which gets solved when computing the steady
 * state solution.
 */
template <typename F>
struct SS_system {
  F f_;
  double ii_;
  int cmt_;  // dosing compartment
  integrator_structure integrator_;

  SS_system() { };

  SS_system(const F& f,
            const double& ii,
            int cmt,
            const integrator_structure& integrator)
   : f_(f), ii_(ii), cmt_(cmt), integrator_(integrator) { }

  /**
   *  Case where amt is fixed data.
   *  dat contains the rates in each compartment followed
   *  by the adjusted amount.
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
    using Eigen::Matrix;
    using Eigen::Dynamic;
    using std::vector;
    using stan::math::to_vector;
    typedef typename boost::math::tools::promote_args<T0, T1>::type scalar;
    typedef typename stan::return_type<T0, T1>::type T_deriv;

    double t0 = 0;
    vector<double> ts(1);
    vector<double> rate_v(dat.size() - 1, 0);
    for (size_t i = 0; i < rate_v.size(); i++) rate_v[i] = dat[i];

    vector<scalar> x0(x.size());
    for (size_t i = 0; i < x0.size(); i++) x0[i] = x(i);
    double amt = dat[dat.size() - 1];
    double rate  = dat[cmt_ - 1];

    Matrix<scalar, Dynamic, 1> result(x.size());

    if (rate == 0) {  // bolus dose
      x0[cmt_ - 1] += amt;
      ts[0] = ii_;
      vector<scalar> pred = integrator_(f_, x0, t0, ts, to_array_1d(y),
                                        dat, dat_int)[0];

      for (int i = 0; i < result.size(); i++)
        result(i) = x(i) - pred[i];

    } else if (ii_ > 0) {  // multiple truncated infusions
      double delta = amt / rate;
      vector<scalar> pred;
      
      if (delta < ii_) {
        // In the case where the duration of the infusion is less
        // than the dosing interval, we can do the calculation without
        // using any discrete variables.
        ts[0] = delta;  // time at which infusion stops
        x0 = integrator_(f_, to_array_1d(x), t0, ts, to_array_1d(y),
                         dat, dat_int)[0];
        ts[0] = ii_ - delta;
        vector<double> rate_v(dat.size(), 0);
        pred = integrator_(f_, x0, t0, ts, to_array_1d(y), rate_v, dat_int)[0];
      } else {
        int N = trunc(delta / ii_) + 1;  // number of overlapping rates
        ts[0] = delta - (N - 1) * ii_;  // time at which the oldest infusion dies
        vector<double> rate_v(dat.size());
        for (size_t i = 0; i < rate_v.size(); i++)
          rate_v[i] = N * dat[i];  // compute superposition of rates

        x0 = integrator_(f_, to_array_1d(x), t0, ts, to_array_1d(y),
                         rate_v, dat_int)[0];

        ts[0] = ii_ - ts[0];
        for (size_t i = 0; i < rate_v.size(); i++)
          rate_v[i] = (N - 1) * dat[i];
  
        pred = integrator_(f_, to_array_1d(x), t0, ts,
                           to_array_1d(y), rate_v, dat_int)[0];

      }
      for (int i = 0; i < result.size(); i++)
        result(i) = x(i) - pred[i];

    } else {  // constant infusion
      vector<T_deriv> derivative = f_(0, to_array_1d(x), to_array_1d(y),
                                      dat, dat_int, 0);
      result = to_vector(derivative);
    }

    return result;
  }
};

#endif
