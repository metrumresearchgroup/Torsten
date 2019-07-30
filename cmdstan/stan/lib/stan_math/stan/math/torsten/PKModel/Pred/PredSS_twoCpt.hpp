#ifndef STAN_MATH_TORSTEN_PKMODEL_PRED_PREDSS_TWOCPT_HPP
#define STAN_MATH_TORSTEN_PKMODEL_PRED_PREDSS_TWOCPT_HPP

#include <stan/math/torsten/PKModel/Pred/PolyExp.hpp>
#include <stan/math/torsten/PKModel/ModelParameters.hpp>
#include <stan/math/prim/mat/fun/Eigen.hpp>
#include <boost/math/tools/promotion.hpp>
#include <iostream>
#include <vector>
#include <limits>

namespace torsten {

struct PredSS_twoCpt {
  PredSS_twoCpt() { }

  /**
   * One compartment model with first-order absorption.
   * Calculate amount in each compartment at the end of a steady-state dosing interval
   * or during a steady-state constant input (if ii=0)
   *
   *  If the initial time equals the time of the event, then the code does
   *	not run the ode integrator, and sets the predicted amount equal to the
   *	initial condition. This can happen when we are dealing with events that
   *	occur simultaneously. The change to the predicted amount caused by bolus
   *	dosing events is handled later in the main Pred function.
   *
   *	 @tparam T_time type of scalar for time
   *	 @tparam T_amt type of scalar for amount
   *	 @tparam T_rate type of scalar for rate
   *	 @tparam T_ii type of scalar for interdose interval
   *	 @tparam T_parameters type of scalar for model parameters
   *	 @tparam T_addParm type of scalar for additional model parameters
   *	 @param[in] parameter model parameters at current event
   *	 @param[in] rate
   *	 @param[in] ii interdose interval
   *	 @param[in] cmt compartment in which the event occurs
   *   @return an eigen vector that contains predicted amount in each compartment
   *           at the current event.
   */
  template<typename T_time, typename T_amt, typename T_rate, typename T_ii,
           typename T_parameters, typename T_biovar, typename T_tlag>
  Eigen::Matrix<typename boost::math::tools::promote_args<T_amt, T_rate,
    T_ii, T_parameters>::type, 1, Eigen::Dynamic>
  operator() (const ModelParameters<T_time, T_parameters, T_biovar,
                T_tlag>& parameter,
              const T_amt& amt,
              const T_rate& rate,
              const T_ii& ii,
              const int& cmt) const {
    using std::vector;
    using Eigen::Matrix;
    using Eigen::Dynamic;
    using boost::math::tools::promote_args;

    typedef typename boost::math::tools::promote_args<T_amt, T_rate,
      T_ii, T_parameters>::type scalar;

    double inf = std::numeric_limits<double>::max();  // "infinity"

    T_parameters CL = parameter.get_RealParameters(false)[0],
      Q = parameter.get_RealParameters(false)[1],
      V2 = parameter.get_RealParameters(false)[2],
      V3 = parameter.get_RealParameters(false)[3],
      ka = parameter.get_RealParameters(false)[4],
      k10 = CL / V2,
      k12 = Q / V2,
      k21 = Q / V3,
      ksum = k10 + k12 + k21;

    vector<T_parameters> alpha(3, 0);
    alpha[0] = (ksum + sqrt(ksum * ksum - 4 * k10 * k21)) / 2,
      alpha[1] = (ksum - sqrt(ksum * ksum - 4 * k10 * k21)) / 2,
      alpha[2] = ka;

    vector<scalar> a(3, 0);
    Matrix<scalar, 1, Dynamic> pred = Matrix<scalar, 1, Dynamic>::Zero(3);
    if (rate == 0) {  // bolus dose
      if (cmt == 1) {
        pred(0, 0) = PolyExp(ii, amt, 0, 0, ii, true, a, alpha, 3);
        a[0] = ka * (k21 - alpha[0]) / ((ka - alpha[0])
                                          * (alpha[1] - alpha[0]));
        a[1] = ka * (k21 - alpha[1]) / ((ka - alpha[1])
                                          * (alpha[0] - alpha[1]));
        a[2] = -(a[0] + a[1]);
        pred(0, 1) = PolyExp(ii, amt, 0, 0, ii, true, a, alpha, 3);
        a[0] = ka * k12 / ((ka - alpha[0]) * (alpha[1] - alpha[0]));
        a[1] = ka * k12 / ((ka - alpha[1]) * (alpha[0] - alpha[1]));
        a[2] = -(a[0] + a[1]);
        pred(0, 2) = PolyExp(ii, amt, 0, 0, ii, true, a, alpha, 3);
      } else if (cmt == 2) {
        a[0] = (k21 - alpha[0]) / (alpha[1] - alpha[0]);
        a[1] = (k21 - alpha[1]) / (alpha[0] - alpha[1]);
        pred(0, 1) = PolyExp(ii, amt, 0, 0, ii, true, a, alpha, 2);
        a[0] = ka * k12 / ((ka - alpha[0]) * (alpha[1] - alpha[0]));
        a[1] = ka * k12 / ((ka - alpha[1]) * (alpha[0] - alpha[1]));
        pred(0, 2) = PolyExp(ii, amt, 0, 0, ii, true, a, alpha, 3);
      } else {  // cmt=3
        a[0] = k21 / (alpha[1] - alpha[0]);
        a[1] = -a[0];
        pred(0, 1) = PolyExp(ii, amt, 0, 0, ii, true, a, alpha, 2);
        a[0] = (k10 + k12 - alpha[0]) / (alpha[1] - alpha[0]);
        a[1] = (k10 + k12 - alpha[1]) / (alpha[0] - alpha[1]);
        pred(0, 2) = PolyExp(ii, amt, 0, 0, ii, true, a, alpha, 2);
      }
    } else if (ii > 0) {  // multiple truncated infusions
      double delta = unpromote(amt / rate);
      static const char* function("Steady State Event");
      check_mti(amt, delta, ii, function);

      if (cmt == 1) {
        a[2] = 1;
        pred(0, 0) = PolyExp(ii, 0, rate, amt / rate, ii, true, a, alpha, 3);
        a[0] = ka * (k21 - alpha[0]) / ((ka - alpha[0])
                                          * (alpha[1] - alpha[0]));
        a[1] = ka * (k21 - alpha[1]) / ((ka - alpha[1])
                                          * (alpha[0] - alpha[1]));
        a[2] = - (a[0] + a[1]);
        pred(0, 1) = PolyExp(ii, 0, rate, amt / rate, ii, true, a, alpha, 3);
        a[0] = ka * k12 / ((ka - alpha[0]) * (alpha[1] - alpha[0]));
        a[1] = ka * k12 / ((ka - alpha[1]) * (alpha[0] - alpha[1]));
        a[2] = -(a[0] + a[1]);
        pred(0, 2) = PolyExp(ii, 0, rate, amt / rate, ii, true, a, alpha, 3);
      } else if (cmt == 2) {
        a[0] = (k21 - alpha[0]) / (alpha[1] - alpha[0]);
        a[1] = (k21 - alpha[1]) / (alpha[0] - alpha[1]);
        pred(0, 1) = PolyExp(ii, 0, rate, amt / rate, ii, true, a, alpha, 2);
        a[0] = k12 / (alpha[1] - alpha[0]);
        a[1] = -a[0];
        pred(0, 2) = PolyExp(ii, 0, rate, amt / rate, ii, true, a, alpha, 2);
      } else {  // cmt=3
        a[0] = k21 / (alpha[1] - alpha[0]);
        a[1] = -a[0];
        pred(0, 1) = PolyExp(ii, 0, rate, amt / rate, ii, true, a, alpha, 2);
        a[0] = (k10 + k12 - alpha[0]) / (alpha[1] - alpha[0]);
        a[1] = (k10 + k12 - alpha[1]) / (alpha[0] - alpha[1]);
        pred(0, 2) = PolyExp(ii, 0, rate, amt / rate, ii, true, a, alpha, 2);
      }
    } else {  // constant infusion
      if (cmt == 1) {
        a[2] = 1;
        pred(0, 0) = PolyExp(0, 0, rate, inf, 0, true, a, alpha, 3);
        a[0] = ka * (k21 - alpha[0]) / ((ka - alpha[0])
                                          * (alpha[1] - alpha[0]));
        a[1] = ka * (k21 - alpha[1]) / ((ka - alpha[1])
                                          * (alpha[0] - alpha[1]));
        a[2] = -(a[0] + a[1]);
        pred(0, 1) = PolyExp(0, 0, rate, inf, 0, true, a, alpha, 3);
        a[0] = ka * k12 / ((ka - alpha[0]) * (alpha[1] - alpha[0]));
        a[1] = ka * k12 / ((ka - alpha[1]) * (alpha[0] - alpha[1]));
        a[2] = -(a[0] + a[1]);
        pred(0, 2) = PolyExp(0, 0, rate, inf, 0, true, a, alpha, 3);
      } else if (cmt == 2) {
        a[0] = (k21 - alpha[0]) / (alpha[1] - alpha[0]);
        a[1] = (k21 - alpha[1]) / (alpha[0] - alpha[1]);
        pred(0, 1) = PolyExp(0, 0, rate, inf, 0, true, a, alpha, 2);
        a[0] = k12 / (alpha[1] - alpha[0]);
        a[1] = -a[0];
        pred(0, 2) = PolyExp(0, 0, rate, inf, 0, true, a, alpha, 2);
      } else {  // cmt=3
        a[0] = k21 / (alpha[1] - alpha[0]);
        a[1] = -a[0];
        pred(0, 1) = PolyExp(0, 0, rate, inf, 0, true, a, alpha, 2);
        a[0] = (k10 + k12 - alpha[0]) / (alpha[1] - alpha[0]);
        a[1] = (k10 + k12 - alpha[1]) / (alpha[0] - alpha[1]);
        pred(0, 2) = PolyExp(0, 0, rate, inf, 0, true, a, alpha, 2);
      }
    }
    return pred;
  }
};

}
#endif
