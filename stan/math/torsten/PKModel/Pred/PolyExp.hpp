#ifndef STAN_MATH_TORSTEN_PKMODEL_POLYEXP_HPP
#define STAN_MATH_TORSTEN_PKMODEL_POLYEXP_HPP

#include <math.h>
#include <iostream>
#include <limits>
#include <vector>

namespace torsten {
// DEV
// Code is not optimized. Instead of using FOR loops for a and alpha,
// we coud use matrix operations. On the other hand, seeing most models
// have a relatively low number of compartments, efficiency gain would
// probably be very small.

/**
 *	PolyExp calculates portions of analytical solutions to certain ODEs.
 *	It is called by the pred and predSS functions for the linear one and two compartment
 *	models.
 *
 *	In particular, it calculates the sums of exponentials of the form:
 *		sum(i=1 to n) a[i]*exp(-alpha[i]*t)
 *	or the convolution of that with a step function.
 *
 *	First, update dose.
 *	Three cases:
 *		a. bolus dose
 *		b. no bolus dose, unsteady state
 *		c. no bolus dose, steady state
 *
 *	Second, update rate.
 * 	Two cases:
 *		a. truncated infusion
 *		b. continuous infusion
 *
 * @tparam T_x type of scalar for independent variable (often time)
 * @tparam T_dose type of scalar for dose
 * @tparam T_rate type of scalar for rate
 * @tparam T_xinf type of scalar for duration of infusion time
 * @tparam T_tau type of scalar for dosing interval (for multiple dosing case)
 * @tparam T_a type of scalar for scale relative to unit bolus
 * @tparam T_alpha type of scalar for unit
 * @param[in] x independent variable (often time)
 * @param[in] dose
 * @param[in] rate
 * @param[in] xinf duration of infusion time
 * @param[in] tau dosing interval (for multiple dosing case)
 * @param[in] ss steady state approximation (0:no, 1:yes)
 * @param[in] a relative unit to bolus
 * @param[in] alpha unit
 * @param[in] number of terms in polyexponential
 * @return sum of exponentials or convolution of sum with a step function
 *
 */
template<typename T_x, typename T_dose, typename T_rate, typename T_xinf,
  typename T_tau, typename T_a, typename T_alpha>
typename boost::math::tools::promote_args<T_x, T_dose, T_rate,
  typename boost::math::tools::promote_args<T_xinf, T_tau, T_a,
  T_alpha>::type>::type
PolyExp(const T_x& x,
        const T_dose& dose,
        const T_rate& rate,
        const T_xinf& xinf,
        const T_tau& tau,
        const bool& ss,
        const std::vector<T_a>& a,
        const std::vector<T_alpha>& alpha,
        const int& n) {
  using std::vector;
  using boost::math::tools::promote_args;

  typedef typename promote_args<T_x, T_dose, T_rate,
    typename promote_args<T_xinf, T_tau, T_a, T_alpha>::type>::type scalar;

  scalar result = 0, bolusResult, dx, nlntv;
  double inf = std::numeric_limits<double>::max();  // "infinity"

  assert((alpha.size() >= (size_t) n) && (a.size() >= (size_t) n));

  // UPDATE DOSE
  if (dose > 0) {  // bolus dose
    if ((tau <= 0) && (x >= 0)) {
      for (int i = 0; i < n; i++) result += a[i] * exp(-alpha[i] * x);
    } else {
      if (!ss) {
        nlntv = x / tau + 1;
        dx = x - trunc(x / tau) * tau;
        for (int i = 0; i < n; i++)
          result += a[i] * exp(-alpha[i] * x) *
            (1 - exp(-nlntv * alpha[i] * tau))
            / (1 - exp(-alpha[i] * tau));
      } else {
        dx = x - trunc(x / tau) * tau;
        for (int i = 0; i < n; i++)
          result += a[i] * exp(-alpha[i] * x) / (1 - exp(-alpha[i] * tau));
        }
    }
  }
  bolusResult = dose * result;

  // UPDATE RATE
  result = 0;
  if ((rate > 0) && (xinf < inf)) {  // truncated infusion
    if (tau <= 0) {
      if (x >= 0) {
        if (x <= xinf) {
           for (int i = 0; i < n; i++) {
            result += a[i] * (1 - exp(-alpha[i] * x)) / alpha[i];
           }
        } else {
          for (int i = 0; i < n; i++)
            result += a[i] * (1 - exp(-alpha[i] * xinf))
              * exp(-alpha[i] * (x - xinf)) / alpha[i];
        }
      }
    } else if (!ss) {
      assert(xinf <= tau);  // and "other case later", says Bill
      dx = x - trunc(x / tau) * tau;
      nlntv = trunc(x / tau) + 1;
      if (dx <= xinf) {
        for (int i = 0; i < n; i++) {
          if (n > 1)
            result += a[i] * (1 - exp(-alpha[i] * xinf)) * exp(-alpha[i] *
              (dx - xinf + tau)) * (1 - exp(-(nlntv - 1) * alpha[i] * tau)) /
              (1 - exp(-alpha[i] * tau)) / alpha[i];
          result += a[i] * (1 - exp(-alpha[i] * dx)) / alpha[i];
        }
      } else {
        for (int i = 0; i < n; i++)
          result += a[i] * (1 - exp(-alpha[i] * xinf)) * exp(-alpha[i]
            * (dx - xinf)) * (1 - exp(-nlntv * alpha[i] * tau))
            / (1 - exp(-alpha[i] * tau)) / alpha[i];
      }
    } else {
      assert(xinf <= tau);
      dx = x - trunc(x / tau) * tau;
      nlntv = trunc(x / tau) + 1;
      if (dx <= xinf) {
        for (int i = 0; i < n; i++)
          result += a[i] * (1 - exp(-alpha[i] * xinf)) * exp(-alpha[i]
            * (dx - xinf + tau)) / (1 - exp(-alpha[i] * tau)) / alpha[i] + a[i]
            * (1 - exp(-alpha[i] * dx)) / alpha[i];
      } else {
        for (int i = 0; i < n; i++)
          result += a[i] * (1 - exp(-alpha[i] * xinf))
            * exp(-alpha[i] * (dx - xinf))
            / (1 - exp(-alpha[i] * tau)) / alpha[i];
      }
    }
  } else {  // continuous infusion (xinf = inf). tau is ignored.
    if (!ss) {
      if (x >= 0) {
        for (int i = 0; i < n; i++)
          result += a[i] * (1 - exp(-alpha[i] * x)) / alpha[i];
      }
    } else {
      for (int i = 0; i < n; i++)
        result += a[i] / alpha[i];
    }
  }

  return bolusResult + rate * result;
}

}

#endif
