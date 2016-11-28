#ifndef STAN_MATH_TORSTEN_PKMODEL_PRED_PRED1_LINCPT_HPP
#define STAN_MATH_TORSTEN_PKMODEL_PRED_PRED1_LINCPT_HPP

#include <stan/math/prim/mat.hpp>
#include <stan/math/prim/mat/fun/matrix_exp.hpp>
#include <iostream>
#include <vector>

/**
 * Linear compartment model.
 * Calculates the amount in each compartment at dt time units after the time
 * of the initial condition.
 *
 * If the initial time equals the time of the event, than the code does
 * not run the ode integrator, and sets the predicted amount equal to the
 * initial condition. This can happen when we are dealing with events that
 * occur simultaneously. The change to the predicted amount caused by bolus
 * dosing events is handled later in the main Pred function.
 *
 * @tparam T_time type of scalar for time
 * @tparam T_rate type of scalar for rate
 * @tparam T_parameters type of scalar for model parameters
 * @tparam T_system type of element in system matrix
 * @param[in] dt time between current and previous event
 * @param[in] parameter model parameters at current event
 * @param[in] init amount in each compartment at previous event
 * @param[in] rate rate in each compartment
 * @param[in] system matrix describing linear ODE system
 * @return an eigen vector that contains predicted amount in each compartment
 *   at the current event.
 */
template<typename T_time, typename T_parameters, typename T_system,
  typename T_rate>
Eigen::Matrix<typename boost::math::tools::promote_args< T_time, T_system,
  T_rate>::type, 1, Eigen::Dynamic>
Pred1_linCpt(const T_time& dt,
             const ModelParameters<T_time, T_parameters, T_system>& parameter,
             const Eigen::Matrix<typename boost::math::tools::promote_args<
               T_time, T_system, T_rate>::type, 1, Eigen::Dynamic>& init,
             const std::vector<T_rate>& rate) {
  using boost::math::tools::promote_args;
  using Eigen::Matrix;
  using Eigen::Dynamic;
  using stan::math::matrix_exp;
  using stan::math::mdivide_left;
  using stan::math::multiply;

  typedef typename promote_args<T_time, T_system, T_rate>::type scalar;

  if (dt == 0) { return init;
  } else {
    Matrix<T_system, Dynamic, Dynamic> system = parameter.get_K();

    bool rate_zeros = true;
    for (int i = 0; i < rate.size(); i++)
    if (rate[i] != 0) rate_zeros = false;

    if (rate_zeros) {
      Matrix<scalar, Dynamic, Dynamic> dt_system = multiply(dt, system);
      Matrix<scalar, Dynamic, 1> pred = matrix_exp(dt_system)
        * init.transpose();
      return pred.transpose();
    } else {
      int nCmt = system.cols();
      Matrix<scalar, Dynamic, 1> rate_vec(rate.size()), x(nCmt), x2(nCmt);
      for (int i = 0; i < rate.size(); i++) rate_vec(i) = rate[i];
      x = mdivide_left(system, rate_vec);
      x2 = x + init.transpose();
      Matrix<scalar, Dynamic, Dynamic> dt_system = multiply(dt, system);
      Matrix<scalar, Dynamic, 1> pred = matrix_exp(dt_system) * x2;
      pred -= x;
      return pred.transpose();
    }
  }
}

#endif
