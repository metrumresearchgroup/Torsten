#ifndef STAN_MATH_TORSTEN_PKMODEL_PRED_PRED1_LINCPT_HPP
#define STAN_MATH_TORSTEN_PKMODEL_PRED_PRED1_LINCPT_HPP

#include <iostream>
#include <stan/math/prim/mat/fun/matrix_exp.hpp>

/**
 * Linear compartment model.
 * Calculates the amount in each compartment at dt time units after the time
 * of the initial condition.
 * 
 *  If the initial time equals the time of the event, than the code does 
 *	not run the ode integrator, and sets the predicted amount equal to the 
 *	initial condition. This can happen when we are dealing with events that 
 *	occur simultaneously. The change to the predicted amount caused by bolus 
 *	dosing events is handled later in the main Pred function. 
 *	
 *	 @tparam T_time type of scalar for time
 *	 @tparam T_rate type of scalar for rate
 *	 @tparam T_parameters type of scalar for model parameters
 *   @tparam T_system type of element in system matrix
 *	 @param[in] dt time between current and previous event
 *	 @param[in] parameter model parameters at current event
 *	 @param[in] init amount in each compartment at previous event
 *	 @param[in] rate rate in each compartment
 *   @param[in] system matrix describing linear ODE system
 *   @return an eigen vector that contains predicted amount in each compartment 
 *           at the current event. 
 */
template<typename T_time, typename T_rate, typename T_system>
Matrix<typename promote_args< T_time, T_rate, T_system>::type, 1, Dynamic> 
Pred1_linCpt(const T_time& dt,
		     const Eigen::Matrix<typename 
		       promote_args<T_time, T_rate, T_system>::type, 1,
		       Eigen::Dynamic>& init,
		     const vector<T_rate>& rate,
		     const Eigen::Matrix<T_system, Eigen::Dynamic, Eigen::Dynamic> system) {

  using boost::math::tools::promote_args;
  using Eigen::Matrix;
  using Eigen::Dynamic;

  typedef typename promote_args<T_time, T_rate, T_system>::type scalar;

  Matrix<scalar, 1, Dynamic> dt_rate(rate.size());
  for(int i=0; i < rate.size(); i++) dt_rate(i) = dt * rate[i];
  Matrix<scalar, Dynamic, Dynamic> dt_system = dt * system;
  Matrix<scalar, Dynamic, 1> pred = stan::math::matrix_exp(dt_system)
    * init.transpose();
 
  return pred.transpose() + dt_rate;
}

#endif
