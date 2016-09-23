#ifndef STAN_MATH_TORSTEN_PKMODEL_PRED_PREDSS_LINCPT_HPP
#define STAN_MATH_TORSTEN_PKMODEL_PRED_PREDSS_LINCPT_HPP

#include <iostream>

using std::vector;
using boost::math::tools::promote_args;
using Eigen::Matrix;
using Eigen::Dynamic; 

/**
 * General compartment model using built-in ODE solver. 
 * Calculate amount in each compartment at the end of a steady-state dosing interval
 * or during a steady-state constant input (if ii=0)
 * 
 * 
 * Model using numerical ODE solver are not available for steady state. If data contains a 
 * steady state event, abort. DEV - use invalid / error message 
 * 
 *	 @tparam T_time type of scalar for time
 *	 @tparam T_amt type of scalar for amount
 *	 @tparam T_rate type of scalar for rate
 *	 @tparam T_ii type of scalar for interdose interval
 *	 @tparam T_system type of elements of Matrix that describes ODE 
 *	 @param[in] parameter model parameters at current event
 *	 @param[in] rate
 *	 @param[in] ii interdose interval
 *	 @param[in] cmt compartment in which the event occurs 
 *	 @param[in] f functor for base ordinary differential equation that defines 
 *              compartment model 
 *   @return an eigen vector that contains predicted amount in each compartment 
 *           at the current event. 
 */
template<typename T_amt, typename T_rate, typename T_ii, typename T_system>
Matrix<typename promote_args<T_amt, T_rate, T_ii, T_system>::type, 1, Dynamic>
PredSS_linCpt(const T_amt& amt, 
		   	  const T_rate& rate,
		   	  const T_ii& ii, 
		   	  const int& cmt,
		   	  const Eigen::Matrix<T_system, Eigen::Dynamic,
		   	    Eigen::Dynamic> system) {
  typedef typename promote_args<T_amt, T_rate, T_ii, T_system>::type scalar; 

  std::cout << "ERROR: LinCptModel function does not currently" 
            << " handle Steady State events." 
             << std::endl;
  abort();

  Matrix<scalar, 1, Dynamic> pred = Matrix<scalar, 1, Dynamic>::Zero(3);
  return pred;

}

#endif
