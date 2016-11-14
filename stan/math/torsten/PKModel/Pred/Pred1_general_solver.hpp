#ifndef STAN_MATH_TORSTEN_PKMODEL_PRED_PRED1_GENERAL_SOLVER_HPP
#define STAN_MATH_TORSTEN_PKMODEL_PRED_PRED1_GENERAL_SOLVER_HPP

#include <iostream>
#include <stan/math/torsten/PKModel/Pred/unpromote.hpp>

/**
 *	General compartment model using the built-in ODE solver. 
 *	Calculates the amount in each compartment at dt time units after the time
 *	of the initial condition.
 * 	
 *	If the initial time equals the time of the event, than the code does 
 *	not run the ode integrator, and sets the predicted amount equal to the 
 *	initial condition. This can happen when we are dealing with events that 
 *	occur simultaneously. The change to the predicted amount caused by bolus 
 *	dosing events is handled later in the main Pred function. 
 *	
 *	 @tparam T_time type of scalar for time
 *	 @tparam T_rate type of scalar for rate
 *	 @tparam T_parameters type of scalar for model parameters
 *	 @tparam F type of ODE system function
 *	 @param[in] dt time between current and previous event
 *	 @param[in] parameter model parameters at current event
 *	 @param[in] init amount in each compartment at previous event
 *	 @param[in] rate rate in each compartment
 *	 @param[in] f functor for base ordinary differential equation that defines 
 *              compartment model.
 *   @return an eigen vector that contains predicted amount in each compartment 
 *           at the current event. 
 *	
 *	DEV - figure out how to handle when rate is data or autodiff  
 */
template<typename T_time, typename T_rate, typename T_parameters, 
  typename T_system, typename F>
Matrix<typename promote_args< T_time, T_rate, T_parameters>::type, 1, Dynamic> 
Pred1_general_solver(const T_time& dt,
		  			 const ModelParameters<T_time, T_parameters, T_system>& parameter, 
		 			 const Matrix<typename promote_args<T_time, T_rate,
		 			   T_parameters>::type, 1, Dynamic>& init, 
		  			 const vector<T_rate>& rate,
		  			 const F& f) {
  using std::vector;

  typedef typename promote_args<T_time, T_rate, T_parameters>::type scalar;
  assert(init.cols() == rate.size());

  T_time EventTime = parameter.time; // time of current event
  T_time InitTime = EventTime - dt;  // time of previous event	
  
  // Convert time parameters to fixed data for ODE integrator
  vector<double> EventTime_d = vector<double>(1, double(0));
  EventTime_d[0] = unpromote(EventTime);
  double InitTime_d = unpromote(InitTime);
  vector<double> rate_d = vector<double>(rate.size(), double(0));
  for(int i = 0; i < rate.size(); i++) rate_d[i] = unpromote(rate[i]);
		
  vector<T_parameters> theta = parameter.RealParameters;
  vector<scalar> init_vector = vector<scalar>(init.cols(), scalar(0));
  for(int i = 0; i < init_vector.size(); i++) init_vector[i] = init(0, i);

  Matrix<scalar, 1, Dynamic> pred;
  if(EventTime_d[0] == InitTime_d) pred = init;
  else {
    vector< vector<scalar> > pred_V;
	vector<int> idummy;
	pred_V = pmetrics_solver(f, init_vector, InitTime_d,  
							 EventTime_d, theta, rate_d, 
							 idummy);     
						
	//Convert vector in row-major vector (eigen Matrix)
	pred.resize(pred_V[0].size());
	for(int i=0; i < pred_V[0].size(); i++) pred(0,i) = pred_V[0][i];
  }
  return pred;
}

#endif
