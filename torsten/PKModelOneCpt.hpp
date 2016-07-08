// version 0.8

/**
 * Computes the predicted amounts in each compartment at each event
 * for a one compartment model with first oder absorption. 
 * *
 * @tparam T0 type of scalars for the model parameters.
 * @tparam T1 type of scalar for time of events. 
 * @tparam T2 type of scalar for amount at each event.
 * @tparam T3 type of scalar for rate at each event.
 * @tparam T4 type of scalar for inter-dose inteveral at each event.
 * @param[in] pMatrix parameters at each event
 * @param[in] time times of events  
 * @param[in] amt amount at each event
 * @param[in] rate rate at each event
 * @param[in] ii inter-dose interval at each event
 * @param[in] evid event identity: 
 *                    (0) observation 
 *                    (1) dosing
 *                    (2) other 
 *                    (3) reset 
 *                    (4) reset AND dosing 
 * @param[in] cmt compartment number at each event 
 * @param[in] addl additional dosing at each event 
 * @param[in] ss steady state approximation at each event (0: no, 1: yes)
 * @return a matrix with predicted amount in each compartment 
 *         at each event. 
 */

#include <stan/model/model_header.hpp>
#include <Eigen/Dense>
#include "PKModel/PKModel.hpp"
#include <boost/math/tools/promotion.hpp>

using std::vector;
using Eigen::Dynamic;
using Eigen::Matrix;
using boost::math::tools::promote_args;


template <typename T0, typename T1, typename T2, typename T3, typename T4> 
Matrix <typename promote_args<T0, T1, T2, T3, T4>::type, Dynamic, Dynamic> 
PKModelOneCpt(const vector< Matrix<T0, Dynamic, 1> >& pMatrix, 
			  const vector<T1>& time,
			  const vector<T2>& amt,
			  const vector<T3>& rate,
			  const vector<T4>& ii,
			  const vector<int>& evid,
			  const vector<int>& cmt,
			  const vector<int>& addl,
			  const vector<int>& ss) 
{  
	PKModel model("OneCptModel"); 
  static const char* function("PKModelOneCpt"); 
  Matrix <typename promote_args<T0, T1, T2, T3, T4>::type, Dynamic, Dynamic> pred;
  
	pmetricsCheck(pMatrix, time, amt, rate, ii, evid, cmt, addl, ss, function, model);
	for(int i=0; i<pMatrix.size(); i++) {
	  stan::math::check_positive_finite(function, "PK parameter CL", pMatrix[i](0,0));
	  stan::math::check_positive_finite(function, "PK parameter V2", pMatrix[i](1,0));
	  stan::math::check_positive_finite(function, "PK parameter ka", pMatrix[i](2,0));
	}
	
  //Construct Pred functions for the model.
  Pred1_structure new_Pred1("OneCptModel");
  PredSS_structure new_PredSS("OneCptModel");
  Pred1 = new_Pred1;
  PredSS = new_PredSS; 

  pred = Pred(pMatrix, time, amt, rate, ii, evid, cmt, addl, ss, model, dummy_ode());
    
  return pred;
}