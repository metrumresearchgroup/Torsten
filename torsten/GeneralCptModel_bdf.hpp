// version 0.8

/**
 * Computes the predicted amounts in each compartment at each event
 * for a general compartment model, defined by a system of ordinary
 * differential equations. Uses the stan::math::integrate_ode_bdf 
 * function. 
 *
 * <b>Warning:</b> This prototype does not handle steady state events. 
 *
 * @tparam T0 type of scalars for the model parameters.
 * @tparam T1 type of scalar for time of events. 
 * @tparam T2 type of scalar for amount at each event.
 * @tparam T3 type of scalar for rate at each event.
 * @tparam T4 type of scalar for inter-dose inteveral at each event.
 * @tparam F type of ODE system function.
 * @param[in] f functor for base ordinary differential equation that defines 
 *            compartment model.
 * @param[in] nCmt number of compartments in model
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
 * @param[in] rel_tol relative tolerance passed to CVODE 
 * @param[in] abs_tol absolute tolerance passed to CVODE
 * @param[in] max_num_steps maximal number of admissable steps 
 *            between time-points
 * @return a matrix with predicted amount in each compartment 
 *         at each event. 
 */

#include <stan/model/model_header.hpp>
#include <Eigen/Dense>
#include <boost/math/tools/promotion.hpp>
#include "PKModel/PKModel.hpp"
using std::vector;
using Eigen::Dynamic;
using Eigen::Matrix;
using boost::math::tools::promote_args;


template <typename T0, typename T1, typename T2, typename T3, typename T4, typename F> 
Matrix <typename promote_args<T0, T1, T2, T3, T4>::type, Dynamic, Dynamic> 
generalCptModel_bdf(const F& f,
                    const int nCmt,
			              const vector< Matrix<T0, Dynamic, 1> >& pMatrix, 
			              const vector<T1>& time,
			              const vector<T2>& amt,
			              const vector<T3>& rate,
			              const vector<T4>& ii,
			              const vector<int>& evid,
			              const vector<int>& cmt,
			              const vector<int>& addl,
			              const vector<int>& ss,
			              double rel_tol = 1e-10,
                    double abs_tol = 1e-10,
                    long int max_num_steps = 1e8)
{
	int nParameters, F1Index, tlag1Index;
	nParameters = pMatrix[0].rows();
	F1Index = nParameters - 2*nCmt;
	tlag1Index = nParameters - nCmt;
	PKModel model(nParameters, F1Index, tlag1Index, nCmt);
	static const char* function("generalCptModel_bdf");
	Matrix <typename promote_args<T0, T1, T2, T3, T4>::type, Dynamic, Dynamic> pred;
	
	pmetricsCheck(pMatrix, time, amt, rate, ii, evid, cmt, addl, ss, function, model);

  // define functors used in Pred()
  Pred1_structure new_Pred1("GeneralCptModel");
  PredSS_structure new_PredSS("GeneralCptModel");
  Pred1 = new_Pred1;
  PredSS = new_PredSS; // WARNING: PredSS returns an error message and aborts program
  pmetrics_solver_structure new_pmetrics_solver(rel_tol, abs_tol, max_num_steps, "bdf");
  pmetrics_solver = new_pmetrics_solver;
        
  pred = Pred(pMatrix, time, amt, rate, ii, evid, cmt, addl, ss, model, f);
        
  return pred;
}

