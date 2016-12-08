#ifndef STAN_MATH_TORSTEN_LINCPTMODEL_HPP
#define STAN_MATH_TORSTEN_LINCPTMODEL_HPP

#include <Eigen/Dense>
#include <boost/math/tools/promotion.hpp>
#include <stan/math/torsten/PKModel/PKModel.hpp>
#include <stan/math/prim/mat/err/check_square.hpp>
#include <vector>

/**
 * Computes the predicted amounts in each compartment at each event
 * for a compartment model, described by a linear system of ordinary
 * differential equations. Uses the stan::math::matrix_exp 
 * function.
 *
 * <b>Warning:</b> This prototype does not handle steady state events. 
 *
 * @tparam T0 type of scalar for matrix describing linear ODE system.
 * @tparam T1 type of scalars for the model parameters (F1 and time lag).
 * @tparam T2 type of scalar for time of events. 
 * @tparam T3 type of scalar for amount at each event.
 * @tparam T4 type of scalar for rate at each event.
 * @tparam T5 type of scalar for inter-dose inteveral at each event.
 * @param[in] system square matrix describing the linear system of ODEs
 * @param[in] pMatrix parameters (F1 and time lag) at each event
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
 * between time-points
 * @return a matrix with predicted amount in each compartment 
 * at each event.
 */
template <typename T0, typename T1, typename T2, typename T3,
  typename T4, typename T5>
Eigen::Matrix <typename boost::math::tools::promote_args<T0, T1, T2, T3,
  T4>::type, Eigen::Dynamic, Eigen::Dynamic>
linCptModel(const std::vector< Eigen::Matrix<T0, Eigen::Dynamic,
              Eigen::Dynamic> >& system,
            const std::vector<std::vector<T1> >& pMatrix,
            const std::vector<T2>& time,
            const std::vector<T3>& amt,
            const std::vector<T4>& rate,
            const std::vector<T5>& ii,
            const std::vector<int>& evid,
            const std::vector<int>& cmt,
            const std::vector<int>& addl,
            const std::vector<int>& ss) {
  using std::vector;
  using Eigen::Dynamic;
  using Eigen::Matrix;
  using boost::math::tools::promote_args;

  static const char* function("linCptModel");
  for (int i = 0; i < system.size(); i++)
    stan::math::check_square(function, "system matrix", system[i]);
  int nCmt = system[0].cols();

  int nParameters = pMatrix[0].size(),
    F1Index = nParameters - 2*nCmt,
    tlag1Index = nParameters - nCmt;
  PKModel model(nParameters, F1Index, tlag1Index, nCmt);

  // Check arguments
  pmetricsCheck(pMatrix, time, amt, rate, ii, evid, cmt, addl, ss,
    function, model);

  // define functors used in Pred()
  Pred1_structure new_Pred1("linCptModel");
  PredSS_structure new_PredSS("linCptModel");
  Pred1 = new_Pred1;
  PredSS = new_PredSS;

  Matrix <typename promote_args<T0, T1, T2, T3, T4>::type, Dynamic,
    Dynamic> pred;
  pred = Pred(pMatrix, time, amt, rate, ii, evid, cmt, addl, ss, model,
    dummy_ode(), system);

  return pred;
}

/*
 * Overload function to allow user to pass an std::vector for 
 * pMatrix.
 */
template <typename T0, typename T1, typename T2, typename T3,
  typename T4, typename T5>
Eigen::Matrix <typename boost::math::tools::promote_args<T0, T1, T2, T3,
  T4>::type, Eigen::Dynamic, Eigen::Dynamic>
linCptModel(const std::vector< Eigen::Matrix<T0, Eigen::Dynamic,
              Eigen::Dynamic> >& system,
            const std::vector<T1>& pMatrix,
            const std::vector<T2>& time,
            const std::vector<T3>& amt,
            const std::vector<T4>& rate,
            const std::vector<T5>& ii,
            const std::vector<int>& evid,
            const std::vector<int>& cmt,
            const std::vector<int>& addl,
            const std::vector<int>& ss) {
  std::vector<std::vector<T1> > vec_pMatrix(1);
  vec_pMatrix[0] = pMatrix;

  return linCptModel(system,
    vec_pMatrix, time, amt, rate, ii, evid, cmt, addl, ss);
}

template <typename T0, typename T1, typename T2, typename T3,
  typename T4, typename T5>
Eigen::Matrix <typename boost::math::tools::promote_args<T0, T1, T2, T3,
  T4>::type, Eigen::Dynamic, Eigen::Dynamic>
linCptModel(const Eigen::Matrix<T0, Eigen::Dynamic,
              Eigen::Dynamic>& system,
            const std::vector<std::vector<T1> >& pMatrix,
            const std::vector<T2>& time,
            const std::vector<T3>& amt,
            const std::vector<T4>& rate,
            const std::vector<T5>& ii,
            const std::vector<int>& evid,
            const std::vector<int>& cmt,
            const std::vector<int>& addl,
            const std::vector<int>& ss) {
  std::vector<Eigen::Matrix<T0, Eigen::Dynamic,
              Eigen::Dynamic> > vec_system(1);
  vec_system[0] = system;

  return linCptModel(vec_system,
    pMatrix, time, amt, rate, ii, evid, cmt, addl, ss);
}

template <typename T0, typename T1, typename T2, typename T3,
  typename T4, typename T5>
Eigen::Matrix <typename boost::math::tools::promote_args<T0, T1, T2, T3,
  T4>::type, Eigen::Dynamic, Eigen::Dynamic>
linCptModel(const Eigen::Matrix<T0, Eigen::Dynamic,
              Eigen::Dynamic>& system,
            const std::vector<T1>& pMatrix,
            const std::vector<T2>& time,
            const std::vector<T3>& amt,
            const std::vector<T4>& rate,
            const std::vector<T5>& ii,
            const std::vector<int>& evid,
            const std::vector<int>& cmt,
            const std::vector<int>& addl,
            const std::vector<int>& ss) {
  std::vector<Eigen::Matrix<T0, Eigen::Dynamic,
              Eigen::Dynamic> > vec_system(1);
  vec_system[0] = system;

  std::vector<std::vector<T1> > vec_pMatrix(1);
  vec_pMatrix[0] = pMatrix;

  return linCptModel(vec_system,
    vec_pMatrix, time, amt, rate, ii, evid, cmt, addl, ss);
}

#endif
