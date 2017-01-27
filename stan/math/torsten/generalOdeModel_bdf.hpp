#ifndef STAN_MATH_TORSTEN_GENERALODEMODEL_BDF_HPP
#define STAN_MATH_TORSTEN_GENERALODEMODEL_BDF_HPP

#include <Eigen/Dense>
#include <stan/math/torsten/PKModel/PKModel.hpp>
#include <boost/math/tools/promotion.hpp>
#include <vector>

/**
 * Computes the predicted amounts in each compartment at each event
 * for a general compartment model, defined by a system of ordinary
 * differential equations. Uses the stan::math::integrate_ode_bdf 
 * function. 
 *
 * <b>Warning:</b> This prototype does not handle steady state events. 
 *
 * @tparam T0 type of scalar for time of events. 
 * @tparam T1 type of scalar for amount at each event.
 * @tparam T2 type of scalar for rate at each event.
 * @tparam T3 type of scalar for inter-dose inteveral at each event.
 * @tparam T4 type of scalars for the model parameters.
 * @tparam T5 type of scalars for the bio-variability parameters.
 * @tparam T6 type of scalars for the model tlag parameters.
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
 * @param[in] rel_tol relative tolerance for the Boost ode solver 
 * @param[in] abs_tol absolute tolerance for the Boost ode solver
 * @param[in] max_num_steps maximal number of steps to take within 
 *            the Boost ode solver 
 * @return a matrix with predicted amount in each compartment 
 *         at each event. 
 */
template <typename T0, typename T1, typename T2, typename T3, typename T4,
          typename T5, typename T6, typename F>
Eigen::Matrix <typename boost::math::tools::promote_args<T0, T1, T2, T3,
  typename boost::math::tools::promote_args<T4, T5, T6>::type>::type,
  Eigen::Dynamic, Eigen::Dynamic>
generalOdeModel_bdf(const F& f,
                    const int nCmt,
                    const std::vector<T0>& time,
                    const std::vector<T1>& amt,
                    const std::vector<T2>& rate,
                    const std::vector<T3>& ii,
                    const std::vector<int>& evid,
                    const std::vector<int>& cmt,
                    const std::vector<int>& addl,
                    const std::vector<int>& ss,
                    const std::vector<std::vector<T4> >& pMatrix,
                    const std::vector<std::vector<T5> >& biovar,
                    const std::vector<std::vector<T6> >& tlag,
                    double rel_tol = 1e-10,
                    double abs_tol = 1e-10,
                    long int max_num_steps = 1e8) {  // NOLINT(runtime/int)
  using std::vector;
  using Eigen::Dynamic;
  using Eigen::Matrix;
  using boost::math::tools::promote_args;

  PKModel model(pMatrix[0].size(), nCmt);

  // check arguments
  static const char* function("generalOdeModel_bdf");
  pmetricsCheck(time, amt, rate, ii, evid, cmt, addl, ss,
                pMatrix, biovar, tlag, function, model);

  // define functors used in Pred()
  Pred1_structure new_Pred1("generalOdeModel");
  PredSS_structure new_PredSS("error");
  Pred1 = new_Pred1;
  PredSS = new_PredSS;  // WARNING: PredSS returns an error
  pmetrics_solver_structure new_pmetrics_solver(rel_tol, abs_tol,
    max_num_steps, "bdf");
  pmetrics_solver = new_pmetrics_solver;

  // Construct dummy matrix for last argument of pred
  Matrix<double, Dynamic, Dynamic> dummy_system;
  vector<Matrix<double, Dynamic, Dynamic> >
    dummy_systems(1, dummy_system);

  Matrix <typename boost::math::tools::promote_args<T0, T1, T2, T3,
    typename boost::math::tools::promote_args<T4, T5, T6>::type>::type,
    Dynamic, Dynamic> pred;
  pred = Pred(time, amt, rate, ii, evid, cmt, addl, ss,
              pMatrix, biovar, tlag, model, f, dummy_systems);

  return pred;
}

/**
 * Overload function to allow user to pass an std::vector for 
 * pMatrix.
 */
template <typename T0, typename T1, typename T2, typename T3, typename T4,
          typename T5, typename T6, typename F>
Eigen::Matrix <typename boost::math::tools::promote_args<T0, T1, T2, T3,
  typename boost::math::tools::promote_args<T4, T5, T6>::type>::type,
  Eigen::Dynamic, Eigen::Dynamic>
generalOdeModel_bdf(const F& f,
                     const int nCmt,
                     const std::vector<T0>& time,
                     const std::vector<T1>& amt,
                     const std::vector<T2>& rate,
                     const std::vector<T3>& ii,
                     const std::vector<int>& evid,
                     const std::vector<int>& cmt,
                     const std::vector<int>& addl,
                     const std::vector<int>& ss,
                     const std::vector<T4>& pMatrix,
                     const std::vector<std::vector<T5> >& biovar,
                     const std::vector<std::vector<T6> >& tlag,
                     double rel_tol = 1e-10,
                     double abs_tol = 1e-10,
                     long int max_num_steps = 1e8) {  // NOLINT(runtime/int)
  std::vector<std::vector<T4> > vec_pMatrix(1, pMatrix);

  return generalOdeModel_bdf(f, nCmt,
                              time, amt, rate, ii, evid, cmt, addl, ss,
                              vec_pMatrix, biovar, tlag);
}

/**
* Overload function to allow user to pass an std::vector for 
* pMatrix and biovar.
*/
template <typename T0, typename T1, typename T2, typename T3, typename T4,
          typename T5, typename T6, typename F>
Eigen::Matrix <typename boost::math::tools::promote_args<T0, T1, T2, T3,
  typename boost::math::tools::promote_args<T4, T5, T6>::type>::type,
  Eigen::Dynamic, Eigen::Dynamic>
generalOdeModel_bdf(const F& f,
                     const int nCmt,
                     const std::vector<T0>& time,
                     const std::vector<T1>& amt,
                     const std::vector<T2>& rate,
                     const std::vector<T3>& ii,
                     const std::vector<int>& evid,
                     const std::vector<int>& cmt,
                     const std::vector<int>& addl,
                     const std::vector<int>& ss,
                     const std::vector<T4>& pMatrix,
                     const std::vector<T5>& biovar,
                     const std::vector<std::vector<T6> >& tlag,
                     double rel_tol = 1e-10,
                     double abs_tol = 1e-10,
                     long int max_num_steps = 1e8) {  // NOLINT(runtime/int)
  std::vector<std::vector<T4> > vec_pMatrix(1, pMatrix);
  std::vector<std::vector<T5> > vec_biovar(1, biovar);

  return generalOdeModel_bdf(f, nCmt,
                              time, amt, rate, ii, evid, cmt, addl, ss,
                              vec_pMatrix, vec_biovar, tlag);
}

/**
* Overload function to allow user to pass an std::vector for 
* pMatrix, biovar, and tlag.
*/
template <typename T0, typename T1, typename T2, typename T3, typename T4,
          typename T5, typename T6, typename F>
Eigen::Matrix <typename boost::math::tools::promote_args<T0, T1, T2, T3,
  typename boost::math::tools::promote_args<T4, T5, T6>::type>::type,
  Eigen::Dynamic, Eigen::Dynamic>
generalOdeModel_bdf(const F& f,
                     const int nCmt,
                     const std::vector<T0>& time,
                     const std::vector<T1>& amt,
                     const std::vector<T2>& rate,
                     const std::vector<T3>& ii,
                     const std::vector<int>& evid,
                     const std::vector<int>& cmt,
                     const std::vector<int>& addl,
                     const std::vector<int>& ss,
                     const std::vector<T4>& pMatrix,
                     const std::vector<T5>& biovar,
                     const std::vector<T6>& tlag,
                     double rel_tol = 1e-10,
                     double abs_tol = 1e-10,
                     long int max_num_steps = 1e8) {  // NOLINT(runtime/int)
  std::vector<std::vector<T4> > vec_pMatrix(1, pMatrix);
  std::vector<std::vector<T5> > vec_biovar(1, biovar);
  std::vector<std::vector<T5> > vec_tlag(1, tlag);

  return generalOdeModel_bdf(f, nCmt,
                              time, amt, rate, ii, evid, cmt, addl, ss,
                              vec_pMatrix, vec_biovar, vec_tlag);
}

/**
* Overload function to allow user to pass an std::vector for 
* pMatrix and tlag.
*/
template <typename T0, typename T1, typename T2, typename T3, typename T4,
          typename T5, typename T6, typename F>
Eigen::Matrix <typename boost::math::tools::promote_args<T0, T1, T2, T3,
  typename boost::math::tools::promote_args<T4, T5, T6>::type>::type,
  Eigen::Dynamic, Eigen::Dynamic>
generalOdeModel_bdf(const F& f,
                     const int nCmt,
                     const std::vector<T0>& time,
                     const std::vector<T1>& amt,
                     const std::vector<T2>& rate,
                     const std::vector<T3>& ii,
                     const std::vector<int>& evid,
                     const std::vector<int>& cmt,
                     const std::vector<int>& addl,
                     const std::vector<int>& ss,
                     const std::vector<T4>& pMatrix,
                     const std::vector<std::vector<T5> >& biovar,
                     const std::vector<T6>& tlag,
                     double rel_tol = 1e-10,
                     double abs_tol = 1e-10,
                     long int max_num_steps = 1e8) {  // NOLINT(runtime/int)
  std::vector<std::vector<T4> > vec_pMatrix(1, pMatrix);
  std::vector<std::vector<T6> > vec_tlag(1, tlag);

  return generalOdeModel_bdf(f, nCmt,
                              time, amt, rate, ii, evid, cmt, addl, ss,
                              vec_pMatrix, biovar, vec_tlag);
}

/**
* Overload function to allow user to pass an std::vector for 
* biovar.
*/
template <typename T0, typename T1, typename T2, typename T3, typename T4,
          typename T5, typename T6, typename F>
Eigen::Matrix <typename boost::math::tools::promote_args<T0, T1, T2, T3,
  typename boost::math::tools::promote_args<T4, T5, T6>::type>::type,
  Eigen::Dynamic, Eigen::Dynamic>
generalOdeModel_bdf(const F& f,
                     const int nCmt,
                     const std::vector<T0>& time,
                     const std::vector<T1>& amt,
                     const std::vector<T2>& rate,
                     const std::vector<T3>& ii,
                     const std::vector<int>& evid,
                     const std::vector<int>& cmt,
                     const std::vector<int>& addl,
                     const std::vector<int>& ss,
                     const std::vector<std::vector<T4> >& pMatrix,
                     const std::vector<T5>& biovar,
                     const std::vector<std::vector<T6> >& tlag,
                     double rel_tol = 1e-10,
                     double abs_tol = 1e-10,
                     long int max_num_steps = 1e8) {  // NOLINT(runtime/int)
  std::vector<std::vector<T5> > vec_biovar(1, biovar);

  return generalOdeModel_bdf(f, nCmt,
                              time, amt, rate, ii, evid, cmt, addl, ss,
                              pMatrix, vec_biovar, tlag);
}

/**
* Overload function to allow user to pass an std::vector for 
* biovar and tlag.
*/
template <typename T0, typename T1, typename T2, typename T3, typename T4,
          typename T5, typename T6, typename F>
Eigen::Matrix <typename boost::math::tools::promote_args<T0, T1, T2, T3,
  typename boost::math::tools::promote_args<T4, T5, T6>::type>::type,
  Eigen::Dynamic, Eigen::Dynamic>
generalOdeModel_bdf(const F& f,
                     const int nCmt,
                     const std::vector<T0>& time,
                     const std::vector<T1>& amt,
                     const std::vector<T2>& rate,
                     const std::vector<T3>& ii,
                     const std::vector<int>& evid,
                     const std::vector<int>& cmt,
                     const std::vector<int>& addl,
                     const std::vector<int>& ss,
                     const std::vector<std::vector<T4> >& pMatrix,
                     const std::vector<T5>& biovar,
                     const std::vector<T6>& tlag,
                     double rel_tol = 1e-10,
                     double abs_tol = 1e-10,
                     long int max_num_steps = 1e8) {  // NOLINT(runtime/int)
  std::vector<std::vector<T5> > vec_biovar(1, biovar);
  std::vector<std::vector<T5> > vec_tlag(1, tlag);

  return generalOdeModel_bdf(f, nCmt,
                              time, amt, rate, ii, evid, cmt, addl, ss,
                              pMatrix, vec_biovar, vec_tlag);
}

/**
 * Overload function to allow user to pass an std::vector for 
 * tlag.
 */
template <typename T0, typename T1, typename T2, typename T3, typename T4,
          typename T5, typename T6, typename F>
Eigen::Matrix <typename boost::math::tools::promote_args<T0, T1, T2, T3,
  typename boost::math::tools::promote_args<T4, T5, T6>::type>::type,
  Eigen::Dynamic, Eigen::Dynamic>
generalOdeModel_bdf(const F& f,
                     const int nCmt,
                     const std::vector<T0>& time,
                     const std::vector<T1>& amt,
                     const std::vector<T2>& rate,
                     const std::vector<T3>& ii,
                     const std::vector<int>& evid,
                     const std::vector<int>& cmt,
                     const std::vector<int>& addl,
                     const std::vector<int>& ss,
                     const std::vector<std::vector<T4> >& pMatrix,
                     const std::vector<std::vector<T5> >& biovar,
                     const std::vector<T6>& tlag,
                     double rel_tol = 1e-10,
                     double abs_tol = 1e-10,
                     long int max_num_steps = 1e8) {  // NOLINT(runtime/int)
  std::vector<std::vector<T5> > vec_tlag(1, tlag);

  return generalOdeModel_bdf(f, nCmt,
                              time, amt, rate, ii, evid, cmt, addl, ss,
                              pMatrix, biovar, vec_tlag);
}

#endif
