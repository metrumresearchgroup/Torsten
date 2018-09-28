#ifndef STAN_MATH_TORSTEN_REFACTOR_GENERALODEMODEL_RK45_HPP
#define STAN_MATH_TORSTEN_REFACTOR_GENERALODEMODEL_RK45_HPP

#include <Eigen/Dense>
#include <stan/math/torsten/PKModel/functors/general_functor.hpp>
#include <stan/math/torsten/Pred2.hpp>
#include <stan/math/torsten/pk_ode_model.hpp>
#include <stan/math/torsten/PKModel/PKModel.hpp>
#include <stan/math/torsten/PKModel/Pred/Pred1_general.hpp>
#include <stan/math/torsten/PKModel/Pred/PredSS_general.hpp>
#include <boost/math/tools/promotion.hpp>
#include <vector>

namespace torsten {

/**
 * Computes the predicted amounts in each compartment at each event
 * for a general compartment model, defined by a system of ordinary
 * differential equations. Uses the stan::math::integrate_ode_rk45 
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
 *
 * FIX ME: currently have a dummy msgs argument. Makes it easier
 * to expose to stan grammar files, because I can follow more closely
 * what was done for the ODE integrator. Not ideal.
 */
template <typename T0, typename T1, typename T2, typename T3, typename T4,
  typename T5, typename T6, typename F>
Eigen::Matrix <typename boost::math::tools::promote_args<T0, T1, T2, T3,
  typename boost::math::tools::promote_args<T4, T5, T6>::type>::type,
  Eigen::Dynamic, Eigen::Dynamic>
generalOdeModel_rk45(const F& f,
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
                     std::ostream* msgs = 0,
                     double rel_tol = 1e-6,
                     double abs_tol = 1e-6,
                     long int max_num_steps = 1e6) {  // NOLINT(runtime/int)
  using std::vector;
  using Eigen::Dynamic;
  using Eigen::Matrix;
  using boost::math::tools::promote_args;

  // check arguments
  static const char* function("generalOdeModel_rk45");
  torsten::pmetricsCheck(time, amt, rate, ii, evid, cmt, addl, ss,
    pMatrix, biovar, tlag, function);

  // Construct dummy matrix for last argument of pred
  Matrix<T4, Dynamic, Dynamic> dummy_system;
  vector<Matrix<T4, Dynamic, Dynamic> >
    dummy_systems(1, dummy_system);

  typedef general_functor<F> F0;

  PredWrapper<refactor::PKODEModel, StanRk45> pr(rel_tol, abs_tol, max_num_steps, msgs);

  const Pred1_general<F0> pred1(F0(f), rel_tol, abs_tol,
                                max_num_steps, msgs, "rk45");
  const PredSS_general<F0> predss (F0(f), rel_tol, abs_tol,
                                   max_num_steps, msgs, "rk45", nCmt);

#ifdef OLD_TORSTEN
  return Pred(time, amt, rate, ii, evid, cmt, addl, ss,
              pMatrix, biovar, tlag, nCmt, dummy_systems,
              pred1, predss);

#else
  return pr.Pred2(time, amt, rate, ii, evid, cmt, addl, ss,
                  pMatrix, biovar, tlag, nCmt, dummy_systems,
                  pred1, predss,
                  f, nCmt);
#endif
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
generalOdeModel_rk45(const F& f,
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
                     std::ostream* msgs = 0,
                     double rel_tol = 1e-6,
                     double abs_tol = 1e-6,
                     long int max_num_steps = 1e6) {  // NOLINT(runtime/int)
  std::vector<std::vector<T4> > vec_pMatrix(1, pMatrix);

  return generalOdeModel_rk45(f, nCmt,
                              time, amt, rate, ii, evid, cmt, addl, ss,
                              vec_pMatrix, biovar, tlag,
                              msgs, rel_tol, abs_tol, max_num_steps);
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
generalOdeModel_rk45(const F& f,
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
                     std::ostream* msgs = 0,
                     double rel_tol = 1e-6,
                     double abs_tol = 1e-6,
                     long int max_num_steps = 1e6) {  // NOLINT(runtime/int)
  std::vector<std::vector<T4> > vec_pMatrix(1, pMatrix);
  std::vector<std::vector<T5> > vec_biovar(1, biovar);

  return generalOdeModel_rk45(f, nCmt,
                              time, amt, rate, ii, evid, cmt, addl, ss,
                              vec_pMatrix, vec_biovar, tlag,
                              msgs, rel_tol, abs_tol, max_num_steps);
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
generalOdeModel_rk45(const F& f,
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
                     std::ostream* msgs = 0,
                     double rel_tol = 1e-6,
                     double abs_tol = 1e-6,
                     long int max_num_steps = 1e6) {  // NOLINT(runtime/int)
  std::vector<std::vector<T4> > vec_pMatrix(1, pMatrix);
  std::vector<std::vector<T5> > vec_biovar(1, biovar);
  std::vector<std::vector<T6> > vec_tlag(1, tlag);

  return generalOdeModel_rk45(f, nCmt,
                              time, amt, rate, ii, evid, cmt, addl, ss,
                              vec_pMatrix, vec_biovar, vec_tlag,
                              msgs, rel_tol, abs_tol, max_num_steps);
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
generalOdeModel_rk45(const F& f,
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
                     std::ostream* msgs = 0,
                     double rel_tol = 1e-6,
                     double abs_tol = 1e-6,
                     long int max_num_steps = 1e6) {  // NOLINT(runtime/int)
  std::vector<std::vector<T4> > vec_pMatrix(1, pMatrix);
  std::vector<std::vector<T6> > vec_tlag(1, tlag);

  return generalOdeModel_rk45(f, nCmt,
                              time, amt, rate, ii, evid, cmt, addl, ss,
                              vec_pMatrix, biovar, vec_tlag,
                              msgs, rel_tol, abs_tol, max_num_steps);
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
generalOdeModel_rk45(const F& f,
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
                     std::ostream* msgs = 0,
                     double rel_tol = 1e-6,
                     double abs_tol = 1e-6,
                     long int max_num_steps = 1e6) {  // NOLINT(runtime/int)
  std::vector<std::vector<T5> > vec_biovar(1, biovar);

  return generalOdeModel_rk45(f, nCmt,
                              time, amt, rate, ii, evid, cmt, addl, ss,
                              pMatrix, vec_biovar, tlag,
                              msgs, rel_tol, abs_tol, max_num_steps);
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
generalOdeModel_rk45(const F& f,
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
                     std::ostream* msgs = 0,
                     double rel_tol = 1e-6,
                     double abs_tol = 1e-6,
                     long int max_num_steps = 1e6) {  // NOLINT(runtime/int)
  std::vector<std::vector<T5> > vec_biovar(1, biovar);
  std::vector<std::vector<T6> > vec_tlag(1, tlag);

  return generalOdeModel_rk45(f, nCmt,
                              time, amt, rate, ii, evid, cmt, addl, ss,
                              pMatrix, vec_biovar, vec_tlag,
                              msgs, rel_tol, abs_tol, max_num_steps);
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
generalOdeModel_rk45(const F& f,
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
                     std::ostream* msgs = 0,
                     double rel_tol = 1e-6,
                     double abs_tol = 1e-6,
                     long int max_num_steps = 1e6) {  // NOLINT(runtime/int)
  std::vector<std::vector<T6> > vec_tlag(1, tlag);

  return generalOdeModel_rk45(f, nCmt,
                              time, amt, rate, ii, evid, cmt, addl, ss,
                              pMatrix, biovar, vec_tlag,
                              msgs, rel_tol, abs_tol, max_num_steps);
}

}
#endif
