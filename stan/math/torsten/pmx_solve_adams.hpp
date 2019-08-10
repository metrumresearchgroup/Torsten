#ifndef STAN_MATH_TORSTEN_REFACTOR_GENERALODEMODEL_ADAMS_HPP
#define STAN_MATH_TORSTEN_REFACTOR_GENERALODEMODEL_ADAMS_HPP

#include <Eigen/Dense>
#include <stan/math/torsten/to_array_2d.hpp>
#include <stan/math/torsten/events_manager.hpp>
#include <stan/math/torsten/pmx_population_check.hpp>
#include <stan/math/torsten/PKModel/functors/general_functor.hpp>
#include <stan/math/torsten/PKModel/PKModel.hpp>
#include <stan/math/torsten/PKModel/Pred/Pred1_general.hpp>
#include <stan/math/torsten/PKModel/Pred/PredSS_general.hpp>
#include <stan/math/torsten/pmx_ode_model.hpp>
#include <stan/math/torsten/event_solver.hpp>
#include <boost/math/tools/promotion.hpp>
#include <vector>

namespace torsten {

/**
 * Computes the predicted amounts in each compartment at each event
 * for a general compartment model, defined by a system of ordinary
 * differential equations. Uses the stan::math::integrate_ode_adams 
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
Eigen::Matrix <typename torsten::return_t<T0, T1, T2, T3, T4, T5, T6>::type,
               Eigen::Dynamic, Eigen::Dynamic>
pmx_solve_adams(const F& f,
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
                    double rel_tol = 1e-10,
                    double abs_tol = 1e-10,
                    long int max_num_steps = 1e8) {  // NOLINT(runtime/int)
  using std::vector;
  using Eigen::Dynamic;
  using Eigen::Matrix;
  using boost::math::tools::promote_args;
  using refactor::PKRec;

  // check arguments
  static const char* function("pmx_solve_adams");
  torsten::pmx_check(time, amt, rate, ii, evid, cmt, addl, ss,
                pMatrix, biovar, tlag, function);

  // Construct dummy matrix for last argument of pred
  Matrix<T4, Dynamic, Dynamic> dummy_system;
  vector<Matrix<T4, Dynamic, Dynamic> >
    dummy_systems(1, dummy_system);

  typedef general_functor<F> F0;

  const Pred1_general<F0> pred1(F0(f), rel_tol, abs_tol,
                                max_num_steps, msgs, "adams");
  const PredSS_general<F0> predss (F0(f), rel_tol, abs_tol,
                                   max_num_steps, msgs, "adams", nCmt);

#ifdef OLD_TORSTEN
  return Pred(time, amt, rate, ii, evid, cmt, addl, ss,
              pMatrix, biovar, tlag, nCmt, dummy_systems,
              pred1, predss);
#else
  using ER = NONMENEventsRecord<T0, T1, T2, T3, std::vector<T4>, T5, T6>;
  using EM = EventsManager<ER>;
  const ER events_rec(nCmt, time, amt, rate, ii, evid, cmt, addl, ss, pMatrix, biovar, tlag);

  Matrix<typename EM::T_scalar, Dynamic, Dynamic> pred =
    Matrix<typename EM::T_scalar, Dynamic, Dynamic>::Zero(events_rec.num_event_times(), EM::nCmt(events_rec));

  using model_type = refactor::PKODEModel<typename EM::T_time, typename EM::T_scalar, typename EM::T_rate, typename EM::T_par, F>;

#ifdef TORSTEN_USE_STAN_ODE
  PMXOdeIntegrator<StanAdams> integrator(rel_tol, abs_tol, max_num_steps, msgs);
  EventSolver<model_type, PMXOdeIntegrator<StanAdams>&> pr;
#else
  PMXOdeIntegrator<PkAdams> integrator(rel_tol, abs_tol, max_num_steps, msgs);
  EventSolver<model_type, PMXOdeIntegrator<PkAdams>&> pr;
#endif

  pr.pred(0, events_rec, pred, integrator, f);
  return pred;

#endif

}

/**
 * Overload function to allow user to pass an std::vector for 
 * pMatrix/bioavailability/tlag
 */
template <typename T0, typename T1, typename T2, typename T3,
          typename T_par, typename T_biovar, typename T_tlag, typename F,
          typename std::enable_if_t<!(torsten::is_std_vector<T_par, T_biovar, T_tlag>::value)>* = nullptr> //NOLINT
Eigen::Matrix <typename torsten::return_t<T0, T1, T2, T3,
                                          typename torsten::value_type<T_par>::type,
                                          typename torsten::value_type<T_biovar>::type,
                                          typename torsten::value_type<T_tlag>::type>::type,
               Eigen::Dynamic, Eigen::Dynamic>
pmx_solve_adams(const F& f,
                const int nCmt,
                const std::vector<T0>& time,
                const std::vector<T1>& amt,
                const std::vector<T2>& rate,
                const std::vector<T3>& ii,
                const std::vector<int>& evid,
                const std::vector<int>& cmt,
                const std::vector<int>& addl,
                const std::vector<int>& ss,
                const std::vector<T_par>& pMatrix,
                const std::vector<T_biovar>& biovar,
                const std::vector<T_tlag>& tlag,
                std::ostream* msgs = 0,
                double rel_tol = 1e-6,
                double abs_tol = 1e-6,
                long int max_num_steps = 1e6) {
  auto param_ = torsten::to_array_2d(pMatrix);
  auto biovar_ = torsten::to_array_2d(biovar);
  auto tlag_ = torsten::to_array_2d(tlag);

  return pmx_solve_adams(f, nCmt,
                         time, amt, rate, ii, evid, cmt, addl, ss,
                         param_, biovar_, tlag_,
                         msgs, rel_tol, abs_tol, max_num_steps);
}

  /*
   * For backward compatibility we keep old version of
   * return type using transpose. This is less efficient and
   * will be decomissioned in formal release.
   */
  template <typename T0, typename T1, typename T2, typename T3,
            typename T_par, typename T_biovar, typename T_tlag,
            typename F>
  Eigen::Matrix <typename torsten::return_t<T0, T1, T2, T3,
                                            typename torsten::value_type<T_par>::type,
                                            typename torsten::value_type<T_biovar>::type,
                                            typename torsten::value_type<T_tlag>::type>::type,
                 Eigen::Dynamic, Eigen::Dynamic>
  generalOdeModel_adams(const F& f,
                      const int nCmt,
                      const std::vector<T0>& time,
                      const std::vector<T1>& amt,
                      const std::vector<T2>& rate,
                      const std::vector<T3>& ii,
                      const std::vector<int>& evid,
                      const std::vector<int>& cmt,
                      const std::vector<int>& addl,
                      const std::vector<int>& ss,
                      const std::vector<T_par>& pMatrix,
                      const std::vector<T_biovar>& biovar,
                      const std::vector<T_tlag>& tlag,
                      std::ostream* msgs = 0,
                      double rel_tol = 1e-6,
                      double abs_tol = 1e-6,
                      long int max_num_steps = 1e6) {
    auto x = pmx_solve_adams(f, nCmt,
                           time, amt, rate, ii, evid, cmt, addl, ss,
                           pMatrix, biovar, tlag,
                           msgs, rel_tol, abs_tol, max_num_steps);
    return x.transpose();
  }

  /*
   * For population models, more often we use ragged arrays
   * to describe the entire population, so in addition we need the arrays of
   * the length of each individual's data. The size of that
   * vector is the size of
   * the population.
   */
template <typename T0, typename T1, typename T2, typename T3, typename T4,
          typename T5, typename T6, typename F>
Eigen::Matrix<typename EventsManager<NONMENEventsRecord<T0, T1, T2, T3, std::vector<T4>, T5, T6> >::T_scalar, // NOLINT
              Eigen::Dynamic, Eigen::Dynamic>
pmx_solve_group_adams(const F& f,
                      const int nCmt,
                      const std::vector<int>& len,
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
                      double rel_tol = 1e-10,
                      double abs_tol = 1e-10,
                      long int max_num_steps = 1e8) {
  static const char* caller("pmx_solve_group_adams");
  torsten::pmx_population_check(len, time, amt, rate, ii, evid, cmt, addl, ss,
                                pMatrix, biovar, tlag, caller);

  using ER = NONMENEventsRecord<T0, T1, T2, T3, std::vector<T4>, T5, T6>;
  using EM = EventsManager<ER>;
  using model_type = refactor::PKODEModel<typename EM::T_time, typename EM::T_scalar, typename EM::T_rate, typename EM::T_par, F>;

  ER events_rec(nCmt, len, time, amt, rate, ii, evid, cmt, addl, ss, pMatrix, biovar, tlag);

#ifdef TORSTEN_USE_STAN_ODE
  PMXOdeIntegrator<StanAdams> integrator(rel_tol, abs_tol, max_num_steps, msgs);
  EventSolver<model_type, PMXOdeIntegrator<StanAdams>&> pr;
#else
  PMXOdeIntegrator<PkAdams> integrator(rel_tol, abs_tol, max_num_steps, msgs);
  EventSolver<model_type, PMXOdeIntegrator<PkAdams>&> pr;
#endif

  Eigen::Matrix<typename EM::T_scalar, -1, -1> pred(nCmt, events_rec.total_num_event_times);

  pr.pred(events_rec, pred, integrator, f);

  return pred;
}

}

#endif
