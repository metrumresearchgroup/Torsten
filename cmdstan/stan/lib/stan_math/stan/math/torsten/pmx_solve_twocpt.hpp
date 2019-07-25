#ifndef STAN_MATH_TORSTEN_PKMODELTWOCPT2_HPP
#define STAN_MATH_TORSTEN_PKMODELTWOCPT2_HPP

#include <Eigen/Dense>
#include <stan/math/prim/scal/err/check_greater_or_equal.hpp>
#include <boost/math/tools/promotion.hpp>
#include <stan/math/torsten/event_solver.hpp>
#include <stan/math/torsten/events_manager.hpp>
#include <stan/math/torsten/PKModel/PKModel.hpp>
#include <stan/math/torsten/PKModel/Pred/Pred1_twoCpt.hpp>
#include <stan/math/torsten/PKModel/Pred/PredSS_twoCpt.hpp>
#include <stan/math/torsten/pmx_twocpt_model.hpp>
#include <stan/math/torsten/pmx_population_check.hpp>
#include <stan/math/torsten/return_type.hpp>
#include <string>
#include <vector>

namespace torsten {

/**
 * Computes the predicted amounts in each compartment at each event
 * for a two compartment model with first oder absorption.
 *
 * @tparam T0 type of scalar for time of events.
 * @tparam T1 type of scalar for amount at each event.
 * @tparam T2 type of scalar for rate at each event.
 * @tparam T3 type of scalar for inter-dose inteveral at each event.
 * @tparam T4 type of scalars for the model parameters.
 * @tparam T5 type of scalars for the 
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
template <typename T0, typename T1, typename T2, typename T3, typename T4,
          typename T5, typename T6>
Eigen::Matrix<typename torsten::return_t<T0, T1, T2, T3, T4, T5, T6>::type,
              Eigen::Dynamic, Eigen::Dynamic>
pmx_solve_twocpt(const std::vector<T0>& time,
              const std::vector<T1>& amt,
              const std::vector<T2>& rate,
              const std::vector<T3>& ii,
              const std::vector<int>& evid,
              const std::vector<int>& cmt,
              const std::vector<int>& addl,
              const std::vector<int>& ss,
              const std::vector<std::vector<T4> >& pMatrix,
              const std::vector<std::vector<T5> >& biovar,
              const std::vector<std::vector<T6> >& tlag) {
  using std::vector;
  using Eigen::Dynamic;
  using Eigen::Matrix;
  using boost::math::tools::promote_args;
  using stan::math::check_positive_finite;
  using refactor::PKRec;

  int nCmt = refactor::PMXTwoCptModel<double, double, double, double>::Ncmt;
  int nParms = refactor::PMXTwoCptModel<double, double, double, double>::Npar;
  static const char* function("pmx_solve_twocpt");

  // Check arguments
  torsten::pmx_check(time, amt, rate, ii, evid, cmt, addl, ss,
                pMatrix, biovar, tlag, function);
  for (size_t i = 0; i < pMatrix.size(); i++) {
    check_positive_finite(function, "PK parameter CL", pMatrix[i][0]);
    check_positive_finite(function, "PK parameter Q", pMatrix[i][1]);
    check_positive_finite(function, "PK parameter V2", pMatrix[i][2]);
    check_positive_finite(function, "PK parameter V3", pMatrix[i][3]);
  }
  std::string message4 = ", but must equal the number of parameters in the model: " // NOLINT
    + boost::lexical_cast<std::string>(nParms) + "!";
  const char* length_error4 = message4.c_str();
  if (!(pMatrix[0].size() == (size_t) nParms))
    stan::math::invalid_argument(function,
    "The number of parameters per event (length of a vector in the first argument) is", // NOLINT
    pMatrix[0].size(), "", length_error4);

  // FIX ME - we want to check every array of pMatrix, not
  // just the first one (at index 0)
  std::string message5 = ", but must equal the number of parameters in the model: " // NOLINT
  + boost::lexical_cast<std::string>(nParms) + "!";
  const char* length_error5 = message5.c_str();
  if (!(pMatrix[0].size() == (size_t) nParms))
    stan::math::invalid_argument(function,
    "The number of parameters per event (length of a vector in the ninth argument) is", // NOLINT
    pMatrix[0].size(), "", length_error5);

  std::string message6 = ", but must equal the number of compartments in the model: " // NOLINT
  + boost::lexical_cast<std::string>(nCmt) + "!";
  const char* length_error6 = message6.c_str();
  if (!(biovar[0].size() == (size_t) nCmt))
    stan::math::invalid_argument(function,
    "The number of biovariability parameters per event (length of a vector in the tenth argument) is", // NOLINT
    biovar[0].size(), "", length_error6);

  if (!(tlag[0].size() == (size_t) nCmt))
    stan::math::invalid_argument(function,
    "The number of lag times parameters per event (length of a vector in the eleventh argument) is", // NOLINT
    tlag[0].size(), "", length_error5);

  // Construct dummy matrix for last argument of pred
  Matrix<T4, Dynamic, Dynamic> dummy_system;
  vector<Matrix<T4, Dynamic, Dynamic> >
    dummy_systems(1, dummy_system);

#ifdef OLD_TORSTEN
  return Pred(time, amt, rate, ii, evid, cmt, addl, ss,
              pMatrix, biovar, tlag,
              nCmt, dummy_systems,
              Pred1_twoCpt(), PredSS_twoCpt());
#else
  using ER = NONMENEventsRecord<T0, T1, T2, T3, std::vector<T4>, T5, T6>;
  using EM = EventsManager<ER>;
  const ER events_rec(nCmt, time, amt, rate, ii, evid, cmt, addl, ss, pMatrix, biovar, tlag);

  Matrix<typename EM::T_scalar, Dynamic, Dynamic> pred =
    Matrix<typename EM::T_scalar, Dynamic, Dynamic>::Zero(events_rec.num_event_times(), EM::nCmt(events_rec));

  using model_type = refactor::PMXTwoCptModel<typename EM::T_time, typename EM::T_scalar, typename EM::T_rate, typename EM::T_par>;
  EventSolver<model_type> pr;
  pr.pred(0, events_rec, pred);
  return pred;

#endif
}

/**
 * Overload function to allow user to pass an std::vector for 
 * pMatrix/bioavailability/tlag
 */
  template <typename T0, typename T1, typename T2, typename T3,
            typename T_par, typename T_biovar, typename T_tlag,
               typename
            std::enable_if_t<
              !(torsten::is_std_vector<T_par>::value && torsten::is_std_vector<T_biovar>::value && torsten::is_std_vector<T_tlag>::value)>* = nullptr> //NOLINT
  Eigen::Matrix <typename torsten::return_t<T0, T1, T2, T3,
                                            typename torsten::value_type<T_par>::type,
                                            typename torsten::value_type<T_biovar>::type,
                                            typename torsten::value_type<T_tlag>::type>::type,
                 Eigen::Dynamic, Eigen::Dynamic>
  pmx_solve_twocpt(const std::vector<T0>& time,
                   const std::vector<T1>& amt,
                   const std::vector<T2>& rate,
                   const std::vector<T3>& ii,
                   const std::vector<int>& evid,
                   const std::vector<int>& cmt,
                   const std::vector<int>& addl,
                   const std::vector<int>& ss,
                   const std::vector<T_par>& pMatrix,
                   const std::vector<T_biovar>& biovar,
                   const std::vector<T_tlag>& tlag) {
    auto param_ = torsten::to_array_2d(pMatrix);
    auto biovar_ = torsten::to_array_2d(biovar);
    auto tlag_ = torsten::to_array_2d(tlag);

    return pmx_solve_twocpt(time, amt, rate, ii, evid, cmt, addl, ss,
                            param_, biovar_, tlag_);
  }

  // old version
  template <typename T0, typename T1, typename T2, typename T3, typename T4,
            typename T5, typename T6>
  Eigen::Matrix <typename torsten::return_t<T0, T1, T2, T3, T4, T5, T6>::type,
                 Eigen::Dynamic, Eigen::Dynamic>
  PKModelTwoCpt(const std::vector<T0>& time,
                const std::vector<T1>& amt,
                const std::vector<T2>& rate,
                const std::vector<T3>& ii,
                const std::vector<int>& evid,
                const std::vector<int>& cmt,
                const std::vector<int>& addl,
                const std::vector<int>& ss,
                const std::vector<std::vector<T4> >& pMatrix,
                const std::vector<std::vector<T5> >& biovar,
                const std::vector<std::vector<T6> >& tlag) {
    auto x = pmx_solve_twocpt(time, amt, rate, ii, evid, cmt, addl, ss,
                              pMatrix, biovar, tlag);
    return x.transpose();
  }

  template <typename T0, typename T1, typename T2, typename T3,
            typename T_par, typename T_biovar, typename T_tlag,
               typename
            std::enable_if_t<
              !(torsten::is_std_vector<T_par>::value && torsten::is_std_vector<T_biovar>::value && torsten::is_std_vector<T_tlag>::value)>* = nullptr> //NOLINT
  Eigen::Matrix <typename torsten::return_t<T0, T1, T2, T3,
                                            typename torsten::value_type<T_par>::type,
                                            typename torsten::value_type<T_biovar>::type,
                                            typename torsten::value_type<T_tlag>::type>::type,
                 Eigen::Dynamic, Eigen::Dynamic>
  PKModelTwoCpt(const std::vector<T0>& time,
                   const std::vector<T1>& amt,
                   const std::vector<T2>& rate,
                   const std::vector<T3>& ii,
                   const std::vector<int>& evid,
                   const std::vector<int>& cmt,
                   const std::vector<int>& addl,
                   const std::vector<int>& ss,
                   const std::vector<T_par>& pMatrix,
                   const std::vector<T_biovar>& biovar,
                   const std::vector<T_tlag>& tlag) {
    auto x = pmx_solve_twocpt(time, amt, rate, ii, evid, cmt, addl, ss,
                            pMatrix, biovar, tlag);
    return x.transpose();
  }

  /* 
   * For population models, we follow the call signature
   * but add the arrays of the length of each individual's data. 
   * The size of that vector is the size of
   * the population.
   */
template <typename T0, typename T1, typename T2, typename T3, typename T4,
          typename T5, typename T6>
Eigen::Matrix<typename EventsManager<NONMENEventsRecord<T0, T1, T2, T3, std::vector<T4>, T5, T6> >::T_scalar, // NOLINT
              Eigen::Dynamic, Eigen::Dynamic>
pmx_solve_group_twocpt(const std::vector<int>& len,
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
                       const std::vector<std::vector<T6> >& tlag) {
  using ER = NONMENEventsRecord<T0, T1, T2, T3, std::vector<T4>, T5, T6>;
  using EM = EventsManager<ER>;

  int nCmt = refactor::PMXTwoCptModel<double, double, double, double>::Ncmt;
  ER events_rec(nCmt, len, time, amt, rate, ii, evid, cmt, addl, ss, pMatrix, biovar, tlag);

  static const char* caller("pmx_solve_group_twocpt");
  torsten::pmx_population_check(len, time, amt, rate, ii, evid, cmt, addl, ss,
                                pMatrix, biovar, tlag, caller);

  using model_type = refactor::PMXTwoCptModel<typename EM::T_time, typename EM::T_scalar, typename EM::T_rate, typename EM::T_par>;
  EventSolver<model_type> pr;

  Eigen::Matrix<typename EM::T_scalar, -1, -1> pred(events_rec.total_num_event_times, nCmt);

  pr.pred(events_rec, pred);

  return pred;
}

}
#endif
