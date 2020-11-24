#ifndef STAN_MATH_TORSTEN_TWOCPT_HPP
#define STAN_MATH_TORSTEN_TWOCPT_HPP

#include <Eigen/Dense>
#include <stan/math/prim/err/check_greater_or_equal.hpp>
#include <stan/math/torsten/to_array_2d.hpp>
#include <stan/math/torsten/is_std_vector.hpp>
#include <stan/math/torsten/pmx_solve_ode.hpp>
#include <stan/math/torsten/pmx_solve_cpt.hpp>
#include <stan/math/torsten/ev_solver.hpp>
#include <stan/math/torsten/ev_manager.hpp>
#include <stan/math/torsten/pmx_twocpt_model.hpp>
#include <stan/math/torsten/pmx_population_check.hpp>
#include <stan/math/prim/meta/return_type.hpp>
#include <string>
#include <vector>

namespace torsten {

/**
 * Computes the predicted amounts in each compartment at each event
 * for a two-compartment model with analytical solution.
 *
 * @tparam Ts types of parameters, see <code>pmx_solve_cpt</code> for
 *         details. 
 * @return a matrix with predicted amount in each compartment 
 *         at each event. 
 *
 */
  template <typename... Ts>
  auto pmx_solve_twocpt(Ts... args) {
    return PMXSolveCPT<PMXTwoCptModel>::solve(args...);
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
  Eigen::Matrix <typename stan::return_type_t<T0, T1, T2, T3,
                                            typename torsten::value_type<T_par>::type,
                                            typename torsten::value_type<T_biovar>::type,
                                            typename torsten::value_type<T_tlag>::type>,
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
  Eigen::Matrix <typename stan::return_type_t<T0, T1, T2, T3, T4, T5, T6>,
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
  Eigen::Matrix <typename stan::return_type_t<T0, T1, T2, T3,
                                            typename torsten::value_type<T_par>::type,
                                            typename torsten::value_type<T_biovar>::type,
                                            typename torsten::value_type<T_tlag>::type>,
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
Eigen::Matrix<typename  EventsManager<NONMENEventsRecord<T0, T1, T2, T3>, NonEventParameters<T0, T4, std::vector, std::tuple<T5, T6> > >::T_scalar, // NOLINT
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
  using ER = NONMENEventsRecord<T0, T1, T2, T3>;
  using EM = EventsManager<ER, NonEventParameters<T0, T4, std::vector, std::tuple<T5, T6> >>;

  int nCmt = torsten::PMXTwoCptModel<double>::Ncmt;
  ER events_rec(nCmt, len, time, amt, rate, ii, evid, cmt, addl, ss);

  static const char* caller("pmx_solve_group_twocpt");
  torsten::pmx_population_check(len, time, amt, rate, ii, evid, cmt, addl, ss,
                                pMatrix, biovar, tlag, caller);

  using model_type = torsten::PMXTwoCptModel<typename EM::T_par>;
   EventSolver<model_type, EM> pr;

  Eigen::Matrix<typename EM::T_scalar, -1, -1> pred(events_rec.total_num_event_times, nCmt);

  pr.pred(events_rec, pred, PMXOdeIntegrator<Analytical>(), pMatrix, biovar, tlag);

  return pred;
}

}
#endif
