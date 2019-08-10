#ifndef STAN_MATH_TORSTEN_REFACTOR_LINODEMODEL_HPP
#define STAN_MATH_TORSTEN_REFACTOR_LINODEMODEL_HPP

#include <Eigen/Dense>
#include <stan/math/torsten/events_manager.hpp>
#include <boost/math/tools/promotion.hpp>
#include <stan/math/torsten/event_solver.hpp>
#include <stan/math/torsten/pmx_linode_model.hpp>
#include <stan/math/torsten/PKModel/PKModel.hpp>
#include <stan/math/torsten/PKModel/Pred/Pred1_linOde.hpp>
#include <stan/math/torsten/PKModel/Pred/PredSS_linOde.hpp>
#include <stan/math/prim/mat/err/check_square.hpp>
#include <vector>

namespace torsten {

/**
 * Computes the predicted amounts in each compartment at each event
 * for a compartment model, described by a linear system of ordinary
 * differential equations. Uses the stan::math::matrix_exp 
 * function.
 *
 * @tparam T0 type of scalar for time of events. 
 * @tparam T1 type of scalar for amount at each event.
 * @tparam T2 type of scalar for rate at each event.
 * @tparam T3 type of scalar for inter-dose inteveral at each event.
 * @tparam T4 type of scalar for matrix describing linear ODE system.
 * @tparam T5 type of scalars for bio-variability parameters.
 * @tparam T6 type of scalars for tlag parameters 
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
 * between time-points
 * @param[in] system square matrix describing the linear system of ODEs
 * @param[in] bio-variability at each event
 * @param[in] lag times at each event
 * @return a matrix with predicted amount in each compartment 
 * at each event.
 */
template <typename T0, typename T1, typename T2, typename T3,
          typename T4, typename T5, typename T6>
Eigen::Matrix <typename torsten::return_t<T0, T1, T2, T3, T4, T5, T6>::type,
               Eigen::Dynamic, Eigen::Dynamic>
pmx_solve_linode(const std::vector<T0>& time,
            const std::vector<T1>& amt,
            const std::vector<T2>& rate,
            const std::vector<T3>& ii,
            const std::vector<int>& evid,
            const std::vector<int>& cmt,
            const std::vector<int>& addl,
            const std::vector<int>& ss,
            const std::vector< Eigen::Matrix<T4, Eigen::Dynamic, Eigen::Dynamic> >& system,
            const std::vector<std::vector<T5> >& biovar,
            const std::vector<std::vector<T6> >& tlag) {
  using std::vector;
  using Eigen::Dynamic;
  using Eigen::Matrix;
  using boost::math::tools::promote_args;
  using refactor::PKRec;

  static const char* function("pmx_solve_linode");
  for (size_t i = 0; i < system.size(); i++)
    stan::math::check_square(function, "system matrix", system[i]);
  int nCmt = system[0].cols();

  std::vector<T4> parameters_dummy(0);
  std::vector<std::vector<T4> > pMatrix_dummy(1, parameters_dummy);
  torsten::pmx_check(time, amt, rate, ii, evid, cmt, addl, ss,
                pMatrix_dummy, biovar, tlag, function);

#ifdef OLD_TORSTEN
  return Pred(time, amt, rate, ii, evid, cmt, addl, ss,
              pMatrix_dummy, biovar, tlag, nCmt, system,
              Pred1_linOde(), PredSS_linOde());
#else
  using ER = NONMENEventsRecord<T0, T1, T2, T3, Eigen::Matrix<T4,-1,-1>, T5, T6>;
  using EM = EventsManager<ER>;
  const ER events_rec(nCmt, time, amt, rate, ii, evid, cmt, addl, ss, system, biovar, tlag);

  Matrix<typename EM::T_scalar, Dynamic, Dynamic> pred =
    Matrix<typename EM::T_scalar, Dynamic, Dynamic>::Zero(events_rec.num_event_times(), EM::nCmt(events_rec));

  using model_type = refactor::PMXLinODEModel<typename EM::T_time, typename EM::T_scalar, typename EM::T_rate, typename EM::T_par>;
  EventSolver<model_type> pr;
  pr.pred(0, events_rec, pred);
  return pred;

#endif
}

/**
 * Overload function to allow user to pass a matrix for 
 * system.
 */
template <typename T0, typename T1, typename T2, typename T3,
          typename T4, typename T_biovar, typename T_tlag>
Eigen::Matrix <typename torsten::return_t<T0, T1, T2, T3, T4,
                                          typename torsten::value_type<T_biovar>::type,
                                          typename torsten::value_type<T_tlag>::type>::type,
               Eigen::Dynamic, Eigen::Dynamic>
pmx_solve_linode(const std::vector<T0>& time,
                 const std::vector<T1>& amt,
                 const std::vector<T2>& rate,
                 const std::vector<T3>& ii,
                 const std::vector<int>& evid,
                 const std::vector<int>& cmt,
                 const std::vector<int>& addl,
                 const std::vector<int>& ss,
                 const Eigen::Matrix<T4, -1, -1>& system,
                 const std::vector<T_biovar>& biovar,
                 const std::vector<T_tlag>& tlag) {
  std::vector<Eigen::Matrix<T4, -1, -1> > system_{system};
  auto biovar_ = torsten::to_array_2d(biovar);
  auto tlag_ = torsten::to_array_2d(tlag);

  return pmx_solve_linode(time, amt, rate, ii, evid, cmt, addl, ss,
                          system_, biovar_, tlag_);
}

/**
 * Overload function to allow user to pass a matrix for 
 * system.
 */
template <typename T0, typename T1, typename T2, typename T3,
          typename T4, typename T_biovar, typename T_tlag,
          typename
          std::enable_if_t<
            !(torsten::is_std_vector<T_biovar>::value && torsten::is_std_vector<T_tlag>::value)>* = nullptr> //NOLINT
Eigen::Matrix <typename torsten::return_t<T0, T1, T2, T3, T4,
                                          typename torsten::value_type<T_biovar>::type,
                                          typename torsten::value_type<T_tlag>::type>::type,
               Eigen::Dynamic, Eigen::Dynamic>
pmx_solve_linode(const std::vector<T0>& time,
                 const std::vector<T1>& amt,
                 const std::vector<T2>& rate,
                 const std::vector<T3>& ii,
                 const std::vector<int>& evid,
                 const std::vector<int>& cmt,
                 const std::vector<int>& addl,
                 const std::vector<int>& ss,
                 const std::vector< Eigen::Matrix<T4, -1, -1> >& system,
                 const std::vector<T_biovar>& biovar,
                 const std::vector<T_tlag>& tlag) {
  auto biovar_ = torsten::to_array_2d(biovar);
  auto tlag_ = torsten::to_array_2d(tlag);

  return pmx_solve_linode(time, amt, rate, ii, evid, cmt, addl, ss,
                          system, biovar_, tlag_);
}

  // old version by using transpose
template <typename T0, typename T1, typename T2, typename T3,
  typename T4, typename T5, typename T6>
Eigen::Matrix <typename boost::math::tools::promote_args<T0, T1, T2, T3,
                                                         typename boost::math::tools::promote_args<T4, T5, T6>::type>::type,
               Eigen::Dynamic, Eigen::Dynamic>
linOdeModel(const std::vector<T0>& time,
            const std::vector<T1>& amt,
            const std::vector<T2>& rate,
            const std::vector<T3>& ii,
            const std::vector<int>& evid,
            const std::vector<int>& cmt,
            const std::vector<int>& addl,
            const std::vector<int>& ss,
            const std::vector< Eigen::Matrix<T4, Eigen::Dynamic, Eigen::Dynamic> >& system,
            const std::vector<std::vector<T5> >& biovar,
            const std::vector<std::vector<T6> >& tlag) {
  auto x = pmx_solve_linode(time, amt, rate, ii, evid, cmt, addl, ss, system, biovar, tlag);
  return x.transpose();
}

template <typename T0, typename T1, typename T2, typename T3,
          typename T4, typename T_biovar, typename T_tlag>
Eigen::Matrix <typename torsten::return_t<T0, T1, T2, T3, T4,
                                          typename torsten::value_type<T_biovar>::type,
                                          typename torsten::value_type<T_tlag>::type>::type,
               Eigen::Dynamic, Eigen::Dynamic>
linOdeModel(const std::vector<T0>& time,
                 const std::vector<T1>& amt,
                 const std::vector<T2>& rate,
                 const std::vector<T3>& ii,
                 const std::vector<int>& evid,
                 const std::vector<int>& cmt,
                 const std::vector<int>& addl,
                 const std::vector<int>& ss,
                 const Eigen::Matrix<T4, -1, -1>& system,
                 const std::vector<T_biovar>& biovar,
                 const std::vector<T_tlag>& tlag) {
  auto x = pmx_solve_linode(time, amt, rate, ii, evid, cmt, addl, ss, system, biovar, tlag);
  return x.transpose();
}

  template <typename T0, typename T1, typename T2, typename T3,
            typename T4, typename T_biovar, typename T_tlag,
            typename
            std::enable_if_t<
              !(torsten::is_std_vector<T_biovar>::value && torsten::is_std_vector<T_tlag>::value)>* = nullptr> //NOLINT
  Eigen::Matrix <typename torsten::return_t<T0, T1, T2, T3, T4,
                                            typename torsten::value_type<T_biovar>::type,
                                            typename torsten::value_type<T_tlag>::type>::type,
                 Eigen::Dynamic, Eigen::Dynamic>
  linOdeModel(const std::vector<T0>& time,
              const std::vector<T1>& amt,
              const std::vector<T2>& rate,
              const std::vector<T3>& ii,
              const std::vector<int>& evid,
              const std::vector<int>& cmt,
              const std::vector<int>& addl,
              const std::vector<int>& ss,
              const std::vector< Eigen::Matrix<T4, -1, -1> >& system,
              const std::vector<T_biovar>& biovar,
              const std::vector<T_tlag>& tlag) {
    auto x = pmx_solve_linode(time, amt, rate, ii, evid, cmt, addl, ss,
                            system, biovar, tlag);
    return x.transpose();
  }

}
#endif
