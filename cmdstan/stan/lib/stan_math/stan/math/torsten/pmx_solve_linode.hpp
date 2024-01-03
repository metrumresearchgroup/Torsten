#ifndef STAN_MATH_TORSTEN_LINODE_HPP
#define STAN_MATH_TORSTEN_LINODE_HPP

#include <Eigen/Dense>
#include <stan/math/torsten/ev_manager.hpp>
#include <boost/math/tools/promotion.hpp>
#include <stan/math/torsten/ev_solver.hpp>
#include <stan/math/torsten/to_array_2d.hpp>
#include <stan/math/torsten/pmx_linode_model.hpp>
#include <stan/math/torsten/pmx_check.hpp>
#include <stan/math/prim/err/check_square.hpp>
#include <vector>

namespace torsten {

  namespace {
    template<typename T>
    using torsten_matrix_dyn_t = Eigen::Matrix<T,-1,-1>;    
  }

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
stan::matrix_return_t<T0, T1, T2, T3, T4, T5, T6>
pmx_solve_linode(const std::vector<T0>& time,
                 const std::vector<T1>& amt,
                 const std::vector<T2>& rate,
                 const std::vector<T3>& ii,
                 const std::vector<int>& evid,
                 const std::vector<int>& cmt,
                 const std::vector<int>& addl,
                 const std::vector<int>& ss,
                 const std::vector< Eigen::Matrix<T4, -1, -1> >& system,
                 const std::vector<std::vector<T5> >& biovar,
                 const std::vector<std::vector<T6> >& tlag) {
  using std::vector;
  using Eigen::Dynamic;
  using Eigen::Matrix;
  using boost::math::tools::promote_args;

  static const char* function("pmx_solve_linode");
  for (size_t i = 0; i < system.size(); i++)
    stan::math::check_square(function, "system matrix", system[i]);
  int nCmt = system[0].cols();

  using ER = NONMENEventsRecord<T0, T1, T2, T3>;
  using EM = EventsManager<ER, NonEventParameters<T0, T4, torsten_matrix_dyn_t, std::tuple<T5, T6> >>;
  const ER events_rec(nCmt, time, amt, rate, ii, evid, cmt, addl, ss);

  Matrix<typename EM::T_scalar, Dynamic, Dynamic> pred =
    Matrix<typename EM::T_scalar, Dynamic, Dynamic>::Zero(events_rec.num_event_times(), EM::nCmt(events_rec));

  using model_type = torsten::PMXLinODEModel<typename EM::T_par>;
  EventSolver<model_type, EM> pr;
  pr.pred(0, events_rec, pred, dsolve::PMXAnalyiticalIntegrator(), system, biovar, tlag, nCmt);
  return pred;
}

/**
 * Overload function to allow user to pass a matrix for 
 * system.
 */
template <typename T0, typename T1, typename T2, typename T3,
          typename T4, typename T_biovar, typename T_tlag>
stan::matrix_return_t<T0, T1, T2, T3, T4, T_biovar, T_tlag>
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
          typename = require_any_not_std_vector_t<T_biovar, T_tlag> >
stan::matrix_return_t<T0, T1, T2, T3, T4, T_biovar, T_tlag>
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
stan::matrix_return_t<T0, T1, T2, T3, T4, T5, T6>
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
stan::matrix_return_t<T0, T1, T2, T3, T4, T_biovar, T_tlag>
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
            typename = require_any_not_std_vector_t<T_biovar, T_tlag> >
  stan::matrix_return_t<T0, T1, T2, T3, T4, T_biovar, T_tlag>
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
