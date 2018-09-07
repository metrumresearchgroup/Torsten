#ifndef STAN_MATH_TORSTEN_REFACTOR_LINODEMODEL_HPP
#define STAN_MATH_TORSTEN_REFACTOR_LINODEMODEL_HPP

#include <Eigen/Dense>
#include <boost/math/tools/promotion.hpp>
#include <stan/math/torsten/Pred2.hpp>
#include <stan/math/torsten/pk_linode_model.hpp>
// #include <stan/math/torsten/pk_linode_solver.hpp>
// #include <stan/math/torsten/pk_linode_solver_ss.hpp>
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
            const std::vector< Eigen::Matrix<T4, Eigen::Dynamic,
              Eigen::Dynamic> >& system,
            const std::vector<std::vector<T5> >& biovar,
            const std::vector<std::vector<T6> >& tlag) {
  using std::vector;
  using Eigen::Dynamic;
  using Eigen::Matrix;
  using boost::math::tools::promote_args;

  static const char* function("linOdeModel");
  for (size_t i = 0; i < system.size(); i++)
    stan::math::check_square(function, "system matrix", system[i]);
  int nCmt = system[0].cols();

  std::vector<T4> parameters_dummy(0);
  std::vector<std::vector<T4> > pMatrix_dummy(1, parameters_dummy);
  torsten::pmetricsCheck(time, amt, rate, ii, evid, cmt, addl, ss,
                pMatrix_dummy, biovar, tlag, function);

  PredWrapper<refactor::PKLinODEModel> pr;

#ifdef OLD_TORSTEN
  return Pred(time, amt, rate, ii, evid, cmt, addl, ss,
              pMatrix_dummy, biovar, tlag, nCmt, system,
              Pred1_linOde(), PredSS_linOde());
#else
  return pr.Pred2(time, amt, rate, ii, evid, cmt, addl, ss,
                  pMatrix_dummy, biovar, tlag, nCmt, system,
                  Pred1_linOde(), PredSS_linOde());
#endif
}

/**
 * Overload function to allow user to pass a matrix for 
 * system.
 */
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
            const Eigen::Matrix<T4, Eigen::Dynamic,
              Eigen::Dynamic>& system,
            const std::vector<std::vector<T5> >& biovar,
            const std::vector<std::vector<T6> >& tlag) {
  std::vector<Eigen::Matrix<T4, Eigen::Dynamic,
                            Eigen::Dynamic> > vec_system(1, system);

  return linOdeModel(time, amt, rate, ii, evid, cmt, addl, ss,
                     vec_system, biovar, tlag);
}

/**
 * Overload function to allow user to pass a matrix for 
 * system and a vector for biovar.
 */
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
            const Eigen::Matrix<T4, Eigen::Dynamic,
              Eigen::Dynamic>& system,
            const std::vector<T5>& biovar,
            const std::vector<std::vector<T6> >& tlag) {
  std::vector<Eigen::Matrix<T4, Eigen::Dynamic,
                            Eigen::Dynamic> > vec_system(1, system);
  std::vector<std::vector<T5> > vec_biovar(1, biovar);

  return linOdeModel(time, amt, rate, ii, evid, cmt, addl, ss,
                     vec_system, vec_biovar, tlag);
}

/**
 * Overload function to allow user to pass a matrix for 
 * system and a vector for biovar and tlag.
 */
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
            const Eigen::Matrix<T4, Eigen::Dynamic,
                                Eigen::Dynamic>& system,
            const std::vector<T5>& biovar,
            const std::vector<T6>& tlag) {
  std::vector<Eigen::Matrix<T4, Eigen::Dynamic,
                            Eigen::Dynamic> > vec_system(1, system);
  std::vector<std::vector<T5> > vec_biovar(1, biovar);
  std::vector<std::vector<T6> > vec_tlag(1, tlag);

  return linOdeModel(time, amt, rate, ii, evid, cmt, addl, ss,
                     vec_system, vec_biovar, vec_tlag);
}

/**
 * Overload function to allow user to pass a matrix for 
 * system and a vector for tlag.
 */
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
            const Eigen::Matrix<T4, Eigen::Dynamic,
                                Eigen::Dynamic>& system,
            const std::vector<std::vector<T5> >& biovar,
            const std::vector<T6>& tlag) {
  std::vector<Eigen::Matrix<T4, Eigen::Dynamic,
                            Eigen::Dynamic> > vec_system(1, system);
  std::vector<std::vector<T6> > vec_tlag(1, tlag);

  return linOdeModel(time, amt, rate, ii, evid, cmt, addl, ss,
                     vec_system, biovar, vec_tlag);
}

/**
 * Overload function to allow user to pass a vector for biovar.
 */
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
            const std::vector<Eigen::Matrix<T4, Eigen::Dynamic,
                                Eigen::Dynamic> >& system,
            const std::vector<T5>& biovar,
            const std::vector<std::vector<T6> >& tlag) {
  std::vector<std::vector<T5> > vec_biovar(1, biovar);

  return linOdeModel(time, amt, rate, ii, evid, cmt, addl, ss,
                     system, vec_biovar, tlag);
}

/**
 * Overload function to allow user to pass a vector for biovar
 * and tlag.
 */
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
            const std::vector<Eigen::Matrix<T4, Eigen::Dynamic,
                                            Eigen::Dynamic> >& system,
            const std::vector<T5>& biovar,
            const std::vector<T6>& tlag) {
  std::vector<std::vector<T5> > vec_biovar(1, biovar);
  std::vector<std::vector<T6> > vec_tlag(1, tlag);

  return linOdeModel(time, amt, rate, ii, evid, cmt, addl, ss,
                     system, vec_biovar, vec_tlag);
}

/**
 * Overload function to allow user to pass a vector for tlag.
 */
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
            const std::vector<Eigen::Matrix<T4, Eigen::Dynamic,
                                            Eigen::Dynamic> >& system,
            const std::vector<std::vector<T5> >& biovar,
            const std::vector<T6>& tlag) {
  std::vector<std::vector<T6> > vec_tlag(1, tlag);

  return linOdeModel(time, amt, rate, ii, evid, cmt, addl, ss,
                     system, biovar, vec_tlag);
}

}
#endif
