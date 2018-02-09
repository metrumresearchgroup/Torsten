#ifndef STAN_MATH_TORSTEN_PKMODELTWOCPT_HPP
#define STAN_MATH_TORSTEN_PKMODELTWOCPT_HPP

#include <Eigen/Dense>
#include <boost/math/tools/promotion.hpp>
#include <stan/math/torsten/PKModel/PKModel.hpp>
#include <stan/math/torsten/PKModel/Pred/Pred1_twoCpt.hpp>
#include <stan/math/torsten/PKModel/Pred/PredSS_twoCpt.hpp>
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
Eigen::Matrix <typename boost::math::tools::promote_args<T0, T1, T2, T3,
  typename boost::math::tools::promote_args<T4, T5, T6>::type>::type,
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
  using std::vector;
  using Eigen::Dynamic;
  using Eigen::Matrix;
  using boost::math::tools::promote_args;
  using stan::math::check_positive_finite;

  int nCmt = 3;
  int nParms = 5;
  static const char* function("PKModelTwoCpt");

  // Check arguments
  torsten::pmetricsCheck(time, amt, rate, ii, evid, cmt, addl, ss,
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

  return Pred(time, amt, rate, ii, evid, cmt, addl, ss,
              pMatrix, biovar, tlag,
              nCmt, dummy_systems,
              Pred1_twoCpt(), PredSS_twoCpt());
}

/**
 * Overload function to allow user to pass an std::vector for pMatrix.
 */
template <typename T0, typename T1, typename T2, typename T3, typename T4,
          typename T5, typename T6>
Eigen::Matrix <typename boost::math::tools::promote_args<T0, T1, T2, T3,
  typename boost::math::tools::promote_args<T4, T5, T6>::type>::type,
  Eigen::Dynamic, Eigen::Dynamic>
PKModelTwoCpt(const std::vector<T0>& time,
              const std::vector<T1>& amt,
              const std::vector<T2>& rate,
              const std::vector<T3>& ii,
              const std::vector<int>& evid,
              const std::vector<int>& cmt,
              const std::vector<int>& addl,
              const std::vector<int>& ss,
              const std::vector<T4>& pMatrix,
              const std::vector<std::vector<T5> >& biovar,
              const std::vector<std::vector<T6> >& tlag) {
  std::vector<std::vector<T4> > vec_pMatrix(1, pMatrix);

  return PKModelTwoCpt(time, amt, rate, ii, evid, cmt, addl, ss,
                       vec_pMatrix, biovar, tlag);
}

/**
* Overload function to allow user to pass an std::vector for pMatrix,
* and biovar.
*/
template <typename T0, typename T1, typename T2, typename T3, typename T4,
          typename T5, typename T6>
Eigen::Matrix <typename boost::math::tools::promote_args<T0, T1, T2, T3,
  typename boost::math::tools::promote_args<T4, T5, T6>::type>::type,
  Eigen::Dynamic, Eigen::Dynamic>
PKModelTwoCpt(const std::vector<T0>& time,
              const std::vector<T1>& amt,
              const std::vector<T2>& rate,
              const std::vector<T3>& ii,
              const std::vector<int>& evid,
              const std::vector<int>& cmt,
              const std::vector<int>& addl,
              const std::vector<int>& ss,
              const std::vector<T4>& pMatrix,
              const std::vector<T5>& biovar,
              const std::vector<std::vector<T6> >& tlag) {
  std::vector<std::vector<T4> > vec_pMatrix(1, pMatrix);
  std::vector<std::vector<T5> > vec_biovar(1, biovar);

  return PKModelTwoCpt(time, amt, rate, ii, evid, cmt, addl, ss,
                       vec_pMatrix, vec_biovar, tlag);
}

/**
* Overload function to allow user to pass an std::vector for pMatrix,
* and tlag.
*/
template <typename T0, typename T1, typename T2, typename T3, typename T4,
          typename T5, typename T6>
Eigen::Matrix <typename boost::math::tools::promote_args<T0, T1, T2, T3,
  typename boost::math::tools::promote_args<T4, T5, T6>::type>::type,
  Eigen::Dynamic, Eigen::Dynamic>
PKModelTwoCpt(const std::vector<T0>& time,
              const std::vector<T1>& amt,
              const std::vector<T2>& rate,
              const std::vector<T3>& ii,
              const std::vector<int>& evid,
              const std::vector<int>& cmt,
              const std::vector<int>& addl,
              const std::vector<int>& ss,
              const std::vector<T4>& pMatrix,
              const std::vector<std::vector<T5> >& biovar,
              const std::vector<T6>& tlag) {
  std::vector<std::vector<T4> > vec_pMatrix(1, pMatrix);
  std::vector<std::vector<T6> > vec_tlag(1, tlag);

  return PKModelTwoCpt(time, amt, rate, ii, evid, cmt, addl, ss,
                       vec_pMatrix, biovar, vec_tlag);
}

/**
* Overload function to allow user to pass an std::vector for pMatrix,
* biovar, and tlag.
*/
template <typename T0, typename T1, typename T2, typename T3, typename T4,
          typename T5, typename T6>
Eigen::Matrix <typename boost::math::tools::promote_args<T0, T1, T2, T3,
  typename boost::math::tools::promote_args<T4, T5, T6>::type>::type,
  Eigen::Dynamic, Eigen::Dynamic>
PKModelTwoCpt(const std::vector<T0>& time,
              const std::vector<T1>& amt,
              const std::vector<T2>& rate,
              const std::vector<T3>& ii,
              const std::vector<int>& evid,
              const std::vector<int>& cmt,
              const std::vector<int>& addl,
              const std::vector<int>& ss,
              const std::vector<T4>& pMatrix,
              const std::vector<T5>& biovar,
              const std::vector<T6>& tlag) {
  std::vector<std::vector<T4> > vec_pMatrix(1, pMatrix);
  std::vector<std::vector<T5> > vec_biovar(1, biovar);
  std::vector<std::vector<T6> > vec_tlag(1, tlag);

  return PKModelTwoCpt(time, amt, rate, ii, evid, cmt, addl, ss,
                       vec_pMatrix, vec_biovar, vec_tlag);
}

/**
* Overload function to allow user to pass an std::vector for biovar.
*/
template <typename T0, typename T1, typename T2, typename T3, typename T4,
          typename T5, typename T6>
Eigen::Matrix <typename boost::math::tools::promote_args<T0, T1, T2, T3,
  typename boost::math::tools::promote_args<T4, T5, T6>::type>::type,
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
              const std::vector<T5>& biovar,
              const std::vector<std::vector<T6> >& tlag) {
  std::vector<std::vector<T5> > vec_biovar(1, biovar);

  return PKModelTwoCpt(time, amt, rate, ii, evid, cmt, addl, ss,
                       pMatrix, vec_biovar, tlag);
}

/**
* Overload function to allow user to pass an std::vector for biovar,
* and tlag.
*/
template <typename T0, typename T1, typename T2, typename T3, typename T4,
          typename T5, typename T6>
Eigen::Matrix <typename boost::math::tools::promote_args<T0, T1, T2, T3,
  typename boost::math::tools::promote_args<T4, T5, T6>::type>::type,
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
              const std::vector<T5>& biovar,
              const std::vector<T6>& tlag) {
  std::vector<std::vector<T5> > vec_biovar(1, biovar);
  std::vector<std::vector<T6> >vec_tlag(1, tlag);

  return PKModelTwoCpt(time, amt, rate, ii, evid, cmt, addl, ss,
                       pMatrix, vec_biovar, vec_tlag);
}

/**
* Overload function to allow user to pass an std::vector for tlag.
*/
template <typename T0, typename T1, typename T2, typename T3, typename T4,
          typename T5, typename T6>
Eigen::Matrix <typename boost::math::tools::promote_args<T0, T1, T2, T3,
  typename boost::math::tools::promote_args<T4, T5, T6>::type>::type,
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
              const std::vector<T6>& tlag) {
  std::vector<std::vector<T6> > vec_tlag(1, tlag);

  return PKModelTwoCpt(time, amt, rate, ii, evid, cmt, addl, ss,
                       pMatrix, biovar, vec_tlag);
}

}
#endif
