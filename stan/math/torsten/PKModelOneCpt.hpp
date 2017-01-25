#ifndef STAN_MATH_TORSTEN_PKMODELONECPT_HPP
#define STAN_MATH_TORSTEN_PKMODELONECPT_HPP

#include <Eigen/Dense>
#include <boost/math/tools/promotion.hpp>
#include <stan/math/torsten/PKModel/PKModel.hpp>
#include <stan/math/torsten/PKModel/dummy_ode.hpp>
#include <string>
#include <vector>

/**
 * Computes the predicted amounts in each compartment at each event
 * for a one compartment model with first oder absorption.
 *
 * @tparam T0 type of scalar for time of events.
 * @tparam T1 type of scalar for amount at each event.
 * @tparam T2 type of scalar for rate at each event.
 * @tparam T3 type of scalar for inter-dose inteveral at each event.
 * @tparam T4 type of scalars for the model parameters.
 * @tparam T5 type of scalar for bio-variability F.
 * @tparam T6 type of scalar for lag times.
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
PKModelOneCpt(const std::vector<T0>& time,
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

  PKModel model("OneCptModel");

  // Check arguments -- FIX ME: handle the new parameter arguments
/*  static const char* function("PKModelOneCpt");
  pmetricsCheck(time, amt, rate, ii, evid, cmt, addl, ss,
                pMatrix, addParm, function, model);
  for (size_t i = 0; i < pMatrix.size(); i++) {
    check_positive_finite(function, "PK parameter CL", pMatrix[i][0]);
    check_positive_finite(function, "PK parameter V2", pMatrix[i][1]);
    check_positive_finite(function, "PK parameter ka", pMatrix[i][2]);
  }

  std::string message4 = ", but must equal the number of parameters in the model: " // NOLINT
    + boost::lexical_cast<std::string>(model.GetNParameter()) + "!";
  const char* length_error4 = message4.c_str();
  if (!(pMatrix[0].size() == (size_t) model.GetNParameter()))
    stan::math::invalid_argument(function,
    "The number of parameters per event (length of a vector in the ninth argument) is", // NOLINT
    pMatrix[0].size(), "", length_error4);

  std::string message5 = ", but must equal the number of additional parameters in the model: " // NOLINT
    + boost::lexical_cast<std::string>(model.GetNAddParm()) + "!";
  const char* length_error5 = message5.c_str();
  if (!(addParm[0].size() == (size_t) model.GetNAddParm()))
    stan::math::invalid_argument(function,
    "The number of additional parameters per event (length of a vector in the tenth argument) is", // NOLINT
    addParm[0].size(), "", length_error5); */

  // Construct Pred functions for the model.
  Pred1_structure new_Pred1("OneCptModel");
  PredSS_structure new_PredSS("OneCptModel");
  Pred1 = new_Pred1;
  PredSS = new_PredSS;

  // Construct dummy matrix for last argument of pred
  Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> dummy_system;
  std::vector<Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> >
    dummy_systems(1, dummy_system);

  Matrix <typename boost::math::tools::promote_args<T0, T1, T2, T3,
    typename boost::math::tools::promote_args<T4, T5, T6>::type>::type,
    Dynamic, Dynamic> pred;
  pred = Pred(time, amt, rate, ii, evid, cmt, addl, ss, pMatrix, biovar,
              tlag, model, dummy_ode(), dummy_systems);

  return pred;
}

/*
 * Overload function to allow user to pass an std::vector for pMatrix.
 */ /*
template <typename T0, typename T1, typename T2, typename T3, typename T4>
Eigen::Matrix <typename boost::math::tools::promote_args<T0, T1, T2, T3,
  T4>::type, Eigen::Dynamic, Eigen::Dynamic>
PKModelOneCpt(const std::vector<T0>& pMatrix,
              const std::vector<T1>& time,
              const std::vector<T2>& amt,
              const std::vector<T3>& rate,
              const std::vector<T4>& ii,
              const std::vector<int>& evid,
              const std::vector<int>& cmt,
              const std::vector<int>& addl,
              const std::vector<int>& ss) {
  std::vector<std::vector<T0> > vec_pMatrix(1);
  vec_pMatrix[0] = pMatrix;

  return PKModelOneCpt(vec_pMatrix, time, amt, rate, ii, evid, cmt, addl, ss);
} */

#endif
