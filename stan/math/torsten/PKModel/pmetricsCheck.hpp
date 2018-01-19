#ifndef STAN_MATH_TORSTEN_PKMODEL_PMETRICSCHECK_HPP
#define STAN_MATH_TORSTEN_PKMODEL_PMETRICSCHECK_HPP

#include <Eigen/Dense>
#include <stan/math/torsten/PKModel/pmxModel.hpp>
#include <stan/math/prim/scal/err/invalid_argument.hpp>
#include <boost/lexical_cast.hpp>
#include <vector>
#include <string>
#include <iostream>

namespace torsten {
/**
 * Checks that the arguments the user inputs in the Torsten
 * functions are valid.
 *
 * @tparam T0 type of scalar for time of events.
 * @tparam T1 type of scalar for amount at each event.
 * @tparam T2 type of scalar for rate at each event.
 * @tparam T3 type of scalar for inter-dose inteveral at each event.
 * @tparam T4 type of scalar for model parameters
 * @tparam T5 type of scalar for bio-variability
 * @tparam T6 type of scalar for lag times
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
 * @param[in] pMatrix parameters at each event
 * @param[in] bio-variability at each event
 * @param[in] lag times at each event
 * @param[in] function The name of the function for which the check is being
 *                     performed.
 * @param[in] model object that contains basic structural information about
 *                  a compartment model.
 * @return void
 *
 */
template <typename T0, typename T1, typename T2, typename T3, typename T4,
  typename T5, typename T6>
void pmetricsCheck(const std::vector<T0>& time,
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
                   const char* function) {
  using std::vector;
  using std::string;
  using Eigen::Dynamic;
  using stan::math::invalid_argument;

  if (!(time.size() > 0)) invalid_argument(function,
    "length of time vector,", time.size(), "",
    "needs to be positive and greater than 0!");

  std::string message = ", but must be the same as the length of the time array: "  // NOLINT
    + boost::lexical_cast<string>(time.size()) + "!";
    const char* length_error = message.c_str();

  if (!(amt.size() == time.size())) invalid_argument(function,
    "the length of the amount (amt) array is", amt.size(), "",
    length_error);
  if (!(rate.size() == time.size())) invalid_argument(function,
    "the length of the rate array is", rate.size(), "",
    length_error);
  if (!(evid.size() == time.size())) invalid_argument(function,
    "the length of the event ID (evid) array is", evid.size(), "",
    length_error);
  if (!(cmt.size() == time.size())) invalid_argument(function,
    "the length of the compartment (cmt) array is", cmt.size(), "",
    length_error);

  std::string message2 = ", but must be either 1 or the same as the length of the time array: "  // NOLINT
    + boost::lexical_cast<string>(time.size()) + "!";
    const char* length_error2 = message2.c_str();

  if (!(ii.size() == time.size()) || (ii.size() == 1)) invalid_argument(
    function,
    "the length of the interdose interval (ii) array is", ii.size(), "",
    length_error2);
  if (!(addl.size() == time.size()) || (addl.size() == 1)) invalid_argument(
    function,
    "the length of the additional dosing (addl) array is", ii.size(), "",
    length_error2);
  if (!(ss.size() == time.size()) || (ss.size() == 1)) invalid_argument(
    function,
    "the length of the steady state approximation (ss) array is", ss.size(),
    "", length_error2);

  std::string message3 = ", but must be the same as the length of the additional dosing (addl) array: " // NOLINT
    + boost::lexical_cast<string>(addl.size()) + "!";
    const char* length_error3 = message3.c_str();
    if (!(ss.size() == time.size()) || (ss.size() == 1)) invalid_argument(
      function,
      "the length of steady state approximation (ss) array is", ss.size(), "",
      length_error3);

  // TEST ARGUMENTS FOR PARAMETERS
  static const char* noCheck("linOdeModel");
  if (strcmp(function, noCheck) != 0) {
    if (!((pMatrix.size() == time.size()) || (pMatrix.size() == 1)))
      invalid_argument(function, "length of the parameter (2d) array,",
        pMatrix.size(), "", length_error2);
    if (!(pMatrix[0].size() > 0)) invalid_argument(function,
      "the number of parameters per event is", pMatrix[0].size(),
      "", " but must be greater than 0!");
  }

  if (!((biovar.size() == time.size()) || (biovar.size() == 1)))
    invalid_argument(function, "length of the biovariability parameter (2d) array,",  // NOLINT
      biovar.size(), "", length_error2);
  if (!(biovar[0].size() > 0)) invalid_argument(function,
    "the number of biovariability parameters per event is", biovar[0].size(),
    "", " but must be greater than 0!");

  if (!((tlag.size() == time.size()) || (tlag.size() == 1)))
    invalid_argument(function, "length of the lag times (2d) array,",  // NOLINT
                     tlag.size(), "", length_error2);
  if (!(tlag[0].size() > 0)) invalid_argument(function,
      "the number of lagtimes parameters per event is", tlag[0].size(),
      "", " but must be greater than 0!");
}

}    // torsten namespace

#endif
