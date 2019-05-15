#ifndef STAN_MATH_TORSTEN_MPI_PMX_POPULATION_CHECK_HPP
#define STAN_MATH_TORSTEN_MPI_PMX_POPULATION_CHECK_HPP

#include <stan/math/prim/scal/err/check_greater_or_equal.hpp>
#include <stan/math/prim/scal/err/check_finite.hpp>
#include <stan/math/prim/scal/err/check_not_nan.hpp>
#include <stan/math/prim/scal/err/check_nonnegative.hpp>
#include <stan/math/prim/scal/err/check_consistent_sizes.hpp>
#include <stan/math/prim/arr/err/check_nonzero_size.hpp>
#include <stan/math/prim/arr/err/check_ordered.hpp>
#include <stan/math/torsten/PKModel/pmxModel.hpp>
#include <stan/math/prim/scal/err/invalid_argument.hpp>
#include <Eigen/Dense>
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
void pmx_population_check(const std::vector<int>& len,
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
                          const char* function) {
  using std::vector;
  using std::string;
  using Eigen::Dynamic;
  using stan::math::invalid_argument;
  using stan::math::check_greater_or_equal;
  using stan::math::check_consistent_sizes;
  using stan::math::check_finite;
  using stan::math::check_not_nan;
  using stan::math::check_nonnegative;

  check_greater_or_equal(function, "time size", time.size(), size_t(1));
  check_consistent_sizes(function, "amt",  amt  , "time", time);
  check_consistent_sizes(function, "rate", rate , "time", time);
  check_consistent_sizes(function, "evid", evid , "time", time);
  check_consistent_sizes(function, "cmt",  cmt  , "time", time);
  check_consistent_sizes(function, "ii",   ii   , "time", time);
  check_consistent_sizes(function, "addl", addl , "time", time);
  check_consistent_sizes(function, "ss",   ss   , "time", time);

  stan::math::check_nonzero_size(function, "times", time);

  check_finite(function, "time", time);
  check_finite(function, "amt",  amt );
  check_finite(function, "rate", rate);
  check_finite(function, "ii",   ii  );

  check_not_nan(function, "time", time);
  check_not_nan(function, "amt",  amt );
  check_not_nan(function, "rate", rate);
  check_not_nan(function, "ii",   ii  );

  check_nonnegative(function, "amt",  amt );
  check_nonnegative(function, "rate", rate);
  check_nonnegative(function, "ii",   ii  );

  // total size should be match
  size_t s = 0;
  for (auto const & i : len) s += i;
  check_greater_or_equal(function, "time size", time.size(), s);  

  if (pMatrix.size() > len.size()) {
    check_consistent_sizes(function, "parameters", pMatrix, "time", time);
  } else {
    check_consistent_sizes(function, "parameters", pMatrix, "len", len);
  }
  for (auto const & p : pMatrix) {
    check_finite(function, "parameters", p);
    check_not_nan(function, "parameters", p);
  }

  if (biovar.size() > len.size()) {
    check_consistent_sizes(function, "bioavailability", biovar, "time", time);
  } else {
    check_consistent_sizes(function, "bioavailability", biovar, "len", len);
  }
  for (auto const & p : biovar) {
    check_nonnegative(function, "bioavailability", p);
    check_finite(function, "bioavailability", p);
    check_not_nan(function, "bioavailability", p);
  }

  if (tlag.size() > len.size()) {
    check_consistent_sizes(function, "lag time", tlag, "time", time);
  } else {
    check_consistent_sizes(function, "lag time", tlag, "len", len);
  }
  for (auto const & p : tlag) {
    check_nonnegative(function, "lag time", p);
    check_finite(function, "lag time", p);
    check_not_nan(function, "lag time", p);
  }
}

}    // torsten namespace

#endif
