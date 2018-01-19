#ifndef STAN_MATH_TORSTEN_PKMODEL_UNPROMOTE_HPP
#define STAN_MATH_TORSTEN_PKMODEL_UNPROMOTE_HPP

#include <stan/math/torsten/PKModel/ModelParameters.hpp>
#include <stan/math/rev/core.hpp>
#include <stan/math/fwd/core.hpp>
#include <vector>
#include <iostream>

namespace torsten {

/**
 * Functions that converts an autodiff variable into a double.
 * The variable will either be a stan::math::var
 * or a double (in which case it will not be modified).
 *
 * @param[in] x the real to unpromote
 * @return the unpromoted real
 */
inline double unpromote(const stan::math::var& x) { return x.val(); }
inline double unpromote(const double& x) { return x; }

/**
 * Unpromote a vector.
 *
 * @param[in] x the vector of real to unpromote
 * return the unpromoted vector of reals
 */
inline
std::vector<double>
unpromote(const std::vector<stan::math::var>& x) {
  using std::vector;
  using stan::math::var;

  size_t size_x = x.size();
  vector<double> x_dbl(size_x);
  for (size_t i = 0; i < size_x; i++)
    x_dbl[i] = x[i].val();

  return x_dbl;
}

inline
std::vector<double>
unpromote(const std::vector<double>& x) {
  return x;
}

/**
 * Unpromote an object of the class model parameters.
 * Specifically, unpromote the RealParameters vector in the object.
 * @tparam T0 scalar type for time in the object
 * @tparam T1 scalar type for bio-availability in the object
 * @tparam T2 scalar type for time lag in the object
 * @param[in] object we wish to unpromote.
 * @return the same object with all members of the class unpromoted.
 */
template <typename T0, typename T1, typename T2, typename T3>
inline
ModelParameters<double, double, double, double>
unpromote(const ModelParameters<T0, T1, T2, T3>& parameters) {
  return ModelParameters<double, double, double, double>
    (unpromote(parameters.get_time()),
     unpromote(parameters.get_RealParameters()),
     unpromote(parameters.get_biovar()),
     unpromote(parameters.get_tlag()));
}

}

#endif
