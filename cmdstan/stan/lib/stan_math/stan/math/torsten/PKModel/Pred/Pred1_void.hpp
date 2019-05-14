#ifndef STAN_MATH_TORSTEN_PKMODEL_PRED_PRED1_VOID_HPP
#define STAN_MATH_TORSTEN_PKMODEL_PRED_PRED1_VOID_HPP

#include <stan/math/torsten/PKModel/ModelParameters.hpp>
#include <stan/math/prim/mat/fun/Eigen.hpp>
#include <boost/math/tools/promotion.hpp>
#include <iostream>
#include <vector>

namespace torsten {

struct Pred1_void {
  Pred1_void() { }

  /**
   * Returns an empty matrix. Useful when we need to call a pred
   * operator, but don't need to compute anything. See SS_system
   * for generalODEModel_*.
   *
   * @tparam T_time type of scalar for time
   * @tparam T_rate type of scalar for rate
   * @tparam T_parameters type of scalar for model parameters
   * @tparam T_addParm type of scalar for model parameters
   * @param[in] dt time between current and previous event
   * @param[in] parameter model parameters at current event
   * @param[in] init amount in each compartment at previous event
   * @param[in] rate rate in each compartment
   * @return an eigen vector that contains predicted amount in each compartment
   *         at the current event.
   */
  template<typename T_time, typename T_rate, typename T_parameters,
           typename T_biovar, typename T_tlag, typename T_init>
  Eigen::Matrix<typename boost::math::tools::promote_args<T_time, T_rate,
    T_init, T_parameters>::type, Eigen::Dynamic, 1>
  operator()(const T_time& dt,
             const ModelParameters<T_time, T_parameters, T_biovar,
                                  T_tlag>& parameter,
             const Eigen::Matrix<T_init, 1, Eigen::Dynamic>& init,
             const std::vector<T_rate>& rate) const {
    using boost::math::tools::promote_args;
    typedef typename promote_args<T_time, T_rate, T_parameters,
                                  T_init>::type scalar;

    Eigen::Matrix<scalar, Eigen::Dynamic, 1> pred;

    return pred;
  }
};

}
#endif
