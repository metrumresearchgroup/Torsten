#ifndef STAN_MATH_TORSTEN_PKMODEL_PRED_PREDSS_ERR_HPP
#define STAN_MATH_TORSTEN_PKMODEL_PRED_PREDSS_ERR_HPP

#include <stan/math/prim/scal/err/invalid_argument.hpp>

namespace torsten {

struct PredSS_err {
  const char* function_;

  explicit PredSS_err(const char* function)
    : function_(function) { }

  /**
   * Use for function which do not handle Steady States.
   */
  template<typename T_time,
           typename T_ii,
           typename T_parameters,
           typename T_biovar,
           typename T_tlag,
           typename T_amt,
           typename T_rate>
  Eigen::Matrix<typename boost::math::tools::promote_args<T_ii,
    T_parameters>::type, Eigen::Dynamic, 1>
  operator() (const ModelParameters<T_time,
                                    T_parameters,
                                    T_biovar,
                                    T_tlag>& parameter,
              const T_amt& amt,
              const T_rate& rate,
              const T_ii& ii,
              const int& cmt) const {
    stan::math::invalid_argument(function_,
      "This function does not handle the Steady State case!",
      "", "", "");

    // return dummy matrix to suppress warning.
    typedef typename boost::math::tools::promote_args<T_ii,
      T_parameters>::type scalar;
    Eigen::Matrix<scalar, Eigen::Dynamic, 1> dummy_matrix;
    return dummy_matrix;
  }
};

}

#endif
