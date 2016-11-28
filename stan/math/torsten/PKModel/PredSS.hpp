#ifndef STAN_MATH_TORSTEN_PKMODEL_PREDSS_HPP
#define STAN_MATH_TORSTEN_PKMODEL_PREDSS_HPP

#include <Eigen/Dense>
#include <stan/math/torsten/PKModel/Pred/PredSS_oneCpt.hpp>
#include <stan/math/torsten/PKModel/Pred/PredSS_twoCpt.hpp>
#include <stan/math/torsten/PKModel/Pred/PredSS_general_solver.hpp>
#include <stan/math/torsten/PKModel/Pred/PredSS_linCpt.hpp>
#include <iostream>
#include <string>

/**
 * The Functor of PredSS, which predicts amounts in each compartment
 * at one event, where the system is approximated to be at a steady 
 * state.
 * Defines a class of PredSS functions (functors). The key components
 * is the constructors used to create the pred function. Uses the 
 * different versions of PredSS stored in the Pred directory.
 *
 * Built-in Model types:
 *   1 - One Compartment Model with first-order absorption
 *   2 - Two Compartment Model with first-order absorption
 *   3 - Linear Compartment Model
 *   4 - General Compartment Model using numerical ODE solver
 *
 * @tparam T_time type of scalar for time
 * @tparam T_amt type of scalar for amount
 * @tparam T_rate type of scalar for rate
 * @tparam T_ii type of scalar for interdose interval
 * @tparam T_parameters type of scalar for model parameters
 * @tparam F type of ODE system function
 * @param[in] parameter model parameters at current event
 * @param[in] amt amount in specified compartment (cmt)
 * @param[in] rate rate in each compartment
 * @param[in] ii interdose interval
 * @param[in] cmt compartment number
 * @param[in] f functor for base ordinary differential equation that defines
 *   compartment model.
 * @return an eigen vector that contains predicted amount in each compartment
 *   at the current event.
 */
struct PredSS_structure {
private:
  std::string modeltype;

public:
  explicit PredSS_structure(std::string p_modeltype) {
    modeltype = p_modeltype;
  }

  // constructor for operator
  template<typename T_time, typename T_amt, typename T_rate, typename T_ii,
    typename T_parameters, typename F, typename T_system>
    Eigen::Matrix<typename boost::math::tools::promote_args< T_time, T_amt,
      T_rate, typename boost::math::tools::promote_args< T_ii, T_parameters,
      T_system>::type>::type, Eigen::Dynamic, 1>
  operator()(const ModelParameters<T_time, T_parameters, T_system>& parameter,
             const T_amt& amt,
             const T_rate& rate,
             const T_ii& ii,
             const int& cmt,
             const F& f) {
    typedef typename boost::math::tools::promote_args< T_time, T_rate,
      T_parameters, T_system>::type scalar;

    if (modeltype == "OneCptModel")
      return PredSS_one(parameter, amt, rate, ii, cmt);
    else if (modeltype == "TwoCptModel")
      return PredSS_two(parameter, amt, rate, ii, cmt);
    else if (modeltype == "GeneralCptModel_solver")
      return PredSS_general_solver(parameter, amt, rate, ii, cmt, f);
    else
      if (modeltype == "linCptModel") {
        return PredSS_linCpt(parameter, amt, rate, ii, cmt);
    } else {
      Eigen::Matrix<scalar, 1, Eigen::Dynamic> default_pred
        = Eigen::Matrix<scalar, 1, Eigen::Dynamic>::Zero(1);
      return default_pred;
    }
  }
};

#endif
