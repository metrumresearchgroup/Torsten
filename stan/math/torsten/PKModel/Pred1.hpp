#ifndef STAN_MATH_TORSTEN_PKMODEL_PRED1_HPP
#define STAN_MATH_TORSTEN_PKMODEL_PRED1_HPP

#include <Eigen/Dense>
#include <stan/math/prim/scal/err/check_finite.hpp>
#include <stan/math/torsten/PKModel/Pred/Pred1_oneCpt.hpp>
#include <stan/math/torsten/PKModel/Pred/Pred1_twoCpt.hpp>
#include <stan/math/torsten/PKModel/Pred/Pred1_general_solver.hpp>
#include <stan/math/torsten/PKModel/Pred/Pred1_linOde.hpp>
#include <stan/math/torsten/PKModel/Pred/Pred1_mix1.hpp>
#include <stan/math/torsten/PKModel/Pred/Pred1_mix2.hpp>
#include <iostream>
#include <string>
#include <vector>

/**
 *   The Functor of Pred1, which predicts amounts in each compartment
 *   at one event. Defines a class of Pred1 functions (functors). The
 *   key components is the constructors used to create the pred operator.
 *   Calls the different versions of Pred1 stored in the Pred directory.
 *
 *   Built-in Model types:
 *     1 - One Compartment Model with first-order absorption
 *     2 - Two Compartment Model with first-order absorption
 *		 3 - General Compartment Model using numerical ODE solver
 *		 4 - EXPERIMENTAL: PKPD model using semi-analytical solver
 *
 *	 @tparam T_time type of scalar for time
 *	 @tparam T_rate type of scalar for rate
 *	 @tparam T_parameters type of scalar for model parameters
 *	 @tparam T_addParm type of scalar for additional parameters
 *	 @tparam F type of ODE system function
 *	 @param[in] dt time between current and previous event
 *	 @param[in] parameter model parameters at current event
 *	 @param[in] init amount in each compartment at previous event
 *	 @param[in] rate rate in each compartment
 *	 @param[in] f functor for base ordinary differential equation that defines
 *              compartment model.
 *   @return an eigen vector that contains predicted amount in each compartment
 *           at the current event.
 */
struct Pred1_structure {
private:
  std::string modeltype;
  double rel_tol;
  double abs_tol;
  long int max_num_steps;  // NOLINT
  std::ostream* msgs;
  std::string integratorType;

public:
  Pred1_structure(const std::string& p_modelType,
                  const double& rel_tol_p,
                  const double& abs_tol_p,
                  const long int& max_num_steps_p,  // NOLINT
                  std::ostream* msgs_p,
                  const std::string& integratorType_p)
  : modeltype(p_modelType),
    rel_tol(rel_tol_p),
    abs_tol(abs_tol_p),
    max_num_steps(max_num_steps_p),
    msgs(msgs_p),
    integratorType(integratorType_p)
  { }

  template <typename T_time, typename T_parameters, typename T_biovar,
            typename T_tlag, typename T_rate, typename F>
  Eigen::Matrix<typename boost::math::tools::promote_args<T_time,
    T_rate, T_parameters>::type, 1, Eigen::Dynamic>
    operator()(const T_time& dt,
               const ModelParameters<T_time, T_parameters, T_biovar,
                                     T_tlag>& parameter,
               const Eigen::Matrix<typename boost::math::tools::
                 promote_args<T_time, T_rate, T_parameters,
                   typename boost::math::tools::promote_args<T_biovar,
                   T_tlag>::type>::type, 1, Eigen::Dynamic>& init,
               const std::vector<T_rate>& rate,
               const F& f) {
    typedef typename boost::math::tools::promote_args<T_time, T_rate,
      T_parameters, typename boost::math::tools::promote_args<T_biovar, T_tlag>
      ::type>::type scalar;

    if (modeltype == "OneCptModel")
      return Pred1_one(dt, parameter, init, rate);
    else if (modeltype == "TwoCptModel")
      return Pred1_two(dt, parameter, init, rate);
    else if (modeltype == "generalOdeModel")
      return Pred1_general_solver(dt, parameter, init, rate, f,
                                  integrator_structure(rel_tol, abs_tol,
                                                       max_num_steps,
                                                       msgs,
                                                       integratorType));
    else if (modeltype == "linOdeModel")
      return Pred1_linOde(dt, parameter, init, rate);
    else if (modeltype == "mixOde1CptModel")
      return Pred1_mix1(dt, parameter, init, rate,
                        f, integrator_structure(rel_tol, abs_tol,
                                                max_num_steps,
                                                msgs,
                                                integratorType));
    else
      if (modeltype == "mixOde2CptModel") {
      return Pred1_mix2(dt, parameter, init, rate,
                        f, integrator_structure(rel_tol, abs_tol,
                                                max_num_steps,
                                                msgs,
                                                integratorType));
    } else {
      Eigen::Matrix<scalar, 1, Eigen::Dynamic> default_pred =
        Eigen::Matrix<scalar, 1, Eigen::Dynamic>::Zero(1);
      return default_pred;
    }
  }
};

#endif
