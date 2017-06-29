#ifndef STAN_MATH_TORSTEN_PKMODEL_PREDSS_HPP
#define STAN_MATH_TORSTEN_PKMODEL_PREDSS_HPP

#include <Eigen/Dense>
#include <stan/math/torsten/PKModel/Pred/PredSS_oneCpt.hpp>
#include <stan/math/torsten/PKModel/Pred/PredSS_twoCpt.hpp>
#include <stan/math/torsten/PKModel/Pred/PredSS_general_solver.hpp>
#include <stan/math/torsten/PKModel/Pred/PredSS_linOde.hpp>
#include <stan/math/torsten/PKModel/integrator.hpp>
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
 * @tparam T_addParm type of scalar for additional parameters
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
  std::string modeltype_;
  double rel_tol_;
  double abs_tol_;
  long int max_num_steps_;  // NOLINT
  std::ostream* msgs_;
  std::string integratorType_;
  int nCmt_;

public:
  PredSS_structure(const std::string& modeltype,
                   const double& rel_tol,
                   const double& abs_tol,
                   const long int& max_num_steps,  // NOLINT
                   std::ostream* msgs,
                   const std::string& integratorType,
                   const int& nCmt)
    : modeltype_(modeltype),
      rel_tol_(rel_tol),
      abs_tol_(abs_tol),
      max_num_steps_(max_num_steps),
      msgs_(msgs),
      integratorType_(integratorType),
      nCmt_(nCmt)
  { }

  // constructor for operator
  template<typename T_time,
           typename T_amt,
           typename T_rate,
           typename T_ii,
           typename T_parameters,
           typename T_biovar,
           typename T_tlag,
           typename F>
  Eigen::Matrix<typename boost::math::tools::promote_args<T_amt, T_rate,
    T_ii, T_parameters>::type, 1, Eigen::Dynamic> 
  operator()(const ModelParameters<T_time,T_parameters, T_biovar,
                                   T_tlag>& parameter,
             const T_amt& amt,
             const T_rate& rate,
             const T_ii& ii,
             const int& cmt,
             const F& f) {
    typedef typename boost::math::tools::promote_args<T_amt, T_rate,
      T_ii, T_parameters>::type scalar;

    if (modeltype_ == "OneCptModel")
      return PredSS_one(parameter, amt, rate, ii, cmt);
    else if (modeltype_ == "TwoCptModel")
      return PredSS_two(parameter, amt, rate, ii, cmt);
    else if (modeltype_ == "generalOdeModel")
      return PredSS_general_solver(parameter, amt, rate, ii, cmt, f,
                                   integrator_structure(rel_tol_, abs_tol_,
                                                        max_num_steps_,
                                                        msgs_,
                                                        integratorType_), 
                                   nCmt_);
    else if (modeltype_ == "linOdeModel")
      return PredSS_linOde(parameter, amt, rate, ii, cmt);
    else {
      std::cout << "IF YOU SEE THIS REPORT AN ISSUE (PREDSS)" << std::endl;
      Eigen::Matrix<scalar, 1, Eigen::Dynamic> default_pred
        = Eigen::Matrix<scalar, 1, Eigen::Dynamic>::Zero(1);
      return default_pred;
    }
  }
};

#endif
