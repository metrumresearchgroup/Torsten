#ifndef STAN_MATH_TORSTEN_PKMODEL_PRED_PRED1_GENERAL_SOLVER_HPP
#define STAN_MATH_TORSTEN_PKMODEL_PRED_PRED1_GENERAL_SOLVER_HPP

#include <stan/math/torsten/PKModel/Pred/Pred1_oneCpt.hpp>
#include <stan/math/torsten/PKModel/Pred/unpromote.hpp>
#include <iostream>
#include <vector>

/**
 *  EXPERIMENTAL: mixed solver
 */
template<typename T_time,
         typename T_rate,
         typename T_parameters,
         typename T_biovar,
         typename T_tlag,
         typename T_system,
         typename F>
Eigen::Matrix<typename boost::math::tools::promote_args< T_time, T_rate,
  T_parameters>::type, 1, Eigen::Dynamic>
Pred1_general_solver(const T_time& dt,
                     const ModelParameters<T_time,
                                           T_parameters,
                                           T_biovar,
                                           T_tlag,
                                           T_system>& parameter,
                     const Eigen::Matrix<typename boost::math::tools::
                       promote_args<T_time,
                                    T_rate,
                                    T_parameters>::type,
                                    1, Eigen::Dynamic>& init,
                     const std::vector<T_rate>& rate,
                     const F& f) {
  using std::vector;
  typedef typename boost::math::tools::promote_args<T_time, T_rate,
    T_parameters>::type scalar;

  // Get times of current and previous events
  T_time EventTime = parameter.get_time();
  T_time InitTime = EventTime - dt;

  // Initialize return object
  Eigen::Matrix<scalar, 1, Eigen::Dynamic> pred;
  pred.resize(init.cols());

  // Compute PK amounts
  int nCmtPK = 2;  // One compartment PK model
  Eigen::Matrix<scalar, 1, Eigen::Dynamic> initPK(nCmtPK);
  for (int i = 0; i < nCmtPK; i++) initPK(i) = init(i);
  Eigen::Matrix<scalar, 1, Eigen::Dynamic> predPK;
  predPK = pred1_one(EventTime - InitTime, parameter, initPK, rate);


  // Convert time and rate to fixed data for ODE integrator
  // FIX ME - user should be able to pass rate as parameters
  vector<double> EventTime_d = vector<double>(1, static_cast<double>(0));
  EventTime_d[0] = unpromote(EventTime);
  double InitTime_d = unpromote(InitTime);
  vector<double> rate_d = vector<double>(rate.size(), static_cast<double>(0));
  for (size_t i = 0; i < rate.size(); i++) rate_d[i] = unpromote(rate[i]);

  // contains both PK and PD parameters
  vector<T_parameters> theta = parameter.get_RealParameters();

  // Create initial vector for PD system (i.e. need to exclude PK)
  vector<scalar> init_vector = vector<scalar>(init.cols() - nCmtPK, scalar(0));
  for (size_t i = 0; i < init_vector.size(); i++)
    init_vector[i] = init(i + nCmtPK);

  // store PK initial states in the parameters theta
  for (int i = 0; i < nCmtPK; i++) {
    theta.pushback(init(i));
  }

  if (EventTime_d[0] == InitTime_d) { pred = init;
  } else {
    vector< vector<scalar> > pred_V;
    vector<int> idummy;
    pred_V = pmetrics_solver(f, init_vector, InitTime_d,
                             EventTime_d, theta,
                             rate_d, idummy);

    for (size_t i = nCmtPK; i < init.cols(); i++)
      pred(i) = pred_V[0][i];
  }

  return pred;
}


#endif
