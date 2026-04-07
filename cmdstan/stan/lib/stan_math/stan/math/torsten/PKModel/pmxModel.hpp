#ifndef STAN_MATH_TORSTEN_PKMODEL_CLASS_HPP
#define STAN_MATH_TORSTEN_PKMODEL_CLASS_HPP

#include <stan/math/torsten/PKModel/functions.hpp>
#include <iostream>
#include <string>

namespace torsten {

/**
 *  Use this class to store information about a
 *  compartment model. On top of basic information,
 *  contains instructions on how to construct the
 *  pred1 and predSS functors.
 *
 *  Members:
 *	nParam_ the number of parameters per event
 *  nCmt_ number of compartments
 *  pred1Type_
 *  predSSType_
 *  integratorType_
 *  rel_tol_ relative tolerance of ODE integrator
 *  abs_tol_ absolute tolerance of ODE integrator
 *  max_num_steps maximum number of steps of ODE integrator.
 */
class pmxModel {
private:
  int nParm_;
  int nCmt_;
  std::string pred1Type_;
  std::string predSSType_;
  std::string integratorType_;
  std::ostream* msgs_;
  double rel_tol_;
  double abs_tol_;
  long int max_num_steps_;  // NOLINT

public:
  pmxModel() {
    nParm_ = 0;
    nCmt_ = 0;
    pred1Type_ = "default";
    predSSType_ = "default";
    integratorType_ = "default";
    msgs_ = 0;
    rel_tol_ = 1e-6;
    abs_tol_ = 1e-6;
    max_num_steps_ = 1e+6;
  }

  pmxModel(const int& nParm,
           const int& nCmt,
           const std::string& predType) {
    nParm_ = nParm;
    nCmt_ = nCmt;
    pred1Type_ = predType;
    predSSType_ = predType;

    integratorType_ = "default";
    rel_tol_ = 0;
    abs_tol_ = 0;
    max_num_steps_ = 0;
  }

  pmxModel(const int& nParm,
           const int& nCmt,
           const std::string& pred1Type,
           const std::string& predSSType,
           const std::string& integratorType,
           std::ostream* msgs,
           const double& rel_tol,
           const double& abs_tol,
           const long int& max_num_steps) {  // NOLINT
    nParm_ = nParm;
    nCmt_ = nCmt;
    pred1Type_ = pred1Type;
    predSSType_ = predSSType;

    integratorType_ = integratorType;
    msgs_ = msgs;
    rel_tol_ = rel_tol;
    abs_tol_ = abs_tol;
    max_num_steps_ = max_num_steps;
  }

  // access function
  int GetNParm() const { return nParm_; }
  int GetNCmt() const { return nCmt_; }
  std::string GetPred1Type() const { return pred1Type_; }
  std::string GetPredSSType() const { return predSSType_; }
  std::string GetIntegratorType() const { return integratorType_; }
  std::ostream* GetMsgs() const { return msgs_; }
  double GetRelTol() const { return rel_tol_; }
  double GetAbsTol() const { return abs_tol_; }
  long int GetMaxNumSteps() const { return max_num_steps_; }  // NOLINT

  void Print() {
    std::cout << nParm_ << " " << nCmt_ << " "
              << pred1Type_ << " " << predSSType_ << " "
              << integratorType_ << " "
              << msgs_ << " "
              << rel_tol_ << " " << abs_tol_ << " "
              << max_num_steps_ << std::endl;
  }
};

}

#endif
