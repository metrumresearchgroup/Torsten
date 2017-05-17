#ifndef STAN_MATH_TORSTEN_PKMODEL_CLASS_HPP
#define STAN_MATH_TORSTEN_PKMODEL_CLASS_HPP

#include <stan/math/torsten/PKModel/functions.hpp>
#include <iostream>
#include <string>

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
    rel_tol_ = 1e-6;
    abs_tol_ = 1e-6;
    max_num_steps_ = 1e+6;
  }

  pmxModel(int nParm, int nCmt, std::string pred1Type,
           std::string predSSType) {
    nParm_ = nParm;
    nCmt_ = nCmt;
    pred1Type_ = pred1Type;
    predSSType_ = predSSType;

    integratorType_ = "default";
    rel_tol_ = 1e-6;
    abs_tol_ = 1e-6;
    max_num_steps_ = 1e+6;
  }

  pmxModel(int nParm, int nCmt, std::string pred1Type,
           std::string predSSType, std::string integratorType,
           double rel_tol, double abs_tol,
           long int max_num_steps) {  // NOLINT
    nParm_ = nParm;
    nCmt_ = nCmt;
    pred1Type_ = pred1Type;
    predSSType_ = predSSType;

    integratorType_ = integratorType;
    rel_tol_ = rel_tol;
    abs_tol_ = abs_tol;
    max_num_steps_ = max_num_steps;
  }

  // model specific constructor
  /*
  explicit PKModel(std::string p_nCmt) {
    // class constructor for one and two compartment(s) model.
    assert((p_nCmt == "OneCptModel") || (p_nCmt == "TwoCptModel"));
    if (p_nCmt == "OneCptModel") {
      nParameter = 3;
      nCmt = 2;
    } else if (p_nCmt == "TwoCptModel") {
      nParameter = 5;
      nCmt = 3;
    } else {
      nParameter = 0;
      nCmt = 0;
    }
  } */

  // access function
  int GetNParm() const { return nParm_; }
  int GetNCmt() const { return nCmt_; }
  std::string GetPred1Type() const { return pred1Type_; }
  std::string GetPredSSType() const { return predSSType_; }
  std::string GetIntegratorType() const { return integratorType_; }
  double GetRelTol() const { return rel_tol_; }
  double GetAbsTol() const { return abs_tol_; }
  long int GetMaxNumSteps() const { return max_num_steps_; }  // NOLINT

  void Print() {
    std::cout << nParm_ << " " << nCmt_ << " "
              << pred1Type_ << " " << predSSType_ << " "
              << integratorType_ << " "
              << rel_tol_ << " " << abs_tol_ << " "
              << max_num_steps_ << std::endl;
  }
};

#endif
