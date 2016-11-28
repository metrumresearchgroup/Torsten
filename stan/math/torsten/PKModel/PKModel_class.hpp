#ifndef STAN_MATH_TORSTEN_PKMODEL_CLASS_HPP
#define STAN_MATH_TORSTEN_PKMODEL_CLASS_HPP

#include <stan/math/torsten/PKModel/functions.hpp>
#include <iostream>
#include <string>

/**
 *	This class contains some basic structural information about a
 *	compartment model.
 *
 *	nParameter: the number of parameters at each event
 *	F1Index: the index of the F1 parameter
 *	tlagIndex: the index of the tlag parameter
 *	nCmt: number of compartments
 *
 */
class PKModel {
private:
  int nParameter, F1Index, tlag1Index, nCmt;

public:
  PKModel() {
    nParameter = 0;
    F1Index = 0;
    tlag1Index = 0;
    nCmt = 0;
  }

  PKModel(int p_nParameter, int p_F1Index, int p_tlag1Index, int p_nCmt) {
    nParameter = p_nParameter;
    F1Index = p_F1Index;
    tlag1Index = p_tlag1Index;
    nCmt = p_nCmt;
  }

  // model specific constructor
  explicit PKModel(std::string p_nCmt) {
    // class constructor for one and two compartment(s) model.
    assert((p_nCmt == "OneCptModel") || (p_nCmt == "TwoCptModel"));
    if (p_nCmt == "OneCptModel") {
      nParameter = 7;
      F1Index = 3;
      tlag1Index = 5;
      nCmt = 2;
    } else if (p_nCmt == "TwoCptModel") {
      nParameter = 11;
      F1Index = 5;
      tlag1Index = 8;
      nCmt = 3;
    } else {
      nParameter = 0;
      F1Index = 0;
      tlag1Index = 0;
      nCmt = 0;
    }
  }

  // access function
  int GetNParameter() { return nParameter; }
  int GetF1Index() { return F1Index; }
  int GetTLagIndex() { return tlag1Index; }
  int GetNCmt() { return nCmt; }

  void Print() {
    std::cout << nParameter << " "
              << F1Index << " "
              << tlag1Index << " "
              << nCmt << std::endl;
  }
};

#endif
