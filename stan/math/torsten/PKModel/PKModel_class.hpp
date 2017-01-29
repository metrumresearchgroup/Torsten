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
  int nParameter, nCmt;

public:
  PKModel() {
    nParameter = 0;
    nCmt = 0;
  }

  PKModel(int p_nParameter, int p_nCmt) {
    nParameter = p_nParameter;
    nCmt = p_nCmt;
  }

  // model specific constructor
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
  }

  // access function
  int GetNParameter() { return nParameter; }
  int GetNCmt() { return nCmt; }

  void Print() {
    std::cout << nParameter << " " << nCmt << std::endl;
  }
};

#endif
