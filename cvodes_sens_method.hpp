#ifndef STAN_MATH_TORSTEN_CVODES_SENS_METHOD_HPP
#define STAN_MATH_TORSTEN_CVODES_SENS_METHOD_HPP

namespace torsten {

  /**
   * Choose among three methods to calculate the
   * sensitivities:
   * CSDA: complex step derivative approximation
   * AD: automatic differentiation by Stan
   * differential quotient provided by CVODES
   **/
  enum PMXCvodesSensMethod {
    CSDA, AD, DQ
  };
}

#endif
