#ifndef STAN_MATH_TORSTEN_PK_NSYS_HPP
#define STAN_MATH_TORSTEN_PK_NSYS_HPP

#include <algorithm>

namespace torsten {
  /*
   * Calculate the size of the ODE system, accounting for sensitivity 
   * equations.
   *
   * @param ncmt size of the original ODE system.
   * @param nvar number of parameters of which forward
   *        sensitivity will be calculated
   * @return the system size for each ODE unknow there will be
   *         one solution value and @c nvar gradients.
   */
  inline int pk_nsys(int ncmt, int nvar) {
    return ncmt * (nvar + 1);
  }

  /*
   * Calculate the size of the ODE system, accounting for sensitivity 
   * equations. When there are transient and steady-state
   * events in the system, the system size is the larger of
   * the two
   *
   * @param ncmt size of the original ODE system.
   * @param nvar1 number of parameters of which forward
   *        sensitivity will be calculated
   * @param nvar2 number of parameters of which forward
   *        sensitivity will be calculated
   * @return the system size for each ODE unknow there will be
   *         one solution value and @c nvar gradients.
   */
  inline int pk_nsys(int ncmt, int nvar1, int nvar2) {
    return std::max(ncmt * (nvar1 + 1), ncmt * (nvar2 + 1));
  }

}

#endif
