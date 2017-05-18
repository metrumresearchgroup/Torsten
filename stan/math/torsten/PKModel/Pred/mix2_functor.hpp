#ifndef STAN_MATH_TORSTEN_PKMODEL_PRED_MIX2_FUNCTOR_HPP
#define STAN_MATH_TORSTEN_PKMODEL_PRED_MIX2_FUNCTOR_HPP

#include <stan/math/rev/core.hpp>
#include <stan/math/fwd/core.hpp>
#include <iostream>

/**
 *  Functor for mix solver with base
 *  one compartment model.
 */
template <typename F0>
struct mix2_functor {
  F0 f0_;

  mix2_functor(const F0& f0) : f0_(f0) { }

  /**
   *  Returns the derivative of the base ODE system. The base 1 PK
   *  component is calculated analytically. We use the last two
   *  elements of theta to store the inital PK states (init_PK) and
   *  the last element of x_r to store the initial time.
   */
  template <typename T0, typename T1, typename T2, typename T3>
  inline
  std::vector<typename boost::math::tools::promote_args<T0, T1, T2, T3>::type>
  operator()(const T0& t,
             const std::vector<T1>& y,
             const std::vector<T2>& theta,
             const std::vector<T3>& x_r,
             const std::vector<int>& x_i,
             std::ostream* pstream_) const {
    using stan::math::to_array_1d;
    using stan::math::to_vector;
    typedef typename boost::math::tools::promote_args<T0, T1, T2, T3>::type
      scalar;
    typedef typename boost::math::tools::promote_args<T0, T2, T3>::type
      T_pk;  // return object of fTwoCpt  doesn't depend on T1

    // Get PK parameters
    int nParmsPK = 5;
    std::vector<T2> thetaPK(nParmsPK);
    thetaPK[0] = theta[0];  // CL
    thetaPK[1] = theta[1];  // Q
    thetaPK[2] = theta[2];  // VC
    thetaPK[3] = theta[3];  // VP
    thetaPK[4] = theta[4];  // ka

    // Get initial PK states
    int nPK = 3;
    std::vector<T2> init_pk(nPK);
    size_t nTheta = theta.size();
    // The last two components of theta should contain the initial PK states
    init_pk[0] = theta[nTheta - 3];
    init_pk[1] = theta[nTheta - 2];
    init_pk[2] = theta[nTheta - 1];
    // Last element of x_r contains the initial time
    T0 dt = t - x_r[x_r.size() - 1];

    std::vector<T_pk> y_pk = fTwoCpt(dt, thetaPK, init_pk, x_r);

    return f0_(dt, y, y_pk, theta, x_r, x_i, pstream_);
  }
};

#endif
