#ifndef STAN_MATH_TORSTEN_PKMODEL_FUNCTORS_MIX2_FUNCTOR_HPP
#define STAN_MATH_TORSTEN_PKMODEL_FUNCTORS_MIX2_FUNCTOR_HPP

#include <stan/math/prim/mat/fun/to_array_1d.hpp>
#include <stan/math/prim/mat/fun/to_vector.hpp>
#include <stan/math/rev/core.hpp>
#include <stan/math/torsten/PKModel/Pred/fTwoCpt.hpp>
#include <vector>
#include <iostream>

namespace torsten {

/**
 *  Functor for mix solver with base
 *  one compartment model.
 */
template <typename F0>
struct mix2_functor {
  F0 f0_;

  mix2_functor() { }

  explicit mix2_functor(const F0& f0) : f0_(f0) { }

  /**
   *  Returns the derivative of the base ODE system. The base 1 PK
   *  component is calculated analytically. We use the last two
   *  elements of theta to store the inital PK states (init_PK) and
   *  the last element of x_r to store the initial time.
   *
   *  Case 1: rate is fixed data and is passed through x_r.
   */
  template <typename T0, typename T1, typename T2, typename T3>
  inline
  std::vector<typename boost::math::tools::promote_args<T0, T1, T2, T3>::type>
  rate_dbl(const T0& t,
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
    // The last three elements of theta should contain the initial PK states
    init_pk[0] = theta[nTheta - 3];
    init_pk[1] = theta[nTheta - 2];
    init_pk[2] = theta[nTheta - 1];

    // Last element of x_r contains the initial time
    T0 dt = t - x_r[x_r.size() - 1];

    std::vector<T_pk> y_pk = fTwoCpt(dt, thetaPK, init_pk, x_r);
    std::vector<scalar> dydt = f0_(dt, y, y_pk, theta, x_r, x_i, pstream_);

    for (size_t i = 0; i < dydt.size(); i++)
      dydt[i] += x_r[nPK + i];

    return dydt;
  }

  /**
   *  Case 2: rate is a parameter, stored in theta.
   *  Theta contains in this order: ODE parameters, rates, and
   *  initial base PK states.
   */
  template <typename T0, typename T1, typename T2, typename T3>
  inline
  std::vector<typename boost::math::tools::promote_args<T0, T1, T2, T3>::type>
  rate_var(const T0& t,
           const std::vector<T1>& y,
           const std::vector<T2>& theta,
           const std::vector<T3>& x_r,
           const std::vector<int>& x_i,
           std::ostream* pstream_) const {
    using std::vector;
    using stan::math::to_array_1d;
    using stan::math::to_vector;
    typedef typename boost::math::tools::promote_args<T0, T1, T2, T3>::type
      scalar;
    typedef typename boost::math::tools::promote_args<T0, T2, T3>::type
      T_pk;  // return object of fTwoCpt  doesn't depend on T1

    size_t
      nTheta = theta.size(),
      nPK = 3,  // number of base PK states (equivalently rates and inits)
      nPD = y.size(),  // number of other states
      nODEparms = nTheta - 2 * nPK - nPD;  // number of ODE parameters

    // Theta first contains the base PK parameters, followed by
    // the other ODE parameters.
    int nParmsPK = 5;
    vector<T2> thetaPK(nParmsPK);
    thetaPK[0] = theta[0];  // CL
    thetaPK[1] = theta[1];  // Q
    thetaPK[2] = theta[2];  // VC
    thetaPK[3] = theta[3];  // VP
    thetaPK[4] = theta[4];  // ka

    // Next theta contains the rates for the base PK compartments...
    vector<T2> ratePK(nPK);
    for (size_t i = 0; i < nPK; i++)
      ratePK[i] = theta[nODEparms + i];

    // followed by the rates in the other compartments.
    vector<T2> ratePD(nPD);
    for (size_t i = 0; i < nPD; i++)
      ratePD[i] = theta[nODEparms + nPK + i];

    // The last elements of theta contain the initial base PK states
    vector<T2> init_pk(nPK);
    init_pk[0] = theta[nTheta - 3];
    init_pk[1] = theta[nTheta - 2];
    init_pk[2] = theta[nTheta - 1];

    // Last element of x_r contains the initial time
    T0 dt = t - x_r[x_r.size() - 1];

    vector<T_pk> y_pk = fTwoCpt(dt, thetaPK, init_pk, ratePK);
    vector<scalar> dydt = f0_(dt, y, y_pk, theta, x_r, x_i, pstream_);

    for (size_t i = 0; i < dydt.size(); i++)
      dydt[i] += ratePD[i];

    return dydt;
  }

  // Dummy operator to shut the compiler up
  // FIX ME - remove / find more elegant solution.
  template <typename T0, typename T1, typename T2, typename T3>
  inline
  std::vector<typename boost::math::tools::promote_args<T0, T1, T2, T3>::type>
  operator()(const T0& t,
             const std::vector<T1>& y,
             const std::vector<T2>& theta,
             const std::vector<T3>& x_r,
             const std::vector<int>& x_i,
             std::ostream* pstream_) const {
    std::cout << "mix2_functor: REPORT A BUG IF YOU SEE THIS." << std::endl;
    std::vector<T2> y_pk;
    return f0_(t, y, y_pk, theta, x_r, x_i, pstream_);
  }
};

}
#endif
