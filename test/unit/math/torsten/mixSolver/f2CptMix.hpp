#include <gtest/gtest.h>

#include <sstream>
#include <vector>

#include <stan/math/prim/arr/functor/integrate_ode_rk45.hpp>
#include <stan/math/prim/mat/fun/to_vector.hpp>
#include <stan/math/torsten/ftwoCpt.hpp>
#include <stan/math/rev/core.hpp>

/**
 * Computes an analytical solution for a mixed solver.
 * This function was specifically written for the case
 * of a Friberg-Karlsson PKPD model, with a base Two Cpt
 * PK component, and no dependence on explicit time.
 * To run fTwoCpt (code from torsten), we need to pass in
 * dt as the first argument. Here, I'll assume t0 = 0 and
 * ts is a vector of length 1 containing dt. 
 */
template <typename F, typename T1, typename T2>
std::vector<std::vector<typename stan::return_type<T1, T2>::type> >
f2CptMix(const F& f,
         const std::vector<T1> y0,
         double t0,
         const std::vector<double>& dt,
         const std::vector<T2>& theta,
         const std::vector<double>& x,
         const std::vector<int>& x_int,
         std::ostream* msgs = 0,
         double relative_tolerance = 1e-6,
         double absolute_tolerance = 1e-6,
         int max_num_steps = 1E6) { 
  using std::vector;
  using stan::math::to_vector;
  theta.push_back(y0[0]);
  theta.push_back(y0[1]);
  theta.push_back(y0[2]);

  vector<double> rate(8, 0);

  int nStates = 5;
  vector<double> init_mix(nStates);
  for (int i = 0; i < nStates; i++) init_mix[i] = y0[3 + i];

  vector<vector<double> >
    temp_PD = integrate_ode_rk45(f, init_mix, 0, dt, theta, 
                                 rate, x_int,
                                 relative_tolerance,
                                 absolute_tolerance,
                                 max_num_steps);

  Eigen::Matrix<double, Eigen::Dynamic, 1>
    tempPK = fTwoCpt(dt[0],
                     to_vector(theta),
                     to_vector(y0),
                     rate);

  vector<vector<double> > temp_mix(1);
  temp_mix[0].resize(8);
  for (int i = 0; i < 3; i++) temp_mix[0][i] = tempPK(i);
  for (int i = 0; i < 5; i++) temp_mix[0][3 + i] = temp_PD[0][i];

  return temp_PD;
}
