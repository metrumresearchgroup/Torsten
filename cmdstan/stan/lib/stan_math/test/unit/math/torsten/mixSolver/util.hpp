#include <gtest/gtest.h>

#include <sstream>
#include <vector>

#include <stan/math/prim/arr/functor/integrate_ode_rk45.hpp>
#include <stan/math/torsten/ftwoCpt.hpp>
#include <stan/math/rev/core.hpp>
#include <test/unit/math/torsten/mixSolver/f2CptMix.hpp>
#include <test/unit/math/rev/arr/functor/util.hpp>
#include <test/unit/util.hpp>

// test integrate_ode_rk45 with initial positions as vars and 
// parameters as vars against finite differences. diff2 (unlike 
// in the util.hpp used to test odes) is used as a relative 
// tolerance. 
template <typename F>
void test_ode_finite_diff_vv_rel(const F& f,
                                 const double& t_in,
                                 const std::vector<double>& ts,
                                 const std::vector<double>& y_in,
                                 const std::vector<double>& theta,
                                 const std::vector<double>& x,
                                 const std::vector<int>& x_int,
                                 const double& diff,
                                 const double& rel_err) {
  std::stringstream msgs;

  std::vector<std::vector<std::vector<double> > > finite_diff_res_y(y_in.size());
  for (size_t i = 0; i < y_in.size(); i++)
    finite_diff_res_y[i] = finite_diff_initial_position(f, t_in, ts, y_in,
                                                        theta, x, x_int, i, diff);

  std::vector<std::vector<std::vector<double> > > finite_diff_res_p(theta.size());
  for (size_t i = 0; i < theta.size(); i++)
    finite_diff_res_p[i] = finite_diff_params(f, t_in, ts, y_in, theta, x,
                                              x_int, i, diff);

  std::vector<double> grads_eff;
  std::vector<stan::math::var> y_in_v;
  for (size_t i = 0; i < y_in.size(); i++)
    y_in_v.push_back(y_in[i]);

  std::vector<stan::math::var> vars = y_in_v;

  std::vector<stan::math::var> theta_v;
  for (size_t i = 0; i < theta.size(); i++)
    theta_v.push_back(theta[i]);

  for (size_t i = 0; i < theta_v.size(); i++)
    vars.push_back(theta_v[i]);

  std::vector<std::vector<stan::math::var> > ode_res;

  ode_res = stan::math::integrate_ode_rk45(f, y_in_v, t_in,
                                           ts, theta_v, x, x_int, &msgs);

  double diff2;

  for (size_t i = 0; i < ts.size(); i++) {
    for (size_t j = 0; j < y_in.size(); j++) {
      grads_eff.clear();
      ode_res[i][j].grad(vars, grads_eff);

      for (size_t k = 0; k < theta.size(); k++) {
        diff2 = std::abs(finite_diff_res_p[k][i][j]) * rel_err;
        EXPECT_NEAR(grads_eff[k+y_in.size()], finite_diff_res_p[k][i][j], diff2)
          << "Gradient of integrate_ode_rk45 failed with initial positions"
          << " unknown and parameters unknown for param at time index " << i
          << ", equation index " << j 
          << ", and parameter index: " << k;
      }
      for (size_t k = 0; k < y_in.size(); k++) {
        diff2 = std::abs(finite_diff_res_y[k][i][j]) * rel_err;
        EXPECT_NEAR(grads_eff[k], finite_diff_res_y[k][i][j], diff2)
          << "Gradient of integrate_ode_rk45 failed with initial positions"
          << " unknown and parameters known for initial position at time index " << i
          << ", equation index " << j 
          << ", and parameter index: " << k;
     }

      stan::math::set_zero_all_adjoints();
    }
  }
}

///////////////////////////////////////////////////////////////////////////////

//calculates finite diffs for mix solver with varying parameters
template <typename F>
std::vector<std::vector<double> > 
finite_diff_params_mix(const F& f,
                       const double& t_in,
                       const std::vector<double>& ts,
                       const std::vector<double>& y_in,
                       const std::vector<double>& theta,
                       const std::vector<double>& x,
                       const std::vector<int>& x_int,
                       const size_t& param_index,
                       const double& diff) {
  std::stringstream msgs;
  std::vector<double> theta_ub(theta.size());
  std::vector<double> theta_lb(theta.size());
  for (size_t i = 0; i < theta.size(); i++) {
    if (i == param_index) {
      theta_ub[i] = theta[i] + diff;
      theta_lb[i] = theta[i] - diff;
    } else {
      theta_ub[i] = theta[i];
      theta_lb[i] = theta[i];
    }
  }

  std::vector<std::vector<double> > ode_res_ub;
  std::vector<std::vector<double> > ode_res_lb;

  ode_res_ub = f2CptMix(f, y_in, t_in,
                        ts, theta_ub, x, x_int, &msgs);
  ode_res_lb = f2CptMix(f, y_in, t_in,
                        ts, theta_lb, x, x_int, &msgs);

  std::vector<std::vector<double> > results(ts.size());

  for (size_t i = 0; i < ode_res_ub.size(); i++) 
    for (size_t j = 0; j < ode_res_ub[j].size(); j++)
      results[i].push_back((ode_res_ub[i][j] - ode_res_lb[i][j]) / (2*diff));
  return results;
}

//calculates finite diffs for mix solver with varying initial positions
template <typename F>
std::vector<std::vector<double> > 
finite_diff_initial_position_mix(const F& f,
                                 const double& t_in,
                                 const std::vector<double>& ts,
                                 const std::vector<double>& y_in,
                                 const std::vector<double>& theta,
                                 const std::vector<double>& x,
                                 const std::vector<int>& x_int,
                                 const size_t& param_index,
                                 const double& diff) {
  std::stringstream msgs;
  std::vector<double> y_in_ub(y_in.size());
  std::vector<double> y_in_lb(y_in.size());
  for (size_t i = 0; i < y_in.size(); i++) {
    if (i == param_index) {
      y_in_ub[i] = y_in[i] + diff;
      y_in_lb[i] = y_in[i] - diff;
    } else {
      y_in_ub[i] = y_in[i];
      y_in_lb[i] = y_in[i];
    }
  }

  std::vector<std::vector<double> > ode_res_ub;
  std::vector<std::vector<double> > ode_res_lb;

  ode_res_ub = fTwoCpt(f, y_in_ub, t_in,
                       ts, theta, x, x_int, &msgs);
  ode_res_lb = fTwoCpt(f, y_in_lb, t_in,
                       ts, theta, x, x_int, &msgs);

  std::vector<std::vector<double> > results(ts.size());

  for (size_t i = 0; i < ode_res_ub.size(); i++) 
    for (size_t j = 0; j < ode_res_ub[j].size(); j++)
      results[i].push_back((ode_res_ub[i][j] - ode_res_lb[i][j]) / (2*diff));
  return results;
}

// test mix solver, using analytical solution to two compartment model
// and integrate_ode_rk45 with initial positions as vars and parameters
// as vars against finite differences.
template <typename F>
void test_mix_ode_finite_diff_vv(const F& f,
                                 const double& t_in,
                                 const std::vector<double>& ts,
                                 const std::vector<double>& y_in,
                                 const std::vector<double>& theta,
                                 const std::vector<double>& x,
                                 const std::vector<int>& x_int,
                                 const double& diff,
                                 const double& rel_err) {
  std::stringstream msgs;

  std::vector<std::vector<std::vector<double> > > finite_diff_res_y(y_in.size());
  for (size_t i = 0; i < y_in.size(); i++)
    finite_diff_res_y[i] = finite_diff_initial_position(f, t_in, ts, y_in,
                                                        theta, x, x_int, i, diff);

  std::vector<std::vector<std::vector<double> > > finite_diff_res_p(theta.size());
  for (size_t i = 0; i < theta.size(); i++)
    finite_diff_res_p[i] = finite_diff_params(f, t_in, ts, y_in, theta, x,
                                              x_int, i, diff);

  std::vector<double> grads_eff;
  std::vector<stan::math::var> y_in_v;
  for (size_t i = 0; i < y_in.size(); i++)
    y_in_v.push_back(y_in[i]);

  std::vector<stan::math::var> vars = y_in_v;

  std::vector<stan::math::var> theta_v;
  for (size_t i = 0; i < theta.size(); i++)
    theta_v.push_back(theta[i]);

  for (size_t i = 0; i < theta_v.size(); i++)
    vars.push_back(theta_v[i]);

  std::vector<std::vector<stan::math::var> > ode_res;

  ode_res = stan::math::integrate_ode_rk45(f, y_in_v, t_in,
                                           ts, theta_v, x, x_int, &msgs);

  double diff2;

  for (size_t i = 0; i < ts.size(); i++) {
    for (size_t j = 0; j < y_in.size(); j++) {
      grads_eff.clear();
      ode_res[i][j].grad(vars, grads_eff);

      for (size_t k = 0; k < theta.size(); k++) {
        diff2 = std::abs(finite_diff_res_p[k][i][j]) * rel_err;
        EXPECT_NEAR(grads_eff[k+y_in.size()], finite_diff_res_p[k][i][j], diff2)
          << "Gradient of integrate_ode_rk45 failed with initial positions"
          << " unknown and parameters unknown for param at time index " << i
          << ", equation index " << j 
          << ", and parameter index: " << k;
      }
      for (size_t k = 0; k < y_in.size(); k++) {
        diff2 = std::abs(finite_diff_res_y[k][i][j]) * rel_err;
        EXPECT_NEAR(grads_eff[k], finite_diff_res_y[k][i][j], diff2)
          << "Gradient of integrate_ode_rk45 failed with initial positions"
          << " unknown and parameters known for initial position at time index " << i
          << ", equation index " << j 
          << ", and parameter index: " << k;
     }

      stan::math::set_zero_all_adjoints();
    }
  }
}
