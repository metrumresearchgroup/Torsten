#include <stan/math.hpp>
#include <stan/math/rev/core.hpp>
#include <test/unit/math/rev/mat/fun/util.hpp>
#include <stan/math/torsten/dsolve/group_functor.hpp>
#include <test/unit/math/torsten/pmx_ode_test_fixture.hpp>
#include <stan/math/rev/mat/functor/integrate_ode_bdf.hpp>
#include <nvector/nvector_serial.h>
#include <boost/mpi.hpp>
#include <test/unit/util.hpp>
#include <gtest/gtest.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <memory>
#include <iostream>
#include <fstream>
#include <limits>
#include <vector>
#include <chrono>
#include <ctime>
#include <random>

struct harm_osc_ode_fun2 {
  template <typename T0, typename T1, typename T2>
  inline std::vector<typename stan::return_type<T1, T2>::type>
  operator()(const T0& t_in, const std::vector<T1>& y_in,
             const std::vector<T2>& theta, const std::vector<double>& x,
             const std::vector<int>& x_int, std::ostream* msgs) const {
    if (y_in.size() != 2)
      throw std::domain_error(
          "this function was called with inconsistent state");

    std::vector<typename stan::return_type<T1, T2>::type> res;
    res.push_back(y_in.at(1));
    res.push_back(-y_in.at(0) - 2.0 * theta.at(0) * y_in.at(1));

    return res;
  }
};

template<typename... Args>
inline auto torsten::dsolve::pmx_ode_group_mpi_functor::operator()(Args&&... args) const {
    if (id == 0) { harm_osc_ode_fun f; return f(std::forward<Args>(args)...); }
    if (id == 1) { harm_osc_ode_fun2 f; return f(std::forward<Args>(args)...); }

    // return default
    harm_osc_ode_fun f;
    return f(std::forward<Args>(args)...);
}

using stan::math::integrate_ode_rk45;
using torsten::pmx_integrate_ode_rk45;
using stan::math::var;
using std::vector;

TEST_F(TorstenOdeTest_sho, fwd_sensitivity_theta_y0) {
  std::vector<var> theta_var1 = stan::math::to_var(theta);
  std::vector<var> y0_var1 = stan::math::to_var(y0);
  std::vector<var> theta_var2 = stan::math::to_var(theta);
  std::vector<var> y0_var2 = stan::math::to_var(y0);

  harm_osc_ode_fun f1;
  harm_osc_ode_fun2 f2;

  {
    torsten::dsolve::pmx_ode_group_mpi_functor f_gen(0);
    std::vector<std::vector<stan::math::var>> y1 = pmx_integrate_ode_rk45(f1, y0_var1, t0, ts, theta_var1, x_r, x_i);
    std::vector<std::vector<stan::math::var>> y2 = pmx_integrate_ode_rk45(f_gen, y0_var2, t0, ts, theta_var2, x_r, x_i);
    torsten::test::test_grad(y0_var1, y0_var2, y1, y2, 1.E-10, 1.E-8);
    torsten::test::test_grad(theta_var1, theta_var2, y1, y2, 1.E-10, 1.E-8);    
  }

  {
    torsten::dsolve::pmx_ode_group_mpi_functor f_gen(1);
    std::vector<std::vector<stan::math::var>> y1 = pmx_integrate_ode_rk45(f2, y0_var1, t0, ts, theta_var1, x_r, x_i);
    std::vector<std::vector<stan::math::var>> y2 = pmx_integrate_ode_rk45(f_gen, y0_var2, t0, ts, theta_var2, x_r, x_i);
    torsten::test::test_grad(y0_var1, y0_var2, y1, y2, 1.E-10, 1.E-8);
    torsten::test::test_grad(theta_var1, theta_var2, y1, y2, 1.E-10, 1.E-8);    
  }

}
