#include <gtest/gtest.h>
#include <stan/math/prim/meta/is_var.hpp>
#include <stan/math/torsten/meta/is_std_ode.hpp>
#include <stan/math/torsten/meta/is_eigen_ode.hpp>
#include <stan/math/torsten/meta/is_nl_system.hpp>
#include <test/unit/math/prim/functor/harmonic_oscillator.hpp>
#include <test/unit/math/rev/functor/util_algebra_solver.hpp>

using torsten::nl_system_adaptor;
TEST(Torsten, nonlinear_system_signature) {
  static_assert(torsten::is_nl_system<nl_system_adaptor<simple_eq_functor>,
                Eigen::VectorXd, std::vector<double>, std::vector<int>>::value,
                "Wrong function signature");
  static_assert(torsten::is_nl_system<nl_system_adaptor<non_linear_eq_functor>,
                Eigen::VectorXd, std::vector<double>, std::vector<int>>::value,
                "Wrong function signature");
  static_assert(!torsten::is_nl_system<simple_eq_functor,
                Eigen::VectorXd, std::vector<double>, std::vector<int>>::value,
                "Wrong function signature");
  static_assert(!torsten::is_nl_system<non_linear_eq_functor,
                Eigen::VectorXd, std::vector<double>, std::vector<int>>::value,
                "Wrong function signature");
}
