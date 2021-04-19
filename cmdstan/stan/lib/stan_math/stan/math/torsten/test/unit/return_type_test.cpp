#include <gtest/gtest.h>
#include <stan/math/prim/meta/is_var.hpp>
#include <stan/math/rev/core.hpp>
#include <stan/math/torsten/meta/is_std_ode.hpp>
#include <stan/math/torsten/meta/is_eigen_ode.hpp>
#include <test/unit/math/prim/functor/harmonic_oscillator.hpp>

using stan::math::var;

// TEST(Torsten, is_std_vector) {
//   EXPECT_FALSE((torsten::is_std_vector<var>::value));
//   EXPECT_FALSE((torsten::is_std_vector<var, double>::value));
//   EXPECT_FALSE((torsten::is_std_vector<double, std::vector<double>>::value));
//   EXPECT_FALSE((torsten::is_std_vector<std::vector<double>, double>::value));
//   EXPECT_FALSE((torsten::is_std_vector<var, std::vector<double>, std::vector<var> >::value));
//   EXPECT_FALSE((torsten::is_std_vector<std::vector<double>, var, std::vector<var> >::value));

//   EXPECT_TRUE((torsten::is_std_vector<std::vector<var>>::value));
//   EXPECT_TRUE((torsten::is_std_vector<std::vector<var>, std::vector<double> >::value));
//   EXPECT_TRUE((torsten::is_std_vector<std::vector<var>, std::vector<double>, std::vector<double> >::value));
// }

// TEST(Torsten, none_std_vector) {
//   EXPECT_FALSE((torsten::none_std_vector<std::vector<double>, std::vector<double>>::value));
//   EXPECT_FALSE((torsten::none_std_vector<std::vector<double>, double>::value));
//   EXPECT_FALSE((torsten::none_std_vector<var, std::vector<double>, std::vector<var> >::value));
//   EXPECT_FALSE((torsten::none_std_vector<std::vector<double>, var, std::vector<var> >::value));
//   EXPECT_FALSE((torsten::none_std_vector<var, std::vector<double>, var>::value));
//   EXPECT_FALSE((torsten::none_std_vector<var, var, var, std::vector<double>, var>::value));

//   EXPECT_TRUE((torsten::none_std_vector<var>::value));
//   EXPECT_TRUE((torsten::none_std_vector<var, double>::value));
//   EXPECT_TRUE((torsten::none_std_vector<var, var, var, var>::value));
//   EXPECT_TRUE((torsten::none_std_vector<double, var, double>::value));
// }

TEST(Torsten, ode_signature) {
  EXPECT_TRUE((torsten::is_std_ode<harm_osc_ode_fun>::value));
  EXPECT_FALSE((torsten::is_std_ode<harm_osc_ode_fun_eigen>::value));
  EXPECT_TRUE((torsten::is_eigen_ode<harm_osc_ode_fun_eigen,
               std::vector<double>, std::vector<double>, std::vector<int>>::value));
  EXPECT_FALSE((torsten::is_eigen_ode<harm_osc_ode_fun,
                std::vector<double>, std::vector<double>, std::vector<int>>::value));
}
