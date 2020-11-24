#include <gtest/gtest.h>
#include <stan/math/torsten/is_var.hpp>
#include <stan/math/torsten/is_std_vector.hpp>

using stan::math::var;

TEST(Torsten, return_t_v) {
  using t1 = typename stan::return_type_t<var>;
  EXPECT_TRUE((std::is_same<t1, var>::value));
  EXPECT_TRUE((torsten::is_var<var>::value));
}

TEST(Torsten, return_t_d) {
  using t1 = typename stan::return_type_t<double>;
  EXPECT_TRUE((std::is_same<t1, double>::value));
  EXPECT_FALSE((torsten::is_var<double>::value));
}

TEST(Torsten, return_t_vd) {
  EXPECT_TRUE((torsten::is_var<var, var, double, double>::value));
}

TEST(Torsten, return_t_dv) {
  EXPECT_TRUE((torsten::is_var<var, double, var, double>::value));
}

TEST(Torsten, return_t_vv) {
  EXPECT_TRUE((torsten::is_var<var, var>::value));
}

TEST(Torsten, return_t_dd) {
  EXPECT_FALSE((torsten::is_var<double, double>::value));
}

TEST(Torsten, return_t_long) {
  EXPECT_TRUE((torsten::is_var<var, double, var, double,
               var, double, var, double,
               var, double, var, double,
               var, double, double>::value));
}

TEST(Torsten, is_std_vector) {
  EXPECT_FALSE((torsten::is_std_vector<var>::value));
  EXPECT_FALSE((torsten::is_std_vector<var, double>::value));
  EXPECT_FALSE((torsten::is_std_vector<double, std::vector<double>>::value));
  EXPECT_FALSE((torsten::is_std_vector<std::vector<double>, double>::value));
  EXPECT_FALSE((torsten::is_std_vector<var, std::vector<double>, std::vector<var> >::value));
  EXPECT_FALSE((torsten::is_std_vector<std::vector<double>, var, std::vector<var> >::value));

  EXPECT_TRUE((torsten::is_std_vector<std::vector<var>>::value));
  EXPECT_TRUE((torsten::is_std_vector<std::vector<var>, std::vector<double> >::value));
  EXPECT_TRUE((torsten::is_std_vector<std::vector<var>, std::vector<double>, std::vector<double> >::value));
}

TEST(Torsten, none_std_vector) {
  EXPECT_FALSE((torsten::none_std_vector<std::vector<double>, std::vector<double>>::value));
  EXPECT_FALSE((torsten::none_std_vector<std::vector<double>, double>::value));
  EXPECT_FALSE((torsten::none_std_vector<var, std::vector<double>, std::vector<var> >::value));
  EXPECT_FALSE((torsten::none_std_vector<std::vector<double>, var, std::vector<var> >::value));
  EXPECT_FALSE((torsten::none_std_vector<var, std::vector<double>, var>::value));
  EXPECT_FALSE((torsten::none_std_vector<var, var, var, std::vector<double>, var>::value));

  EXPECT_TRUE((torsten::none_std_vector<var>::value));
  EXPECT_TRUE((torsten::none_std_vector<var, double>::value));
  EXPECT_TRUE((torsten::none_std_vector<var, var, var, var>::value));
  EXPECT_TRUE((torsten::none_std_vector<double, var, double>::value));
}
