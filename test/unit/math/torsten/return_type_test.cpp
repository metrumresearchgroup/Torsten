#include <gtest/gtest.h>
#include <stan/math/torsten/is_var.hpp>
#include <stan/math/torsten/is_std_vector.hpp>
#include <gtest/gtest.h>

using stan::math::var;

TEST(Torsten, return_t_v) {
  using t1 = typename torsten::return_t<var>::type;
  EXPECT_TRUE((std::is_same<t1, var>::value));
  EXPECT_TRUE((torsten::is_var<var>::value));
}

TEST(Torsten, return_t_d) {
  using t1 = typename torsten::return_t<double>::type;
  EXPECT_TRUE((std::is_same<t1, double>::value));
  EXPECT_FALSE((torsten::is_var<double>::value));
}

TEST(Torsten, return_t_vd) {
  using t1 = typename torsten::return_t<var, var, double, double>::type;
  using t2 = typename stan::return_type<var, var, double, double>::type;
  EXPECT_TRUE((std::is_same<t1, t2>::value));
  EXPECT_TRUE((torsten::is_var<var, var, double, double>::value));
}

TEST(Torsten, return_t_dv) {
  using t1 = typename torsten::return_t<var, double, var, double>::type;
  using t2 = typename stan::return_type<var, double, var, double>::type;
  EXPECT_TRUE((std::is_same<t1, t2>::value));
  EXPECT_TRUE((torsten::is_var<var, double, var, double>::value));
}

TEST(Torsten, return_t_vv) {
  using t1 = typename torsten::return_t<var, var>::type;
  EXPECT_TRUE((std::is_same<t1, var>::value));
  EXPECT_TRUE((torsten::is_var<var, var>::value));
}

TEST(Torsten, return_t_dd) {
  using t1 = typename torsten::return_t<double, double>::type;
  EXPECT_TRUE((std::is_same<t1, double>::value));
  EXPECT_FALSE((torsten::is_var<double, double>::value));
}

TEST(Torsten, return_t_long) {
  using t1 = typename torsten::return_t<var, double, var, double,
                                        var, double, var, double,
                                        var, double, var, double,
                                        var, double, double>::type;
  EXPECT_TRUE((std::is_same<t1, var>::value));
  EXPECT_TRUE((torsten::is_var<var, double, var, double,
               var, double, var, double,
               var, double, var, double,
               var, double, double>::value));
}

TEST(Torsten, is_std_vector) {
  EXPECT_FALSE((torsten::is_std_vector<var>::value));
  EXPECT_FALSE((torsten::is_std_vector<std::vector<double>, double>::value));
  EXPECT_FALSE((torsten::is_std_vector<var>::value));
  EXPECT_FALSE((torsten::is_std_vector<var, std::vector<double>, std::vector<var> >::value));
  EXPECT_FALSE((torsten::is_std_vector<std::vector<double>, var, std::vector<var> >::value));

  EXPECT_TRUE((torsten::is_std_vector<std::vector<var>>::value));
  EXPECT_TRUE((torsten::is_std_vector<std::vector<var>, std::vector<double> >::value));
  EXPECT_TRUE((torsten::is_std_vector<std::vector<var>, std::vector<double>, std::vector<double> >::value));
}
