#include <gtest/gtest.h>
#include <test/unit/lang/utility.hpp>

TEST(lang_parser, univariate_integrate_good) {
  test_parsable("torsten/univariate_integral_good");
}
TEST(lang_parser, univariate_integrate_rk45_bad) {
  test_throws("torsten/univariate_integral/bad_fun_type_rk45",
              "first argument to univariate_integral_rk45 "
              "must be the name of a function with signature");
  test_throws("torsten/univariate_integral/bad_t0_rk45",
              "second argument to univariate_integral_rk45 "
              "must have type real for time limit;");
  test_throws("torsten/univariate_integral/bad_t1_rk45",
              "third argument to univariate_integral_rk45 "
              "must have type real for time limit;");
  test_throws("torsten/univariate_integral/bad_theta_rk45",
              "fourth argument to univariate_integral_rk45 "
              "must have type real[] for parameters;");
  test_throws("torsten/univariate_integral/bad_x_r_rk45",
              "fifth argument to univariate_integral_rk45 "
              "must have type real[] for real data;");
  test_throws("torsten/univariate_integral/bad_x_i_rk45",
              "sixth argument to univariate_integral_rk45 "
              "must have type int[] for integer data;");
}
TEST(lang_parser, univariate_integrate_bdf_bad) {
  test_throws("torsten/univariate_integral/bad_fun_type_bdf",
              "first argument to univariate_integral_bdf "
              "must be the name of a function with signature");
  test_throws("torsten/univariate_integral/bad_t0_bdf",
              "second argument to univariate_integral_bdf "
              "must have type real for time limit;");
  test_throws("torsten/univariate_integral/bad_t1_bdf",
              "third argument to univariate_integral_bdf "
              "must have type real for time limit;");
  test_throws("torsten/univariate_integral/bad_theta_bdf",
              "fourth argument to univariate_integral_bdf "
              "must have type real[] for parameters;");
  test_throws("torsten/univariate_integral/bad_x_r_bdf",
              "fifth argument to univariate_integral_bdf "
              "must have type real[] for real data;");
  test_throws("torsten/univariate_integral/bad_x_i_bdf",
              "sixth argument to univariate_integral_bdf "
              "must have type int[] for integer data;");
}

TEST(lang_parser, generalOdeModel_good) {
  test_parsable("torsten/generalOdeModel_good");
}

TEST(lang_parser, mixOdeModel_good) {
  test_parsable("torsten/mixOdeModel_good");
}

 // old tests
TEST(lang_parser, generalCptModel) {
  test_parsable("torsten/generalCptModel");
}
TEST(lang_parser, PKModelOneCpt_function_signatures) {
    test_parsable("torsten/PKModelOneCpt");
}
TEST(lang_parser, PKModelTwoCpt_function_signatures) {
    test_parsable("torsten/PKModelTwoCpt");
}
TEST(lang_parser, linOdeModel_function_signatures) {
    test_parsable("torsten/linOdeModel");
}
TEST(lang_parser, mixOde1CptModel_function_signatures) {
    test_parsable("torsten/mixOde1CptModel");
}
TEST(lang_parser, mixOde2CptModel_function_signatures) {
    test_parsable("torsten/mixOde2CptModel");
}
