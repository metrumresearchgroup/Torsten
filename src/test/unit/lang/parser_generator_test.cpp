#include <fstream>
#include <iostream>
#include <gtest/gtest.h>

void read_file(const std::string& path, std::string& contents) {
  std::stringstream s;
  std::fstream f(path.c_str());
  while (f.good()) {
    s << static_cast<char>(f.get());
  }
  contents = s.str();
}

void test_pg(const std::string& program_name,
             const std::string& expected_substring) {
  std::string path = "test/test-models/good/parser-generator";
  path += "/";
  path += program_name;
  path += ".hpp";

  std::string hpp_code;
  read_file(path, hpp_code);

  EXPECT_TRUE(hpp_code.find(expected_substring) != std::string::npos)
      << "program_name: " << program_name << std::endl
      << "expected_substring: " << expected_substring;
}

int count_occurrences(const std::string target, const std::string s) {
  size_t count = 0;
  size_t offset = 0;
  while (true) {
    offset = s.find(target, offset);
    if (offset == std::string::npos) {
      break;
    }
    offset += target.size();
    ++count;
  }
  return count;
}

void test_pg_count(const std::string& program_name,
                   const std::string& expected_substring,
                   const int expected_count) {
  std::string path = "test/test-models/good/parser-generator";
  path += "/";
  path += program_name;
  path += ".hpp";

  std::string hpp_code;
  read_file(path, hpp_code);

  EXPECT_EQ(expected_count, count_occurrences(expected_substring, hpp_code));
}

TEST(unitLang, simpleTest) {
  test_pg("user-function-struct-const",
          "operator()(const T0__& x, std::ostream* pstream__) const {");
}

TEST(unitLang, odeTest) {
  std::string expected;
  expected
      = "stan::math::assign(y_hat, "
        "integrate_ode_rk45(sho_functor__(), y0, t0, ts, theta, x, x_int, "
        "pstream__));";
  test_pg("ode", expected);
  test_pg_count("ode", expected, 1);
}

TEST(unitLang, algebra_solverTest) {
  std::string expected;
  expected
      = "stan::math::assign(y_s, "
        "algebra_solver(algebra_system_functor__(), y, theta, x_r, x_i, "
        "pstream__));";
  test_pg("algebra_solver", expected);
  test_pg_count("algebra_solver", expected, 1);
}

// torsten
TEST(unitLang, univariate_integral) {
  std::string expected;
  expected = "stan::math::assign(y_hat, "
    "univariate_integral_rk45(foo_functor__(), t0, t1, theta, x, x_int, pstream__));";
  test_pg("univariate_integral", expected);
  test_pg_count("univariate_integral", expected, 1);
}

TEST(unitLang, pmx_solve_rk45) {
  std::string expected;
  expected = "stan::math::assign(x, pmx_solve_rk45(twoCptModelODE_functor__(), nCmt, time, amt, rate, ii, evid, cmt, addl, ss, theta_data, biovar_data, tlag_data, pstream__));";
  test_pg("pmx_solve_rk45", expected);
  test_pg_count("pmx_solve_rk45", expected, 1);

  expected = "stan::math::assign(x, pmx_solve_rk45(twoCptModelODE_functor__(), nCmt, time, amt, rate, ii, evid, cmt, addl, ss, theta_data, biovar_data, tlag_data, 1.e-5, 1.e-6, 1000, pstream__));";
  test_pg("pmx_solve_rk45", expected);
  test_pg_count("pmx_solve_rk45", expected, 1);

  expected = "stan::math::assign(x, pmx_solve_rk45(twoCptModelODE_functor__(), nCmt, time, amt, rate, ii, evid, cmt, addl, ss, theta_data, biovar_data, tlag_data, 1.e-5, 1.e-6, 1000, 1.e-5, 1.e-6, 100, pstream__));";
  test_pg("pmx_solve_rk45", expected);
  test_pg_count("pmx_solve_rk45", expected, 1);
}

TEST(unitLang, pmx_solve_group) {
  std::string expected;
  expected = "stan::math::assign(x, pmx_solve_group_bdf(twoCptModelODE_functor__(), 3, len, time, amt, rate, ii, evid, cmt, addl, ss, theta, biovar, tlag, pstream__));";
  test_pg("pmx_solve_group", expected);
  test_pg_count("pmx_solve_group", expected, 2);

  expected = "stan::math::assign(x, pmx_solve_group_bdf(twoCptModelODE_functor__(), 3, len, time, amt, rate, ii, evid, cmt, addl, ss, theta, biovar, tlag, 1.e-5, 1.e-5, 1000, pstream__));";
  test_pg("pmx_solve_group", expected);
  test_pg_count("pmx_solve_group", expected, 2);
}
