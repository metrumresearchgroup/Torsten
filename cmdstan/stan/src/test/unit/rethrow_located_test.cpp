#include <stan/model/rethrow_located.hpp>
#include <gtest/gtest.h>

template <typename E, typename E2>
void test_rethrow_located_2() {
  try {
    try {
      throw E("foo");
    } catch (const std::exception& e) {
      stan::lang::rethrow_located(
          e, "(at line 5, file /Users/carp/temp2/bar.stan)");
    }
  } catch (const E2& e) {
    EXPECT_TRUE(std::string(e.what()).find_first_of("5") != std::string::npos);
    EXPECT_TRUE(std::string(e.what()).find_first_of("foo")
                != std::string::npos);
    return;
  } catch (...) {
    FAIL();
  }
  FAIL();
}

template <typename E>
void test_rethrow_located() {
  test_rethrow_located_2<E, E>();
}

template <typename E>
void test_rethrow_located_nullary(const std::string& original_type) {
  try {
    try {
      throw E();
    } catch (const std::exception& e) {
      stan::lang::rethrow_located(
          e, "(at line 5, file /Users/carp/temp2/bar.stan)");
    }
  } catch (const E& e) {
    EXPECT_TRUE(std::string(e.what()).find_first_of("5") != std::string::npos);
    EXPECT_TRUE(std::string(e.what()).find_first_of(original_type)
                != std::string::npos);
    return;
  } catch (...) {
    FAIL();
  }
  FAIL();
}

struct my_test_exception : public std::exception {
  const std::string what_;
  my_test_exception(const std::string& what) throw() : what_(what) {}
  ~my_test_exception() throw() {}
  const char* what() const throw() { return what_.c_str(); }
};

TEST(langRethrowLocated, allExpected) {
  test_rethrow_located_nullary<std::bad_alloc>("bad_alloc");
  test_rethrow_located_nullary<std::bad_cast>("bad_cast");
  test_rethrow_located_nullary<std::bad_exception>("bad_exception");
  test_rethrow_located_nullary<std::bad_typeid>("bad_typeid");

  test_rethrow_located<std::domain_error>();
  test_rethrow_located<std::invalid_argument>();
  test_rethrow_located<std::length_error>();
  test_rethrow_located<std::out_of_range>();
  test_rethrow_located<std::logic_error>();
  test_rethrow_located<std::overflow_error>();
  test_rethrow_located<std::range_error>();
  test_rethrow_located<std::underflow_error>();
  test_rethrow_located<std::runtime_error>();

  test_rethrow_located_nullary<std::exception>("std::exception");
  test_rethrow_located_2<my_test_exception, std::exception>();
}
