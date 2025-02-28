#include <test/unit/math/test_ad.hpp>

TEST(mathMixCore, operatorGreaterThanOrEqual) {
  auto f = [](const auto& x1, const auto& x2) { return x1 >= x2; };
  stan::test::expect_common_comparison(f);
}
