#ifndef STAN_MATH_PRIM_ERR_CHECK_NONNEGATIVE_HPP
#define STAN_MATH_PRIM_ERR_CHECK_NONNEGATIVE_HPP

#include <stan/math/prim/meta.hpp>
#include <stan/math/prim/err/throw_domain_error.hpp>
#include <stan/math/prim/err/throw_domain_error_vec.hpp>
#include <stan/math/prim/fun/get.hpp>
#include <stan/math/prim/fun/size.hpp>
#include <type_traits>

namespace stan {
namespace math {

namespace internal {
template <typename T_y, bool is_vec>
struct nonnegative {
  static void check(const char* function, const char* name, const T_y& y) {
    // have to use not is_unsigned. is_signed will be false for
    // floating point types that have no unsigned versions.
    if (!std::is_unsigned<T_y>::value && !(y >= 0)) {
      throw_domain_error(function, name, y, "is ", ", but must be >= 0!");
    }
  }
};

template <typename T_y>
struct nonnegative<T_y, true> {
  static void check(const char* function, const char* name, const T_y& y) {
    for (size_t n = 0; n < stan::math::size(y); n++) {
      if (!std::is_unsigned<typename value_type<T_y>::type>::value
          && !(stan::get(y, n) >= 0)) {
        throw_domain_error_vec(function, name, y, n, "is ",
                               ", but must be >= 0!");
      }
    }
  }
};
}  // namespace internal

/**
 * Check if <code>y</code> is non-negative.
 * This function is vectorized and will check each element of <code>y</code>.
 * @tparam T_y Type of y
 * @param function Function name (for error messages)
 * @param name Variable name (for error messages)
 * @param y Variable to check
 * @throw <code>domain_error</code> if y is negative or
 *   if any element of y is NaN.
 */
template <typename T_y>
inline void check_nonnegative(const char* function, const char* name,
                              const T_y& y) {
  internal::nonnegative<T_y, is_vector_like<T_y>::value>::check(function, name,
                                                                y);
}
}  // namespace math
}  // namespace stan
#endif
