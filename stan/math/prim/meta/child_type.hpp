#ifndef STAN_MATH_PRIM_META_CHILD_TYPE_HPP
#define STAN_MATH_PRIM_META_CHILD_TYPE_HPP

#include <type_traits>

namespace stan {
namespace math {

/** \ingroup type_trait
 * Primary template class for metaprogram to compute child type of
 * T.
 *
 * See <code>test/unit/math/meta/child_type_test.cpp</code> for
 * intended usage.
 *
 * @tparam T type of container.
 */
template <typename T>
struct child_type {
  using type = double;
};

/** \ingroup type_trait
 * Specialization for template classes / structs.
 *
 * See <code>test/unit/math/meta/child_type_test.cpp</code> for
 * intended usage.
 *
 * @tparam T_struct type of parent.
 * @tparam T_child type of child type.
 */
template <template <typename> class T_struct, typename T_child>
struct child_type<T_struct<T_child>> {
  using type = T_child;
};

}  // namespace math
}  // namespace stan

#endif
