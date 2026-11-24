#ifndef STAN_CALLBACKS_TEE_WRITER_HPP
#define STAN_CALLBACKS_TEE_WRITER_HPP

#include <stan/math/prim/fun/Eigen.hpp>
#include <stan/math/prim/functor/apply.hpp>
#include <stan/math/prim/functor/for_each.hpp>
#include <stan/callbacks/writer.hpp>
#include <ostream>
#include <vector>
#include <string>

namespace stan::callbacks {

/**
 * `tee_writer` is an layer on top of a writer class that
 *  allows for multiple output streams to be written to.
 * @tparam Writers A parameter pack of types that inherit from `writer`
 */
template <typename... Writers>
class tee_writer {
 public:
  /**
   * Constructs a multi stream writer from a parameter pack of writers.
   * @param[in, out] args A parameter pack of writers
   */
  explicit tee_writer(Writers&... args) : output_(args...) {}

  tee_writer() = default;

  /**
   * @tparam T Any type accepted by a `writer` overload
   * @param[in] x A value to write to the output streams
   */
  template <typename T>
  void operator()(T&& x) {
    stan::math::for_each([&](auto&& output) { output(x); }, output_);
  }
  /**
   * Write a comment prefix to each writer
   */
  void operator()() {
    stan::math::for_each([](auto&& output) { output(); }, output_);
  }

  /**
   * Checks if all underlying writers are nonnull.
   */
  inline bool is_valid() const noexcept {
    return stan::math::apply(
        [](auto&&... output) { return (output.is_valid() && ...); }, output_);
  }

  /**
   * Get the tuple of underlying streams
   */
  inline auto& get_stream() noexcept { return output_; }

 private:
  // Output streams
  std::tuple<std::reference_wrapper<Writers>...> output_;
};

namespace internal {
template <typename T>
struct is_tee_writer : std::false_type {};

template <typename... Types>
struct is_tee_writer<tee_writer<Types...>> : std::true_type {};
}  // namespace internal

/**
 * Type trait that checks if a type is a `tee_writer`
 * @tparam T A type to check
 */
template <typename T>
struct is_tee_writer : internal::is_tee_writer<std::decay_t<T>> {};

/**
 * Helper variable template to check if a type is a `tee_writer`
 */
template <typename T>
inline constexpr bool is_tee_writer_v = is_tee_writer<T>::value;

}  // namespace stan::callbacks

#endif
