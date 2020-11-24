#ifndef STAN_MATH_OPENCL_KERNELS_CHECK_NAN_HPP
#define STAN_MATH_OPENCL_KERNELS_CHECK_NAN_HPP
#ifdef STAN_OPENCL

#include <stan/math/opencl/kernel_cl.hpp>
#include <stan/math/opencl/buffer_types.hpp>
#include <string>

namespace stan {
namespace math {
namespace opencl_kernels {
// \cond
static const std::string is_nan_kernel_code = STRINGIFY(
    // \endcond
    /** \ingroup opencl_kernels
     * Check if the <code>matrix_cl</code> has NaN values
     *
     * @param[in] A The matrix to check.
     * @param rows The number of rows in matrix A.
     * @param cols The number of columns in matrix A.
     * @param[out] flag the flag to be written to if any diagonal is zero.
     * @note Code is a <code>const char*</code> held in
     * <code>is_nan_kernel_code.</code>
     *  Kernel for stan/math/opencl/err/check_nan.hpp.
     *  This kernel uses the helper macros available in helpers.cl.
     */
    __kernel void is_nan(__global double *A, __global int *flag,
                         unsigned int rows, unsigned int cols) {
      const int i = get_global_id(0);
      const int j = get_global_id(1);
      if (i < rows && j < cols) {
        if (isnan(A(i, j))) {
          flag[0] = 1;
        }
      }
    }
    // \cond
);
// \endcond

/** \ingroup opencl_kernels
 * See the docs for \link kernels/check_nan.hpp is_nan() \endlink
 */
const kernel_cl<in_buffer, out_buffer, int, int> check_nan(
    "is_nan", {indexing_helpers, is_nan_kernel_code});

}  // namespace opencl_kernels
}  // namespace math
}  // namespace stan
#endif
#endif
