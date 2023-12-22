#ifndef STAN_SERVICES_UTIL_VALIDATE_DIAG_INV_METRIC_HPP
#define STAN_SERVICES_UTIL_VALIDATE_DIAG_INV_METRIC_HPP

#include <stan/callbacks/logger.hpp>
#include <stan/math/prim.hpp>

namespace stan {
namespace services {
namespace util {

/**
 * Validate that diag inverse Euclidean metric is positive definite
 *
 * @param[in] inv_metric  inverse Euclidean metric
 * @param[in,out] logger Logger for messages
 * @throws std::domain_error if matrix is not positive definite
 */
inline void validate_diag_inv_metric(const Eigen::VectorXd& inv_metric,
                                     callbacks::logger& logger) {
  try {
    stan::math::check_finite("check_finite", "inv_metric", inv_metric);
    stan::math::check_positive("check_positive", "inv_metric", inv_metric);
  } catch (const std::domain_error& e) {
    logger.error("Inverse Euclidean metric not positive definite.");
    throw std::domain_error("Initialization failure");
  }
}

}  // namespace util
}  // namespace services
}  // namespace stan

#endif
