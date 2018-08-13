#ifndef TEST_MATH_MATRIX_EXPECT_NEAR_MATRIX_EQ_HPP
#define TEST_MATH_MATRIX_EXPECT_NEAR_MATRIX_EQ_HPP

#include <gtest/gtest.h>

void expect_near_matrix_eq(const
                           Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic>& a,
                           const
                           Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic>& b,
                           const double& rel_err,
                           const double& abs_err = 1e-15) {
  EXPECT_EQ(a.rows(), b.rows());
  EXPECT_EQ(a.cols(), b.cols());
  double err;
  for (int i = 0; i < a.rows(); ++i)
    for (int j = 0; j < a.cols(); ++j) {
     err = std::max(std::abs(rel_err * a(i, j)), abs_err);
     EXPECT_NEAR(a(i, j), b(i, j), err);
    }
}

#endif
