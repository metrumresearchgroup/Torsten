#ifndef STAN_MATH_TORSTEN_TEST_UTIL
#define STAN_MATH_TORSTEN_TEST_UTIL

#include <gtest/gtest.h>
#include <stan/math/rev/core.hpp>
#include <stan/math/prim/meta.hpp>
#include <stan/math/prim/meta/is_eigen_matrix_base.hpp>
#include <stan/math/rev/meta.hpp>
#include <stan/math/prim/fun/Eigen.hpp>
#include <stan/math/torsten/finite_diff_gradient.hpp>
#include <stan/math/torsten/torsten_def.hpp>
#include <stan/math/torsten/mpi/precomputed_gradients.hpp>
#include <stan/math/torsten/to_var.hpp>
#include <iostream>
#include <sstream>
#include <vector>
#include <string>

/**
 * Tests if any elementwise difference of the input matrices
 * of doubles is greater than DELTA. This uses the
 * EXPECT_NEAR macro from GTest.
 *
 * @param A first input matrix to compare
 * @param B second input matrix to compare
 * @param DELTA the maximum allowed difference
 */
#define EXPECT_MAT_VAL_NEAR(A, B, DELTA)                        \
  {                                                             \
    EXPECT_EQ(A.rows(), B.rows());                              \
    EXPECT_EQ(A.cols(), B.cols());                              \
    for (int i = 0; i < A.rows(); ++i) {                        \
      for (int j = 0; j < A.cols(); ++j) {                      \
        double a = stan::math::value_of(A(i, j));               \
        double b = stan::math::value_of(B(i, j));               \
        EXPECT_NEAR(a, b, DELTA) << "as entry (" << i << "," << j << ")";\
      }                                                         \
    }                                                           \
  }

/**
 * Tests if any elementwise difference of the input matrices
 * of doubles is greater than DELTA. This uses the
 * EXPECT_NEAR macro from GTest.
 *
 * @param A first input matrix to compare
 * @param B second input matrix to compare
 * @param DELTA the maximum allowed difference
 */
#define EXPECT_MAT_VAL_FLOAT_EQ(A, B)                           \
  {                                                             \
    EXPECT_EQ(A.rows(), B.rows());                              \
    EXPECT_EQ(A.cols(), B.cols());                              \
    for (int i = 0; i < A.rows(); ++i) {                        \
      for (int j = 0; j < A.cols(); ++j) {                      \
        double a = stan::math::value_of(A(i, j));               \
        double b = stan::math::value_of(B(i, j));               \
        EXPECT_FLOAT_EQ(a, b) << "as entry (" << i << "," << j << ")";\
      }                                                         \
    }                                                           \
  }

/**
 * Tests if any elementwise difference of the input 2d-arrays
 * of doubles is greater than DELTA. This uses the
 * EXPECT_FLOAT_EQ macro from GTest.
 *
 * @param A first input matrix to compare
 * @param B second input matrix to compare
 */
#define EXPECT_ARRAY2D_VAL_FLOAT_EQ(A, B)                               \
  {                                                                     \
  EXPECT_EQ(A.size(), B.size());                                        \
  for (int i = 0; i < A.size(); ++i) {                                  \
    for (int j = 0; j < A[i].size(); ++j) {                             \
      EXPECT_EQ(A[i].size(), B[i].size());                              \
      double a = stan::math::value_of(A[i][j]);                         \
      double b = stan::math::value_of(B[i][j]);                         \
      EXPECT_FLOAT_EQ(a, b) << "as entry (" << i << "," << j << ")";    \
    }                                                                   \
  }                                                                     \
  }

/**
 * Tests if any elementwise difference of the input matrices
 * of doubles is greater than DELTA. This uses the
 * EXPECT_NEAR macro from GTest.
 *
 * @param A first input matrix to compare
 * @param B second input matrix to compare
 * @param DELTA the maximum allowed difference
 */
#define EXPECT_ARRAY2D_VAL_NEAR(A, B, DELTA)                            \
  {                                                                     \
  EXPECT_EQ(A.size(), B.size());                                        \
  for (int i = 0; i < A.size(); ++i) {                                  \
    for (int j = 0; j < A[i].size(); ++j) {                             \
      EXPECT_EQ(A[i].size(), B[i].size());                              \
      double a = stan::math::value_of(A[i][j]);                         \
      double b = stan::math::value_of(B[i][j]);                         \
      EXPECT_NEAR(a, b, DELTA) << "as entry (" << i << "," << j << ")"; \
    }                                                                   \
  }                                                                     \
  }

/** 
 * For MAT vs 2D-array alike, check member value float equal
 * 
 * @param A Matrix
 * @param B 2D array
 * 
 */
#define EXPECT_MAT_ARRAY2D_VAL_FLOAT_EQ(A, B)                   \
  {                                                             \
    EXPECT_EQ(A.cols(), B.size());                              \
    for (int i = 0; i < A.cols(); ++i) {                        \
      for (int j = 0; j < A.rows(); ++j) {                   \
        double a = stan::math::value_of(A(j, i));               \
        double b = stan::math::value_of(B[i][j]);               \
        EXPECT_FLOAT_EQ(a, b) << "at (" << i << "," << j << ")";\
      }                                                         \
    }                                                           \
  }

#define EXPECT_ARRAY2D_ADJ_NEAR(A, B, P, NESTED, DELTA, MSG)                       \
  {                                                                     \
    EXPECT_EQ(A.size(), B.size());                                      \
    auto theta = stan::math::to_array_1d(P);                            \
    std::vector<double> ga, gb;                                         \
    for (auto i = 0; i < A.size(); ++i) {                               \
      for (auto j = 0; j < A[i].size(); ++j) {                          \
        NESTED.set_zero_all_adjoints();                                 \
        A[i][j].grad(theta, ga);                                            \
        NESTED.set_zero_all_adjoints();                                 \
        B[i][j].grad(theta, gb);                                            \
        for (auto k = 0; k < theta.size(); ++k) {                           \
          EXPECT_NEAR(ga[k], gb[k], DELTA) << k << "'th grad at " << i << ": " << MSG; \
        }                                                               \
      }                                                                 \
    }                                                                   \
  }

#define EXPECT_MAT_ARRAY2D_ADJ_NEAR(B, A, P, NESTED, DELTA, MSG)        \
  {                                                                     \
    EXPECT_EQ(A.size(), B.cols());                                      \
    std::vector<double> ga(P.size()), gb(P.size());                     \
    auto theta = stan::math::to_array_1d(P);                            \
    for (int i = 0; i < A.size(); ++i) {                                \
      for (int j = 0; j < A[i].size(); ++j) {                           \
        NESTED.set_zero_all_adjoints();                                 \
        A[i][j].grad(theta, ga);                                        \
        NESTED.set_zero_all_adjoints();                                 \
        B(j, i).grad(theta, gb);                                        \
        for (auto k = 0; k < P.size(); ++k) {                           \
          EXPECT_NEAR(ga[k], gb[k], DELTA) << k << "'th grad at (" << i << ", " << j<< ") of A and B: " << MSG; \
        }                                                               \
      }                                                                 \
    }                                                                   \
  }


/**
 * Tests if any elementwise difference of the input vectors
 * is greater than DELTA. This uses the
 * EXPECT_NEAR macro from GTest.
 *
 * @param A first input matrix to compare
 * @param B second input matrix to compare
 * @param DELTA the maximum allowed difference
 */
#define EXPECT_VEC_VAL_NEAR(A, B, DELTA)                        \
  {                                                             \
    EXPECT_EQ(A.size(), B.size());                              \
    for (int i = 0; i < A.size(); ++i) {                        \
      double a = stan::math::value_of(A[i]);                    \
      double b = stan::math::value_of(B[i]);                    \
      EXPECT_NEAR(a, b, DELTA) << "as entry " << i; \
    }                                                           \
  }

/**
 * Tests if any elementwise difference of the input vectors
 * is greater than DELTA. This uses the
 * EXPECT_NEAR macro from GTest.
 *
 * @param A first input matrix to compare
 * @param B second input matrix to compare
 * @param DELTA the maximum allowed difference
 */
#define EXPECT_VEC_VAL_FLOAT_EQ(A, B)                   \
  {                                                     \
    EXPECT_EQ(A.size(), B.size());                      \
    for (int i = 0; i < A.size(); ++i) {                \
      double a = stan::math::value_of(A[i]);            \
      double b = stan::math::value_of(B[i]);            \
      EXPECT_FLOAT_EQ(a, b) << "as entry " << i;        \
    }                                                   \
  }

/**
 * @param P parameter array of which grad is checked
 * @param NESTED nested_rev_autodiff where A, B and P live
 * @param DELTA the maximum allowed difference
 */
#define EXPECT_VEC_ADJ_NEAR(A, B, P, NESTED, DELTA, MSG)                \
  {                                                                     \
  EXPECT_EQ(A.size(), B.size());                                        \
  std::vector<double> ga(P.size()), gb(P.size());                       \
  auto theta = stan::math::to_array_1d(P);\
  for (int i = 0; i < A.size(); ++i) {                                  \
    nested.set_zero_all_adjoints();                                \
    A[i].grad(theta, ga);                                               \
    nested.set_zero_all_adjoints();                            \
    B[i].grad(theta, gb);                                               \
    for (auto k = 0; k < P.size(); ++k) {                               \
      EXPECT_NEAR(ga[k], gb[k], DELTA) << k << "'th grad at " << i << ": " << MSG; \
    }                                                                   \
  }                                                                     \
  }

/**
 * Tests if any elementwise difference of the input matrices'
 * element adjoint is greater than DELTA. This uses the
 * EXPECT_NEAR macro from GTest.
 *
 * @param A first input matrix to compare
 * @param B second input matrix to compare
 * @param P parameter array of which grad is checked
 * @param NESTED nested_rev_autodiff where A, B and P live
 * @param DELTA the maximum allowed difference
 */
#define EXPECT_MAT_ADJ_NEAR(A, B, P, NESTED, DELTA, MSG)                    \
  {                                                                     \
    EXPECT_EQ(A.rows(), B.rows());                                      \
    EXPECT_EQ(A.cols(), B.cols());                                      \
    std::vector<double> ga(P.size()), gb(P.size());                     \
    for (int i = 0; i < A.rows(); ++i) {                                \
      for (int j = 0; j < A.cols(); ++j) {                              \
        NESTED.set_zero_all_adjoints();                                 \
        A(i, j).grad(P, ga);                                            \
        NESTED.set_zero_all_adjoints();                                 \
        B(i, j).grad(P, gb);                                            \
        for (auto k = 0; k < P.size(); ++k) {                           \
          EXPECT_NEAR(ga[k], gb[k], DELTA) << k << "'th grad at (" << i << ", " << j<< ") of A and B: " << MSG; \
        }                                                               \
      }                                                                 \
    }                                                                   \
  }

/**
 * Tests if a function's autodiff gradient is near finite diff gradient
 * w.r.t given vector, only check positive compnents of the vector.
 *
 * @param FUNC function
 * @param P parameter vector
 * @param NESTED nested_rev_autodiff where P lives
 * @param H finite diff step size
 * @param TOL tolerance
 * @param MSG diagnostic message
 */
#define EXPECT_MAT_FUNC_POSITIVE_PARAM_NEAR_FD(FUN, P, NESTED, H, TOL, MSG) \
  {                                                                     \
    std::vector<stan::math::var> xvar = stan::math::to_var(P);          \
    std::vector<double> x_d = stan::math::value_of(P);                  \
    auto res1 = FUN(xvar);                                              \
    for (auto i = 0; i < x_d.size(); ++i) {                             \
      if (x_d[i] < h) continue;                                         \
      double x_i = x_d[i];                                              \
      x_d[i] += h;                                                      \
      auto res_ud = FUN(x_d);                                           \
      x_d[i] -= 2 * h;                                                  \
      auto res_ld = FUN(x_d);                                           \
      x_d[i] = x_i;                                                     \
      for (auto j = 0; j < res1.rows(); ++j) {                          \
        for (auto k = 0; k < res1.cols(); ++k) {                        \
          nested.set_zero_all_adjoints();                               \
          res1(j, k).grad();                                            \
          double res_fd = 0.5 * (stan::math::value_of(res_ud(j, k)) - stan::math::value_of(res_ld(j, k))) / h; \
          EXPECT_NEAR(xvar[i].adj(), res_fd, tol) << i << "'th param at (" << j << ", " << k << "):" << MSG; \
        }                                                               \
      }                                                                 \
    }                                                                   \
  }

/**
 * Tests if a function's autodiff gradient is near finite diff gradient
 * w.r.t given vector.
 *
 * @param FUNC function
 * @param P parameter vector
 * @param NESTED nested_rev_autodiff where P lives
 * @param H finite diff step size
 * @param TOL tolerance
 * @param MSG diagnostic message
 */
#define EXPECT_MAT_FUNC_PARAM_NEAR_FD(FUN, P, NESTED, H, TOL, MSG) \
  {                                                                     \
    std::vector<stan::math::var> xvar = stan::math::to_var(P);          \
    std::vector<double> x_d = stan::math::value_of(P);                  \
    auto res1 = FUN(xvar);                                              \
    for (auto i = 0; i < x_d.size(); ++i) {                             \
      double x_i = x_d[i];                                              \
      x_d[i] += h;                                                      \
      auto res_ud = FUN(x_d);                                           \
      x_d[i] -= 2 * h;                                                  \
      auto res_ld = FUN(x_d);                                           \
      x_d[i] = x_i;                                                     \
      for (auto j = 0; j < res1.rows(); ++j) {                          \
        for (auto k = 0; k < res1.cols(); ++k) {                        \
          nested.set_zero_all_adjoints();                               \
          res1(j, k).grad();                                            \
          double res_fd = 0.5 * (stan::math::value_of(res_ud(j, k)) - stan::math::value_of(res_ld(j, k))) / h; \
          EXPECT_NEAR(xvar[i].adj(), res_fd, tol) << i << "'th param at (" << j << ", " << k << "):" << MSG; \
        }                                                               \
      }                                                                 \
    }                                                                   \
  }

namespace torsten {
  namespace test {
    /*
     * Test @c std::vector<var> results between two results.
     * An example use would be to have the results coming from torsten
     * and stan, respectively, so ensure the soundness of
     * torsten results.
     *
     * @param y1 one result
     * @param y2 the other result to be compared against
     *              with, must of same shape and size as to @c pk_y
     */
    template<typename T1, typename T2,
             stan::require_all_stan_scalar_t<T1, T2>* = nullptr>
    inline void test_val(const T1& y1, const T2& y2) {
      EXPECT_FLOAT_EQ(stan::math::value_of(y1), stan::math::value_of(y2));
    }

    /*
     * Test @c std::vector<var> results between two results. 
     * An example use would be to have the results coming from torsten
     * and stan, respectively, so ensure the soundness of
     * torsten results.
     *
     * @param y1 one result
     * @param y2 the other result to be compared against
     *              with, must of same shape and size as to @c pk_y
     */
    template<typename T1, typename T2,
             stan::require_all_stan_scalar_t<T1, T2>* = nullptr>
    inline void test_val(const std::vector<T1>& y1,
                         const std::vector<T2>& y2) {
      EXPECT_EQ(y1.size(), y2.size());
      for (size_t i = 0; i < y1.size(); ++i) {
        EXPECT_FLOAT_EQ(stan::math::value_of(y1[i]), stan::math::value_of(y2[i]));
      }
    }

    /*
     * Test @c std::vector<var> results between two results. 
     * An example use would be to have the results coming from torsten
     * and stan, respectively, so ensure the soundness of
     * torsten results.
     *
     * @param y1 one result
     * @param y2 the other result to be compared against
     *              with, must of same shape and size as to @c pk_y
     */
    template<typename T1, typename T2,
             stan::require_all_stan_scalar_t<T1, T2>* = nullptr>
    void test_val(const std::vector<std::vector<T1>>& y1,
                  const std::vector<std::vector<T2>>& y2) {
      EXPECT_EQ(y1.size(), y2.size());
      for (size_t i = 0; i < y1.size(); ++i) {
        EXPECT_EQ(y1[i].size(), y2[i].size());
        for (size_t j = 0; j < y1[i].size(); ++j) {
          EXPECT_FLOAT_EQ(stan::math::value_of(y1[i][j]), stan::math::value_of(y2[i][j]));
        }
      }
    }

    /*
     * Test @c MatrixXd results between two results.
     * An example use would be to have the results coming from torsten
     * and stan, respectively, so ensure the soundness of
     * torsten results.
     *
     * @param y1 one result
     * @param y2 the other result to be compared against
     *              with, must of same shape and size as to @c pk_y
     * @param rtol relative tolerance
     * @param rtol absolute tolerance
     */
    template<typename T1, typename T2, stan::require_all_stan_scalar_t<T1, T2>* = nullptr>
    void test_val(const std::vector<std::vector<T1>>& y1,
                  const std::vector<std::vector<T2>>& y2,
                  double rtol, double atol) {
      using stan::math::value_of;
      EXPECT_EQ(y1.size(), y2.size());
      for (int i = 0; i < y1.size(); ++i) {
        EXPECT_EQ(y1[i].size(), y2[i].size());
        for (size_t j = 0; j < y1[i].size(); ++j) {
          double y1_ij = value_of(y1[i][j]);
          double y2_ij = value_of(y2[i][j]);
          if (abs(y1_ij) < 1e-5 && abs(y2_ij) < 1e-5) {
            EXPECT_NEAR(y1_ij, y2_ij, atol);
          } else {
            EXPECT_NEAR(y1_ij, y2_ij, std::max(abs(y1_ij), abs(y2_ij)) * rtol);
          }
        }
      }
    }

    /*
     * Test @c MatrixXd results between two results.
     * An example use would be to have the results coming from torsten
     * and stan, respectively, so ensure the soundness of
     * torsten results.
     *
     * @param y1 one result
     * @param y2 the other result to be compared against
     *              with, must of same shape and size as to @c pk_y
     */
    template<typename T1, typename T2,
             stan::require_all_eigen_matrix_base_t<T1, T2>* = nullptr>
    void test_val(const T1& y1, const T2& y2) {
      using stan::math::value_of;
      EXPECT_EQ(y1.rows(), y2.rows());
      EXPECT_EQ(y1.cols(), y2.cols());
      for (int i = 0; i < y1.size(); ++i) {
        EXPECT_FLOAT_EQ(value_of(y1(i)), value_of(y2(i)));
      }
    }

    /*
     * Test @c MatrixXd results between two results.
     * An example use would be to have the results coming from torsten
     * and stan, respectively, so ensure the soundness of
     * torsten results.
     *
     * @param y1 one result
     * @param y2 the other result to be compared against
     *              with, must of same shape and size as to @c pk_y
     * @param rtol relative tolerance
     * @param rtol absolute tolerance
     */
    template<typename T1, typename T2,
             stan::require_all_eigen_matrix_base_t<T1, T2>* = nullptr>
    void test_val(const T1& y1, const T2& y2,
                  double rtol, double atol) {
      using stan::math::value_of;
      EXPECT_EQ(y1.rows(), y2.rows());
      EXPECT_EQ(y1.cols(), y2.cols());
      for (int i = 0; i < y1.rows(); ++i) {
        for (int j = 0; j < y1.cols(); ++j) {
          double y1_i = value_of(y1(i, j));
          double y2_i = value_of(y2(i, j));
          if (abs(y1_i) < 1e-5 && abs(y2_i) < 1e-5) {
            EXPECT_NEAR(y1_i, y2_i, atol);
          } else {
            EXPECT_NEAR(y1_i, y2_i, std::max(abs(y1_i), abs(y2_i)) * rtol);
          }
        }
      }
    }

    /*
     * Test vector of vectors results against @c MatrixXd.
     *
     * @param y1 results in the form of vector of vectors
     * @param y2 @c MatrixXd results in column-major format
     */
    void test_val(const std::vector<std::vector<double> > & y1,
                  const Eigen::MatrixXd& y2) {
      EXPECT_EQ(y1.size(), y2.cols());
      for (size_t i = 0; i < y1.size(); ++i) {
        EXPECT_EQ(y1[i].size(), y2.rows());
        for (int j = 0; j < y2.rows(); ++j) {
          EXPECT_FLOAT_EQ(y1[i][j], y2(j, i));
        }
      }
    }

    /*
     * Test vector of vectors results against @c MatrixXd.
     *
     * @param y1 results in the form of vector of vectors
     * @param y2 @c MatrixXd results in column-major format
     */
    void test_val(const Eigen::MatrixXd& y2,
                  const std::vector<std::vector<double> > & y1) {
      test_val(y1, y2);
    }

    /*
     * Test @c VectorXd results between two results.
     * An example use would be to have the results coming from torsten
     * and stan, respectively, so ensure the soundness of
     * torsten results.
     *
     * @param y1 one result
     * @param y2 the other result to be compared against
     *              with, must of same shape and size as to @c pk_y
     */
    void test_val(const Eigen::VectorXd& y1, const Eigen::VectorXd& y2) {
      EXPECT_EQ(y1.size(), y2.size());
      for (int i = 0; i < y1.size(); ++i) {
        EXPECT_FLOAT_EQ(y1(i), y2(i));
      }
    }

    void test_val(const Eigen::VectorXd& y1, const Eigen::VectorXd& y2,
                  double rtol, double atol) {
      EXPECT_EQ(y1.size(), y2.size());
      for (int i = 0; i < y1.size(); ++i) {
        if (abs(y1(i)) < 1e-5 && abs(y2(i)) < 1e-5) {
          EXPECT_NEAR(y1(i), y2(i), atol);
        } else {
          EXPECT_NEAR(y1(i), y2(i), std::max(abs(y1(i)), abs(y2(i))) * rtol);
        }
      }
    }

    /*
     * Test @c VectorXd results between two results.
     * An example use would be to have the results coming from torsten
     * and stan, respectively, so ensure the soundness of
     * torsten results.
     *
     * @param y1 one result
     * @param y2 the other result to be compared against
     *              with, must of same shape and size as to @c pk_y
     */
    void test_val(const Eigen::VectorXd& y1, const std::vector<double>& y2) {
      EXPECT_EQ(y1.size(), y2.size());
      for (int i = 0; i < y1.size(); ++i) {
        EXPECT_FLOAT_EQ(y1(i), y2[i]);
      }
    }

    /*
     * Test @c VectorXd results between two results.
     * An example use would be to have the results coming from torsten
     * and stan, respectively, so ensure the soundness of
     * torsten results.
     *
     * @param y1 one result
     * @param y2 the other result to be compared against
     *              with, must of same shape and size as to @c pk_y
     */
    void test_val(const std::vector<double>& y1, const Eigen::VectorXd& y2) {
      test_val(y2, y1);
    }

    /*
     * Test @c std::vector<var> results between two results. 
     * An example use would be to have the results coming from torsten
     * and stan, respectively, so ensure the soundness of
     * torsten results.
     *
     * @param theta parameters regarding which the gradient
     *              would be taken and checked.
     * @param pk_y one result
     * @param stan_y the other result to be compared against
     *              with, must of same shape and size as to @c pk_y
     * @param fval_esp tolerance of values
     * @param sens_esp tolerance of gradients
     */
    void test_grad(std::vector<stan::math::var>& theta,
                   std::vector<std::vector<stan::math::var>>& pk_y,
                   std::vector<std::vector<stan::math::var>>& stan_y,
                   double fval_eps,
                   double sens_eps) {
      EXPECT_EQ(pk_y.size(), stan_y.size());
      for (size_t i = 0; i < pk_y.size(); ++i) { 
        EXPECT_EQ(pk_y[i].size(), stan_y[i].size());
      }

      for (size_t i = 0; i < pk_y.size(); ++i) {
        for (size_t j = 0; j < pk_y[i].size(); ++j) {
          EXPECT_NEAR(pk_y[i][j].val(), stan_y[i][j].val(), fval_eps);
        }
      }

      std::vector<double> g, g1;
      for (size_t i = 0; i < pk_y.size(); ++i) {
        for (size_t j = 0; j < pk_y[i].size(); ++j) {
          stan::math::set_zero_all_adjoints();
          pk_y[i][j].grad(theta, g);
          stan::math::set_zero_all_adjoints();
          stan_y[i][j].grad(theta, g1);
          for (size_t m = 0; m < theta.size(); ++m) {
            EXPECT_NEAR(g[m], g1[m], sens_eps);
          }
        }
      }
    }

    /*
     * Test @c std::vector<var> results between two results.
     * An example use would be to have the results coming from torsten
     * and stan, respectively, so ensure the soundness of
     * torsten results.
     *
     * @param theta1 parameters regarding which the gradient
     *              would be taken by @c y1 and checked.
     * @param theta2 parameters regarding which the gradient
     *              would be taken by @c y2 and checked.
     * @param pk_y one result
     * @param stan_y the other result to be compared against
     *              with, must of same shape and size as to @c pk_y
     * @param fval_esp tolerance of values
     * @param sens_esp tolerance of gradients
     */
    void test_grad(std::vector<stan::math::var>& theta1,
                   std::vector<stan::math::var>& theta2,
                   std::vector<std::vector<stan::math::var>>& y1,
                   std::vector<std::vector<stan::math::var>>& y2,
                   double fval_eps,
                   double sens_eps) {
      EXPECT_EQ(theta1.size(), theta2.size());
      EXPECT_EQ(y1.size(), y2.size());
      for (size_t i = 0; i < y1.size(); ++i) {
        EXPECT_EQ(y1[i].size(), y2[i].size());
      }

      for (size_t i = 0; i < y1.size(); ++i) {
        for (size_t j = 0; j < y1[i].size(); ++j) {
          EXPECT_NEAR(y1[i][j].val(), y2[i][j].val(), fval_eps);
        }
      }

      std::vector<double> g, g1;
      for (size_t i = 0; i < y1.size(); ++i) {
        for (size_t j = 0; j < y1[i].size(); ++j) {
          stan::math::set_zero_all_adjoints();
          y1[i][j].grad(theta1, g);
          stan::math::set_zero_all_adjoints();
          y2[i][j].grad(theta2, g1);
          for (size_t m = 0; m < theta1.size(); ++m) {
            EXPECT_NEAR(g[m], g1[m], sens_eps);
          }
        }
      }
    }

    /*
     * Test @c vector of @c var vectors against @c var matrix
     *
     * @param theta1 parameters regarding which the gradient
     *              would be taken by @c y1 and checked.
     * @param theta2 parameters regarding which the gradient
     *              would be taken by @c y2 and checked.
     * @param y1 result in form of vector of vectors
     * @param y2 result in form of matrix.
     * @param fval_esp tolerance of values
     * @param sens_esp tolerance of gradients
     */
    void test_grad(std::vector<stan::math::var>& theta1,
                   std::vector<stan::math::var>& theta2,
                   std::vector<std::vector<stan::math::var>>& y1,
                   Eigen::Matrix<stan::math::var, -1, -1>& y2,
                   double fval_eps,
                   double sens_eps) {
      EXPECT_EQ(theta1.size(), theta2.size());
      EXPECT_EQ(y1.size(), y2.cols());
      for (size_t i = 0; i < y1.size(); ++i) {
        EXPECT_EQ(y1[i].size(), y2.rows());
      }

      for (size_t i = 0; i < y1.size(); ++i) {
        for (size_t j = 0; j < y1[i].size(); ++j) {
          EXPECT_NEAR(y1[i][j].val(), y2(j, i).val(), fval_eps);
        }
      }

      std::vector<double> g, g1;
      for (size_t i = 0; i < y1.size(); ++i) {
        for (size_t j = 0; j < y1[i].size(); ++j) {
          stan::math::set_zero_all_adjoints();
          y1[i][j].grad(theta1, g);
          stan::math::set_zero_all_adjoints();
          y2(j, i).grad(theta2, g1);
          for (size_t m = 0; m < theta1.size(); ++m) {
            EXPECT_NEAR(g[m], g1[m], sens_eps);
          }
        }
      }
    }

    /*
     * Test @c vector of @c var vectors against a data matrix
     * containing both value and gradient.
     *
     * @param theta parameters regarding which the gradient
     *              would be taken by @c y1 and checked.
     * @param y1 result in form of vector of vectors
     * @param y2 result in form of dataa matrix, each column
     *              ordered as y1, dy1/dp1, dy1/dp2..., y2, dy2/dp1...
     * @param fval_esp tolerance of values
     * @param sens_esp tolerance of gradients
     */
    void test_grad(std::vector<stan::math::var>& theta,
                   std::vector<std::vector<stan::math::var>>& y1,
                   Eigen::MatrixXd& y2,
                   double fval_eps,
                   double sens_eps) {
      Eigen::Matrix<stan::math::var, -1, -1> y2_v = torsten::precomputed_gradients(y2, theta);
      test_grad(theta, theta, y1, y2_v, fval_eps, sens_eps);
    }

    /*
     * Test @c std::vector<var> results between two results. 
     * An example use would be to have the results coming from torsten
     * and stan, respectively, so ensure the soundness of
     * torsten results.
     *
     * @param theta parameters regarding which the gradient
     *              would be taken and checked.
     * @param pk_y one result
     * @param stan_y the other result to be compared against
     *              with, must of same shape and size as to @c pk_y
     * @param fval_esp tolerance of values
     * @param sens_esp tolerance of gradients
     */
    void test_grad(std::vector<stan::math::var>& theta,
                   Eigen::Matrix<stan::math::var, -1, -1>& pk_y,
                   Eigen::Matrix<stan::math::var, -1, -1>& stan_y,
                   double fval_eps,
                   double sens_eps) {
      EXPECT_EQ(pk_y.rows(), stan_y.rows());
      EXPECT_EQ(pk_y.cols(), stan_y.cols());

      for (int i = 0; i < pk_y.size(); ++i) {
        EXPECT_NEAR(pk_y(i).val(), stan_y(i).val(), fval_eps);
      }

      std::vector<double> g, g1;
      for (int i = 0; i < pk_y.size(); ++i) {
        stan::math::set_zero_all_adjoints();
        pk_y(i).grad(theta, g);
        stan::math::set_zero_all_adjoints();
        stan_y(i).grad(theta, g1);
        for (size_t m = 0; m < theta.size(); ++m) {
          EXPECT_NEAR(g[m], g1[m], sens_eps);
        }
      }
    }

    /*
     * Test @c std::vector<var> results between two results.
     * An example use would be to have the results coming from torsten
     * and stan, respectively, so ensure the soundness of
     * torsten results.
     *
     * @param theta parameters regarding which the gradient
     *              would be taken and checked.
     * @param y1 one result
     * @param y2 the other result to be compared against
     *              with, must of same shape and size as to @c y1
     * @param fval_esp tolerance of values
     * @param sens_esp tolerance of gradients
     */
    void test_grad(std::vector<stan::math::var>& theta,
                   stan::math::vector_v& y1,
                   stan::math::vector_v& y2,
                   double fval_eps,
                   double sens_eps) {
      EXPECT_EQ(y1.size(), y2.size());

      for (int i = 0; i < y1.size(); ++i) {
        EXPECT_NEAR(y1(i).val(), y2(i).val(), fval_eps);
      }

      std::vector<double> g, g1;
      for (int i = 0; i < y1.size(); ++i) {
        stan::math::set_zero_all_adjoints();
        y1(i).grad(theta, g);
        stan::math::set_zero_all_adjoints();
        y2(i).grad(theta, g1);
        for (size_t m = 0; m < theta.size(); ++m) {
          EXPECT_NEAR(g[m], g1[m], sens_eps);
        }
      }
    }

    /*
     * Test @c std::vector<var> results between two results.
     * An example use would be to have the results coming from torsten
     * and stan, respectively, so ensure the soundness of
     * torsten results.
     *
     * @param theta1 parameters regarding which the gradient
     *               would be taken by @c y1 and checked.
     * @param theta2 parameters regarding which the gradient
     *               would be taken by @c y2 and checked.
     * @param y1 one result
     * @param y2 the other result to be compared against
     *              with, must of same shape and size as to @c pk_y
     * @param fval_esp tolerance of values
     * @param sens_esp tolerance of gradients
     */
    void test_grad(std::vector<stan::math::var>& theta1,
                   std::vector<stan::math::var>& theta2,
                   stan::math::vector_v& y1,
                   stan::math::vector_v& y2,
                   double fval_eps,
                   double sens_eps) {
      EXPECT_EQ(y1.size(), y2.size());
      EXPECT_EQ(theta1.size(), theta2.size());

      for (int i = 0; i < y1.size(); ++i) {
        EXPECT_NEAR(y1(i).val(), y2(i).val(), fval_eps);
      }

      std::vector<double> g, g1;
      for (int i = 0; i < y1.size(); ++i) {
        stan::math::set_zero_all_adjoints();
        y1(i).grad(theta1, g);
        stan::math::set_zero_all_adjoints();
        y2(i).grad(theta2, g1);
        for (size_t m = 0; m < theta1.size(); ++m) {
          EXPECT_NEAR(g[m], g1[m], sens_eps);
        }
      }
    }

    /**
     * Test @c std::vector<var> results between two results.
     * An example use would be to have the results coming from torsten
     * and stan, respectively, so ensure the soundness of
     * torsten results.
     *
     * @param theta1 parameters regarding which the gradient
     *               would be taken by @c y1 and checked.
     * @param theta2 parameters regarding which the gradient
     *               would be taken by @c y2 and checked.
     * @param y1 one result
     * @param y2 the other result to be compared against
     *              with, must of same shape and size as to @c pk_y
     * @param fval_esp tolerance of values
     * @param sens_esp tolerance of gradients
     */
    void test_grad(std::vector<stan::math::var>& theta1,
                   torsten::PKRec<stan::math::var>& theta2,
                   stan::math::vector_v& y1,
                   stan::math::vector_v& y2,
                   double fval_eps,
                   double sens_eps) {
      // grad() only accepts std::vector
      std::vector<stan::math::var> theta(theta2.size());
      torsten::PKRec<stan::math::var>::Map(theta.data(), theta2.size()) = theta2;

      test_grad(theta1, theta, y1, y2, fval_eps, sens_eps);
    }

    /*
     * Test @c std::vector<var> results between two results.
     * An example use would be to have the results coming from torsten
     * and stan, respectively, so ensure the soundness of
     * torsten results.
     *
     * @param theta parameters regarding which the gradient
     *              would be taken and checked.
     * @param pk_y one result
     * @param stan_y the other result to be compared against
     *              with, must of same shape and size as to @c pk_y
     * @param fval_esp tolerance of values
     * @param sens_esp tolerance of gradients
     */
    void test_grad(std::vector<stan::math::var>& theta1,
                   std::vector<stan::math::var>& theta2,
                   Eigen::Matrix<stan::math::var, -1, -1>& y1,
                   Eigen::Matrix<stan::math::var, -1, -1>& y2,
                   double fval_eps,
                   double sens_eps) {
      EXPECT_EQ(theta1.size(), theta2.size());
      EXPECT_EQ(y1.rows(), y2.rows());
      EXPECT_EQ(y1.cols(), y2.cols());

      for (int i = 0; i < y1.size(); ++i) {
        EXPECT_NEAR(y1(i).val(), y2(i).val(), fval_eps);
      }

      std::vector<double> g, g1;
      for (int i = 0; i < y1.size(); ++i) {
        stan::math::set_zero_all_adjoints();
        y1(i).grad(theta1, g);
        stan::math::set_zero_all_adjoints();
        y2(i).grad(theta2, g1);
        for (size_t m = 0; m < theta1.size(); ++m) {
          EXPECT_NEAR(g[m], g1[m], sens_eps);
        }
      }
    }

    /*
     * Test @c std::vector<var> results between two results.
     * An example use would be to have the results coming from torsten
     * and stan, respectively, so ensure the soundness of
     * torsten results.
     *
     * @param theta parameters regarding which the gradient
     *              would be taken and checked.
     * @param pk_y one result
     * @param stan_y the other result to be compared against
     *              with, must of same shape and size as to @c pk_y
     */
    void test_grad(std::vector<stan::math::var>& theta1,
                   std::vector<stan::math::var>& theta2,
                   Eigen::Matrix<stan::math::var, -1, -1>& y1,
                   Eigen::Matrix<stan::math::var, -1, -1>& y2) {
      EXPECT_EQ(theta1.size(), theta2.size());
      EXPECT_EQ(y1.rows(), y2.rows());
      EXPECT_EQ(y1.cols(), y2.cols());

      for (int i = 0; i < y1.size(); ++i) {
        EXPECT_FLOAT_EQ(y1(i).val(), y2(i).val());
      }

      std::vector<double> g, g1;
      for (int i = 0; i < y1.size(); ++i) {
        stan::math::set_zero_all_adjoints();
        y1(i).grad(theta1, g);
        stan::math::set_zero_all_adjoints();
        y2(i).grad(theta2, g1);
        for (size_t m = 0; m < theta1.size(); ++m) {
          EXPECT_FLOAT_EQ(g[m], g1[m]);
        }
      }
    }

    /*
     * test gradients against finite difference
     *
     * Given a functor that takes a single parameter vector,
     * compare the gradients of the functor w.r.t. the
     * vector parameter.
     *
     * @tparam F1 functor type that return data, must takes a single vector
     * @tparam F2 functor type that return @c var, must takes a single vector
     * argument and returns a matrix.
     * @param f functor
     * @param theta parameter of the functor's function
     * @param h step size when eval finite difference
     * @param fval_esp tolerance of values
     * @param sens_esp tolerance of gradients
     */
    template<typename F1, typename F2>
    void test_grad(F1& f1, F2& f2,
                   std::vector<std::vector<double> >& theta,
                   double h,
                   double fval_eps,
                   double r_sens_eps,  double a_sens_eps) {
      std::vector<std::vector<stan::math::var> > theta_v(torsten::to_var(theta));
      std::vector<std::vector<double> > theta_1(theta.begin(), theta.end());
      std::vector<std::vector<double> > theta_2(theta.begin(), theta.end());

      const int m = theta.size();
      const int n = theta[0].size();
      Eigen::MatrixXd fd = f1(theta), fd_1, fd_2;
      Eigen::Matrix<stan::math::var, -1, -1> fv = f2(theta_v);

      EXPECT_EQ(fd.rows(), fv.rows());
      EXPECT_EQ(fd.cols(), fv.cols());
      for (int i = 0; i < fd.size(); ++i) {
        EXPECT_NEAR(fv(i).val(), fd(i), fval_eps);
      }

      std::vector<double> g;
      for (size_t i = 0; i < theta_v.size(); ++i) {
        for (size_t j = 0; j < theta_v[i].size(); ++j) {
          std::vector<stan::math::var> p{theta_v[i][j]};
          theta_1[i][j] -= h;
          theta_2[i][j] += h;
          fd_1 = f1(theta_1);
          fd_2 = f1(theta_2);
          theta_1[i][j] = theta[i][j];
          theta_2[i][j] = theta[i][j];
          for (int k = 0; k < fv.size(); ++k) {
            stan::math::set_zero_all_adjoints();
            fv(k).grad(p, g);
            double g_fd = (fd_2(k) - fd_1(k))/(2 * h);
            if (abs(g[0]) < 1e-4 || abs(g_fd) < 1e-4) {
              EXPECT_NEAR(g[0], g_fd, a_sens_eps);
            } else {
              EXPECT_NEAR(g[0], g_fd, r_sens_eps * std::max(abs(g[0]), abs(g_fd)));
            }
          }
        }
      }
    }

    /*
     * test gradients against finite difference
     *
     * Given a functor that takes a single parameter vector,
     * compare the gradients of the functor w.r.t. the
     * vector parameter.
     *
     * @tparam F1 functor type that return data, must takes a single vector
     * @tparam F2 functor type that return @c var, must takes a single vector
     * argument and returns a matrix.
     * @param f functor
     * @param theta parameter of the functor's function
     * @param h step size when eval finite difference
     * @param fval_esp tolerance of values
     * @param sens_esp tolerance of gradients
     */
    template<typename F1, typename F2>
    void test_grad(F1& f1, F2& f2,
                   std::vector<double>& theta,
                   double h,
                   double fval_eps,
                   double r_sens_eps,  double a_sens_eps) {
      using stan::math::value_of;
      std::vector<stan::math::var> theta_v(stan::math::to_var(theta));
      std::vector<double> theta_1(theta), theta_2(theta);

      auto fd = f1(theta);
      Eigen::Matrix<stan::math::var, -1, -1> fv = f2(theta_v);

      EXPECT_EQ(fd.rows(), fv.rows());
      EXPECT_EQ(fd.cols(), fv.cols());

      std::vector<double> g;
      for (int k = 0; k < fv.size(); ++k) {
        stan::math::set_zero_all_adjoints();
        fv(k).grad(theta_v, g);
        auto g_fd = torsten::finite_diff_gradient(f1, theta, k, h);
        auto fx = fd(k);
        EXPECT_NEAR(fv(k).val(), value_of(fx), fval_eps);
        for (size_t i = 0; i < g.size(); ++i) {
          if (abs(g[i]) < 1e-4 || abs(g_fd[i]) < 1e-4) {
            EXPECT_NEAR(g[i], value_of(g_fd[i]), a_sens_eps);
          } else {
            EXPECT_NEAR(g[i], value_of(g_fd[i]), r_sens_eps * std::max(abs(g[i]), abs(value_of(g_fd[i]))));
          }
        }
      }
    }

    /*
     * test gradients against finite difference
     *
     * Given a functor that takes a single parameter vector,
     * compare the gradients of the functor w.r.t. the
     * vector parameter.
     *
     * @tparam F1 functor type that return data, must takes a single vector
     * @tparam F2 functor type that return @c var, must takes a single vector
     * argument and returns a matrix.
     * @param f functor
     * @param theta parameter of the functor's function
     * @param h step size when eval finite difference
     * @param fval_esp tolerance of values
     * @param sens_esp tolerance of gradients
     */
    template<typename F1, typename F2>
    void test_grad(F1& f1, F2& f2,
                   std::vector<Eigen::MatrixXd>& theta,
                   double h,
                   double fval_eps,
                   double r_sens_eps,  double a_sens_eps) {
      using Eigen::MatrixXd;
      using stan::math::matrix_v;
      std::vector<matrix_v> theta_v(torsten::to_var(theta));
      std::vector<MatrixXd> theta_1(theta);
      std::vector<MatrixXd> theta_2(theta);

      Eigen::MatrixXd fd = f1(theta), fd_1, fd_2;
      matrix_v fv = f2(theta_v);

      EXPECT_EQ(fd.rows(), fv.rows());
      EXPECT_EQ(fd.cols(), fv.cols());
      for (int i = 0; i < fd.size(); ++i) {
        EXPECT_NEAR(fv(i).val(), fd(i), fval_eps);
      }

      std::vector<double> g;
      for (size_t i = 0; i < theta_v.size(); ++i) {
        for (int j = 0; j < theta_v[i].size(); ++j) {
          std::vector<stan::math::var> p{theta_v[i](j)};
          theta_1[i](j) -= h;
          theta_2[i](j) += h;
          fd_1 = f1(theta_1);
          fd_2 = f1(theta_2);
          theta_1[i](j) = theta[i](j);
          theta_2[i](j) = theta[i](j);
          for (int k = 0; k < fv.size(); ++k) {
            stan::math::set_zero_all_adjoints();
            fv(k).grad(p, g);
            double g_fd = (fd_2(k) - fd_1(k))/(2 * h);
            if (abs(g[0]) < 1e-4 || abs(g_fd) < 1e-4) {
              EXPECT_NEAR(g[0], g_fd, a_sens_eps);
            } else {
              EXPECT_NEAR(g[0], g_fd, r_sens_eps * std::max(abs(g[0]), abs(g_fd)));
            }
          }
        }
      }
    }

  } // namespace test
}

#endif
