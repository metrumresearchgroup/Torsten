#ifndef STAN_MATH_TORSTEN_TEST_UTIL
#define STAN_MATH_TORSTEN_TEST_UTIL

#include <stan/math/prim/mat/fun/Eigen.hpp>
#include <stan/math/rev/mat.hpp>
#include <gtest/gtest.h>
#include <test/unit/math/prim/arr/functor/harmonic_oscillator.hpp>
#include <test/unit/math/prim/arr/functor/lorenz.hpp>
#include <stan/math/torsten/mpi/precomputed_gradients.hpp>
#include <stan/math/torsten/to_var.hpp>
#include <iostream>
#include <sstream>
#include <vector>
#include <string>

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
    template<typename T1, typename T2>
    void test_val(const std::vector<std::vector<T1>>& y1,
                  const std::vector<std::vector<T2>>& y2) {
      using stan::math::value_of;
      EXPECT_EQ(y1.size(), y2.size());
      for (size_t i = 0; i < y1.size(); ++i) {
        EXPECT_EQ(y1[i].size(), y2[i].size());
        for (size_t j = 0; j < y1[i].size(); ++j) {
          EXPECT_FLOAT_EQ(value_of(y1[i][j]), value_of(y2[i][j]));
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
    template<typename T1, typename T2>
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
    template<typename T1, typename T2>
    void test_val(const Eigen::Matrix<T1, -1, -1>& y1,
                  const Eigen::Matrix<T2, -1, -1>& y2) {
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
    template<typename T1, typename T2>
    void test_val(const Eigen::Matrix<T1, -1, -1>& y1,
                  const Eigen::Matrix<T2, -1, -1>& y2,
                  double rtol, double atol) {
      using stan::math::value_of;
      EXPECT_EQ(y1.rows(), y2.rows());
      EXPECT_EQ(y1.cols(), y2.cols());
      for (int i = 0; i < y1.size(); ++i) {
        double y1_i = value_of(y1(i));
        double y2_i = value_of(y2(i));
        if (abs(y1_i) < 1e-5 && abs(y2_i) < 1e-5) {
          EXPECT_NEAR(y1_i, y2_i, atol);
        } else {
          EXPECT_NEAR(y1_i, y2_i, std::max(abs(y1_i), abs(y2_i)) * rtol);
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
    template<typename T1, typename T2>
    void test_val(const Eigen::Matrix<T1, -1, 1>& y1,
                  const Eigen::Matrix<T2, -1, 1>& y2) {
      using stan::math::value_of;
      EXPECT_EQ(y1.size(), y2.size());
      for (int i = 0; i < y1.size(); ++i) {
        EXPECT_FLOAT_EQ(value_of(y1(i)), value_of(y2(i)));
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
                   Eigen::Matrix<stan::math::var, 1, -1>& theta2,
                   stan::math::vector_v& y1,
                   stan::math::vector_v& y2,
                   double fval_eps,
                   double sens_eps) {
      // grad() only accepts std::vector
      std::vector<stan::math::var> theta(theta2.size());
      Eigen::Matrix<stan::math::var, 1, -1>::Map(theta.data(), theta2.size()) = theta2;

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
      std::vector<stan::math::var> theta_v(stan::math::to_var(theta));
      std::vector<double> theta_1(theta), theta_2(theta);

      Eigen::MatrixXd fd = f1(theta), fd_1, fd_2;
      Eigen::Matrix<stan::math::var, -1, -1> fv = f2(theta_v);

      EXPECT_EQ(fd.rows(), fv.rows());
      EXPECT_EQ(fd.cols(), fv.cols());
      for (int i = 0; i < fd.size(); ++i) {
        EXPECT_NEAR(fv(i).val(), fd(i), fval_eps);
      }

      std::vector<double> g;
      for (size_t i = 0; i < theta_v.size(); ++i) {
        std::vector<stan::math::var> p{theta_v[i]};
        theta_1[i] -= h;
        theta_2[i] += h;
        fd_1 = f1(theta_1);
        fd_2 = f1(theta_2);
        theta_1[i] = theta[i];
        theta_2[i] = theta[i];
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
}   // namespace torsten

/*
 * Macro to test overloaded torsten functions when @c theta,
 * @c biovar and @c tlag that can be constatnt or time-dependent.
 */
#define TORSTEN_CPT_PARAM_OVERLOAD_TEST(FUN, TIME, AMT, RATE, II, EVID, CMT, ADDL, SS, THETA, BIOVAR, TLAG, \
                                    EPS_VAL, EPS_GRAD)                                                      \
  {                                                                                                         \
  std::vector<std::vector<stan::math::var> > theta_v(1, stan::math::to_var(THETA[0]));                      \
  std::vector<std::vector<stan::math::var> > biovar_v(1, stan::math::to_var(BIOVAR[0]));                    \
  std::vector<std::vector<stan::math::var> > tlag_v(1, stan::math::to_var(tlag[0]));                        \
  {                                                                                                         \
      auto x0 = FUN(TIME, AMT, RATE, II, EVID, CMT, ADDL, SS, theta_v,    BIOVAR,    TLAG);                 \
      auto x1 = FUN(TIME, AMT, RATE, II, EVID, CMT, ADDL, SS, theta_v[0], BIOVAR,    TLAG);                 \
      auto x2 = FUN(TIME, AMT, RATE, II, EVID, CMT, ADDL, SS, theta_v[0], BIOVAR[0], TLAG);                 \
      auto x3 = FUN(TIME, AMT, RATE, II, EVID, CMT, ADDL, SS, theta_v[0], BIOVAR[0], TLAG[0]);              \
      auto x4 = FUN(TIME, AMT, RATE, II, EVID, CMT, ADDL, SS, theta_v[0], BIOVAR,    TLAG[0]);              \
      auto x5 = FUN(TIME, AMT, RATE, II, EVID, CMT, ADDL, SS, theta_v,    BIOVAR[0], TLAG);                 \
      auto x6 = FUN(TIME, AMT, RATE, II, EVID, CMT, ADDL, SS, theta_v,    BIOVAR[0], TLAG[0]);              \
      auto x7 = FUN(TIME, AMT, RATE, II, EVID, CMT, ADDL, SS, theta_v,    BIOVAR,    TLAG[0]);              \
      torsten::test::test_grad(theta_v[0],  x0, x1, EPS_VAL, EPS_GRAD);                                     \
      torsten::test::test_grad(theta_v[0],  x0, x2, EPS_VAL, EPS_GRAD);                                     \
      torsten::test::test_grad(theta_v[0],  x0, x3, EPS_VAL, EPS_GRAD);                                     \
      torsten::test::test_grad(theta_v[0],  x0, x4, EPS_VAL, EPS_GRAD);                                     \
      torsten::test::test_grad(theta_v[0],  x0, x5, EPS_VAL, EPS_GRAD);                                     \
      torsten::test::test_grad(theta_v[0],  x0, x6, EPS_VAL, EPS_GRAD);                                     \
      torsten::test::test_grad(theta_v[0],  x0, x7, EPS_VAL, EPS_GRAD);                                     \
  }                                                                                                         \
  {                                                                                                         \
      auto x0 = FUN(TIME, AMT, RATE, II, EVID, CMT, ADDL, SS, THETA,    biovar_v,    TLAG);                 \
      auto x1 = FUN(TIME, AMT, RATE, II, EVID, CMT, ADDL, SS, THETA[0], biovar_v,    TLAG);                 \
      auto x2 = FUN(TIME, AMT, RATE, II, EVID, CMT, ADDL, SS, THETA[0], biovar_v[0], TLAG);                 \
      auto x3 = FUN(TIME, AMT, RATE, II, EVID, CMT, ADDL, SS, THETA[0], biovar_v[0], TLAG[0]);              \
      auto x4 = FUN(TIME, AMT, RATE, II, EVID, CMT, ADDL, SS, THETA[0], biovar_v,    TLAG[0]);              \
      auto x5 = FUN(TIME, AMT, RATE, II, EVID, CMT, ADDL, SS, THETA,    biovar_v[0], TLAG);                 \
      auto x6 = FUN(TIME, AMT, RATE, II, EVID, CMT, ADDL, SS, THETA,    biovar_v[0], TLAG[0]);              \
      auto x7 = FUN(TIME, AMT, RATE, II, EVID, CMT, ADDL, SS, THETA,    biovar_v,    TLAG[0]);              \
      torsten::test::test_grad(biovar_v[0],  x0, x1, EPS_VAL, EPS_GRAD);                                    \
      torsten::test::test_grad(biovar_v[0],  x0, x2, EPS_VAL, EPS_GRAD);                                    \
      torsten::test::test_grad(biovar_v[0],  x0, x3, EPS_VAL, EPS_GRAD);                                    \
      torsten::test::test_grad(biovar_v[0],  x0, x4, EPS_VAL, EPS_GRAD);                                    \
      torsten::test::test_grad(biovar_v[0],  x0, x5, EPS_VAL, EPS_GRAD);                                    \
      torsten::test::test_grad(biovar_v[0],  x0, x6, EPS_VAL, EPS_GRAD);                                    \
      torsten::test::test_grad(biovar_v[0],  x0, x7, EPS_VAL, EPS_GRAD);                                    \
  }                                                                                                         \
  {                                                                                                         \
      auto x0 = FUN(TIME, AMT, RATE, II, EVID, CMT, ADDL, SS, THETA,    BIOVAR,    tlag_v);                 \
      auto x1 = FUN(TIME, AMT, RATE, II, EVID, CMT, ADDL, SS, THETA[0], BIOVAR,    tlag_v);                 \
      auto x2 = FUN(TIME, AMT, RATE, II, EVID, CMT, ADDL, SS, THETA[0], BIOVAR[0], tlag_v);                 \
      auto x3 = FUN(TIME, AMT, RATE, II, EVID, CMT, ADDL, SS, THETA[0], BIOVAR[0], tlag_v[0]);              \
      auto x4 = FUN(TIME, AMT, RATE, II, EVID, CMT, ADDL, SS, THETA[0], BIOVAR,    tlag_v[0]);              \
      auto x5 = FUN(TIME, AMT, RATE, II, EVID, CMT, ADDL, SS, THETA,    BIOVAR[0], tlag_v);                 \
      auto x6 = FUN(TIME, AMT, RATE, II, EVID, CMT, ADDL, SS, THETA,    BIOVAR[0], tlag_v[0]);              \
      auto x7 = FUN(TIME, AMT, RATE, II, EVID, CMT, ADDL, SS, THETA,    BIOVAR,    tlag_v[0]);              \
      torsten::test::test_grad(tlag_v[0],  x0, x1, EPS_VAL, EPS_GRAD);                                      \
      torsten::test::test_grad(tlag_v[0],  x0, x2, EPS_VAL, EPS_GRAD);                                      \
      torsten::test::test_grad(tlag_v[0],  x0, x3, EPS_VAL, EPS_GRAD);                                      \
      torsten::test::test_grad(tlag_v[0],  x0, x4, EPS_VAL, EPS_GRAD);                                      \
      torsten::test::test_grad(tlag_v[0],  x0, x5, EPS_VAL, EPS_GRAD);                                      \
      torsten::test::test_grad(tlag_v[0],  x0, x6, EPS_VAL, EPS_GRAD);                                      \
      torsten::test::test_grad(tlag_v[0],  x0, x7, EPS_VAL, EPS_GRAD);                                      \
  }                                                                                                         \
  {                                                                                                         \
      auto x0 = FUN(TIME, AMT, RATE, II, EVID, CMT, ADDL, SS, theta_v,    biovar_v,    TLAG);               \
      auto x1 = FUN(TIME, AMT, RATE, II, EVID, CMT, ADDL, SS, theta_v[0], biovar_v,    TLAG);               \
      auto x2 = FUN(TIME, AMT, RATE, II, EVID, CMT, ADDL, SS, theta_v[0], biovar_v[0], TLAG);               \
      auto x3 = FUN(TIME, AMT, RATE, II, EVID, CMT, ADDL, SS, theta_v[0], biovar_v[0], TLAG[0]);            \
      auto x4 = FUN(TIME, AMT, RATE, II, EVID, CMT, ADDL, SS, theta_v[0], biovar_v,    TLAG[0]);            \
      auto x5 = FUN(TIME, AMT, RATE, II, EVID, CMT, ADDL, SS, theta_v,    biovar_v[0], TLAG);               \
      auto x6 = FUN(TIME, AMT, RATE, II, EVID, CMT, ADDL, SS, theta_v,    biovar_v[0], TLAG[0]);            \
      auto x7 = FUN(TIME, AMT, RATE, II, EVID, CMT, ADDL, SS, theta_v,    biovar_v,    TLAG[0]);            \
      torsten::test::test_grad(theta_v[0],  x0, x1, EPS_VAL, EPS_GRAD);                                     \
      torsten::test::test_grad(theta_v[0],  x0, x2, EPS_VAL, EPS_GRAD);                                     \
      torsten::test::test_grad(theta_v[0],  x0, x3, EPS_VAL, EPS_GRAD);                                     \
      torsten::test::test_grad(theta_v[0],  x0, x4, EPS_VAL, EPS_GRAD);                                     \
      torsten::test::test_grad(theta_v[0],  x0, x5, EPS_VAL, EPS_GRAD);                                     \
      torsten::test::test_grad(theta_v[0],  x0, x6, EPS_VAL, EPS_GRAD);                                     \
      torsten::test::test_grad(theta_v[0],  x0, x7, EPS_VAL, EPS_GRAD);                                     \
      torsten::test::test_grad(biovar_v[0],  x0, x1, EPS_VAL, EPS_GRAD);                                    \
      torsten::test::test_grad(biovar_v[0],  x0, x2, EPS_VAL, EPS_GRAD);                                    \
      torsten::test::test_grad(biovar_v[0],  x0, x3, EPS_VAL, EPS_GRAD);                                    \
      torsten::test::test_grad(biovar_v[0],  x0, x4, EPS_VAL, EPS_GRAD);                                    \
      torsten::test::test_grad(biovar_v[0],  x0, x5, EPS_VAL, EPS_GRAD);                                    \
      torsten::test::test_grad(biovar_v[0],  x0, x6, EPS_VAL, EPS_GRAD);                                    \
      torsten::test::test_grad(biovar_v[0],  x0, x7, EPS_VAL, EPS_GRAD);                                    \
  }                                                                                                         \
  {                                                                                                         \
      auto x0 = FUN(TIME, AMT, RATE, II, EVID, CMT, ADDL, SS, theta_v,    BIOVAR,    tlag_v);               \
      auto x1 = FUN(TIME, AMT, RATE, II, EVID, CMT, ADDL, SS, theta_v[0], BIOVAR,    tlag_v);               \
      auto x2 = FUN(TIME, AMT, RATE, II, EVID, CMT, ADDL, SS, theta_v[0], BIOVAR[0], tlag_v);               \
      auto x3 = FUN(TIME, AMT, RATE, II, EVID, CMT, ADDL, SS, theta_v[0], BIOVAR[0], tlag_v[0]);            \
      auto x4 = FUN(TIME, AMT, RATE, II, EVID, CMT, ADDL, SS, theta_v[0], BIOVAR,    tlag_v[0]);            \
      auto x5 = FUN(TIME, AMT, RATE, II, EVID, CMT, ADDL, SS, theta_v,    BIOVAR[0], tlag_v);               \
      auto x6 = FUN(TIME, AMT, RATE, II, EVID, CMT, ADDL, SS, theta_v,    BIOVAR[0], tlag_v[0]);            \
      auto x7 = FUN(TIME, AMT, RATE, II, EVID, CMT, ADDL, SS, theta_v,    BIOVAR,    tlag_v[0]);            \
      torsten::test::test_grad(theta_v[0],  x0, x1, EPS_VAL, EPS_GRAD);                                     \
      torsten::test::test_grad(theta_v[0],  x0, x2, EPS_VAL, EPS_GRAD);                                     \
      torsten::test::test_grad(theta_v[0],  x0, x3, EPS_VAL, EPS_GRAD);                                     \
      torsten::test::test_grad(theta_v[0],  x0, x4, EPS_VAL, EPS_GRAD);                                     \
      torsten::test::test_grad(theta_v[0],  x0, x5, EPS_VAL, EPS_GRAD);                                     \
      torsten::test::test_grad(theta_v[0],  x0, x6, EPS_VAL, EPS_GRAD);                                     \
      torsten::test::test_grad(theta_v[0],  x0, x7, EPS_VAL, EPS_GRAD);                                     \
      torsten::test::test_grad(tlag_v[0],  x0, x1, EPS_VAL, EPS_GRAD);                                      \
      torsten::test::test_grad(tlag_v[0],  x0, x2, EPS_VAL, EPS_GRAD);                                      \
      torsten::test::test_grad(tlag_v[0],  x0, x3, EPS_VAL, EPS_GRAD);                                      \
      torsten::test::test_grad(tlag_v[0],  x0, x4, EPS_VAL, EPS_GRAD);                                      \
      torsten::test::test_grad(tlag_v[0],  x0, x5, EPS_VAL, EPS_GRAD);                                      \
      torsten::test::test_grad(tlag_v[0],  x0, x6, EPS_VAL, EPS_GRAD);                                      \
      torsten::test::test_grad(tlag_v[0],  x0, x7, EPS_VAL, EPS_GRAD);                                      \
  }                                                                                                         \
  {                                                                                                         \
      auto x0 = FUN(TIME, AMT, RATE, II, EVID, CMT, ADDL, SS, THETA,    biovar_v,    tlag_v);               \
      auto x1 = FUN(TIME, AMT, RATE, II, EVID, CMT, ADDL, SS, THETA[0], biovar_v,    tlag_v);               \
      auto x2 = FUN(TIME, AMT, RATE, II, EVID, CMT, ADDL, SS, THETA[0], biovar_v[0], tlag_v);               \
      auto x3 = FUN(TIME, AMT, RATE, II, EVID, CMT, ADDL, SS, THETA[0], biovar_v[0], tlag_v[0]);            \
      auto x4 = FUN(TIME, AMT, RATE, II, EVID, CMT, ADDL, SS, THETA[0], biovar_v,    tlag_v[0]);            \
      auto x5 = FUN(TIME, AMT, RATE, II, EVID, CMT, ADDL, SS, THETA,    biovar_v[0], tlag_v);               \
      auto x6 = FUN(TIME, AMT, RATE, II, EVID, CMT, ADDL, SS, THETA,    biovar_v[0], tlag_v[0]);            \
      auto x7 = FUN(TIME, AMT, RATE, II, EVID, CMT, ADDL, SS, THETA,    biovar_v,    tlag_v[0]);            \
      torsten::test::test_grad(biovar_v[0],  x0, x1, EPS_VAL, EPS_GRAD);                                    \
      torsten::test::test_grad(biovar_v[0],  x0, x2, EPS_VAL, EPS_GRAD);                                    \
      torsten::test::test_grad(biovar_v[0],  x0, x3, EPS_VAL, EPS_GRAD);                                    \
      torsten::test::test_grad(biovar_v[0],  x0, x4, EPS_VAL, EPS_GRAD);                                    \
      torsten::test::test_grad(biovar_v[0],  x0, x5, EPS_VAL, EPS_GRAD);                                    \
      torsten::test::test_grad(biovar_v[0],  x0, x6, EPS_VAL, EPS_GRAD);                                    \
      torsten::test::test_grad(biovar_v[0],  x0, x7, EPS_VAL, EPS_GRAD);                                    \
      torsten::test::test_grad(tlag_v[0],  x0, x1, EPS_VAL, EPS_GRAD);                                      \
      torsten::test::test_grad(tlag_v[0],  x0, x2, EPS_VAL, EPS_GRAD);                                      \
      torsten::test::test_grad(tlag_v[0],  x0, x3, EPS_VAL, EPS_GRAD);                                      \
      torsten::test::test_grad(tlag_v[0],  x0, x4, EPS_VAL, EPS_GRAD);                                      \
      torsten::test::test_grad(tlag_v[0],  x0, x5, EPS_VAL, EPS_GRAD);                                      \
      torsten::test::test_grad(tlag_v[0],  x0, x6, EPS_VAL, EPS_GRAD);                                      \
      torsten::test::test_grad(tlag_v[0],  x0, x7, EPS_VAL, EPS_GRAD);                                      \
  }                                                                                                         \
  {                                                                                                         \
      auto x0 = FUN(TIME, AMT, RATE, II, EVID, CMT, ADDL, SS, theta_v,    biovar_v,    tlag_v);             \
      auto x1 = FUN(TIME, AMT, RATE, II, EVID, CMT, ADDL, SS, theta_v[0], biovar_v,    tlag_v);             \
      auto x2 = FUN(TIME, AMT, RATE, II, EVID, CMT, ADDL, SS, theta_v[0], biovar_v[0], tlag_v);             \
      auto x3 = FUN(TIME, AMT, RATE, II, EVID, CMT, ADDL, SS, theta_v[0], biovar_v[0], tlag_v[0]);          \
      auto x4 = FUN(TIME, AMT, RATE, II, EVID, CMT, ADDL, SS, theta_v[0], biovar_v,    tlag_v[0]);          \
      auto x5 = FUN(TIME, AMT, RATE, II, EVID, CMT, ADDL, SS, theta_v,    biovar_v[0], tlag_v);             \
      auto x6 = FUN(TIME, AMT, RATE, II, EVID, CMT, ADDL, SS, theta_v,    biovar_v[0], tlag_v[0]);          \
      auto x7 = FUN(TIME, AMT, RATE, II, EVID, CMT, ADDL, SS, theta_v,    biovar_v,    tlag_v[0]);          \
      torsten::test::test_grad(theta_v[0],  x0, x1, EPS_VAL, EPS_GRAD);                                     \
      torsten::test::test_grad(theta_v[0],  x0, x2, EPS_VAL, EPS_GRAD);                                     \
      torsten::test::test_grad(theta_v[0],  x0, x3, EPS_VAL, EPS_GRAD);                                     \
      torsten::test::test_grad(theta_v[0],  x0, x4, EPS_VAL, EPS_GRAD);                                     \
      torsten::test::test_grad(theta_v[0],  x0, x5, EPS_VAL, EPS_GRAD);                                     \
      torsten::test::test_grad(theta_v[0],  x0, x6, EPS_VAL, EPS_GRAD);                                     \
      torsten::test::test_grad(theta_v[0],  x0, x7, EPS_VAL, EPS_GRAD);                                     \
      torsten::test::test_grad(biovar_v[0],  x0, x1, EPS_VAL, EPS_GRAD);                                    \
      torsten::test::test_grad(biovar_v[0],  x0, x2, EPS_VAL, EPS_GRAD);                                    \
      torsten::test::test_grad(biovar_v[0],  x0, x3, EPS_VAL, EPS_GRAD);                                    \
      torsten::test::test_grad(biovar_v[0],  x0, x4, EPS_VAL, EPS_GRAD);                                    \
      torsten::test::test_grad(biovar_v[0],  x0, x5, EPS_VAL, EPS_GRAD);                                    \
      torsten::test::test_grad(biovar_v[0],  x0, x6, EPS_VAL, EPS_GRAD);                                    \
      torsten::test::test_grad(biovar_v[0],  x0, x7, EPS_VAL, EPS_GRAD);                                    \
      torsten::test::test_grad(tlag_v[0],  x0, x1, EPS_VAL, EPS_GRAD);                                      \
      torsten::test::test_grad(tlag_v[0],  x0, x2, EPS_VAL, EPS_GRAD);                                      \
      torsten::test::test_grad(tlag_v[0],  x0, x3, EPS_VAL, EPS_GRAD);                                      \
      torsten::test::test_grad(tlag_v[0],  x0, x4, EPS_VAL, EPS_GRAD);                                      \
      torsten::test::test_grad(tlag_v[0],  x0, x5, EPS_VAL, EPS_GRAD);                                      \
      torsten::test::test_grad(tlag_v[0],  x0, x6, EPS_VAL, EPS_GRAD);                                      \
      torsten::test::test_grad(tlag_v[0],  x0, x7, EPS_VAL, EPS_GRAD);                                      \
  }                                                                                                         \
  }

/*
 * Macro to test overloaded torsten functions when @c theta,
 * @c biovar and @c tlag that can be constatnt or time-dependent.
 */
#define TORSTEN_ODE_PARAM_OVERLOAD_TEST(FUN, F, NCMT, TIME, AMT, RATE, II, EVID, CMT, ADDL, SS, THETA, BIOVAR, TLAG, \
                                    EPS_VAL, EPS_GRAD)                                                               \
  {                                                                                                                  \
  std::vector<std::vector<stan::math::var> > theta_v(1, stan::math::to_var(THETA[0]));                               \
  std::vector<std::vector<stan::math::var> > biovar_v(1, stan::math::to_var(BIOVAR[0]));                             \
  std::vector<std::vector<stan::math::var> > tlag_v(1, stan::math::to_var(tlag[0]));                                 \
  {                                                                                                                  \
      auto x0 = FUN(F, NCMT, TIME, AMT, RATE, II, EVID, CMT, ADDL, SS, theta_v,    BIOVAR,    TLAG);                 \
      auto x1 = FUN(F, NCMT, TIME, AMT, RATE, II, EVID, CMT, ADDL, SS, theta_v[0], BIOVAR,    TLAG);                 \
      auto x2 = FUN(F, NCMT, TIME, AMT, RATE, II, EVID, CMT, ADDL, SS, theta_v[0], BIOVAR[0], TLAG);                 \
      auto x3 = FUN(F, NCMT, TIME, AMT, RATE, II, EVID, CMT, ADDL, SS, theta_v[0], BIOVAR[0], TLAG[0]);              \
      auto x4 = FUN(F, NCMT, TIME, AMT, RATE, II, EVID, CMT, ADDL, SS, theta_v[0], BIOVAR,    TLAG[0]);              \
      auto x5 = FUN(F, NCMT, TIME, AMT, RATE, II, EVID, CMT, ADDL, SS, theta_v,    BIOVAR[0], TLAG);                 \
      auto x6 = FUN(F, NCMT, TIME, AMT, RATE, II, EVID, CMT, ADDL, SS, theta_v,    BIOVAR[0], TLAG[0]);              \
      auto x7 = FUN(F, NCMT, TIME, AMT, RATE, II, EVID, CMT, ADDL, SS, theta_v,    BIOVAR,    TLAG[0]);              \
      torsten::test::test_grad(theta_v[0],  x0, x1, EPS_VAL, EPS_GRAD);                                              \
      torsten::test::test_grad(theta_v[0],  x0, x2, EPS_VAL, EPS_GRAD);                                              \
      torsten::test::test_grad(theta_v[0],  x0, x3, EPS_VAL, EPS_GRAD);                                              \
      torsten::test::test_grad(theta_v[0],  x0, x4, EPS_VAL, EPS_GRAD);                                              \
      torsten::test::test_grad(theta_v[0],  x0, x5, EPS_VAL, EPS_GRAD);                                              \
      torsten::test::test_grad(theta_v[0],  x0, x6, EPS_VAL, EPS_GRAD);                                              \
      torsten::test::test_grad(theta_v[0],  x0, x7, EPS_VAL, EPS_GRAD);                                              \
  }                                                                                                                  \
  {                                                                                                                  \
      auto x0 = FUN(F, NCMT, TIME, AMT, RATE, II, EVID, CMT, ADDL, SS, THETA,    biovar_v,    TLAG);                 \
      auto x1 = FUN(F, NCMT, TIME, AMT, RATE, II, EVID, CMT, ADDL, SS, THETA[0], biovar_v,    TLAG);                 \
      auto x2 = FUN(F, NCMT, TIME, AMT, RATE, II, EVID, CMT, ADDL, SS, THETA[0], biovar_v[0], TLAG);                 \
      auto x3 = FUN(F, NCMT, TIME, AMT, RATE, II, EVID, CMT, ADDL, SS, THETA[0], biovar_v[0], TLAG[0]);              \
      auto x4 = FUN(F, NCMT, TIME, AMT, RATE, II, EVID, CMT, ADDL, SS, THETA[0], biovar_v,    TLAG[0]);              \
      auto x5 = FUN(F, NCMT, TIME, AMT, RATE, II, EVID, CMT, ADDL, SS, THETA,    biovar_v[0], TLAG);                 \
      auto x6 = FUN(F, NCMT, TIME, AMT, RATE, II, EVID, CMT, ADDL, SS, THETA,    biovar_v[0], TLAG[0]);              \
      auto x7 = FUN(F, NCMT, TIME, AMT, RATE, II, EVID, CMT, ADDL, SS, THETA,    biovar_v,    TLAG[0]);              \
      torsten::test::test_grad(biovar_v[0],  x0, x1, EPS_VAL, EPS_GRAD);                                             \
      torsten::test::test_grad(biovar_v[0],  x0, x2, EPS_VAL, EPS_GRAD);                                             \
      torsten::test::test_grad(biovar_v[0],  x0, x3, EPS_VAL, EPS_GRAD);                                             \
      torsten::test::test_grad(biovar_v[0],  x0, x4, EPS_VAL, EPS_GRAD);                                             \
      torsten::test::test_grad(biovar_v[0],  x0, x5, EPS_VAL, EPS_GRAD);                                             \
      torsten::test::test_grad(biovar_v[0],  x0, x6, EPS_VAL, EPS_GRAD);                                             \
      torsten::test::test_grad(biovar_v[0],  x0, x7, EPS_VAL, EPS_GRAD);                                             \
  }                                                                                                                  \
  {                                                                                                                  \
      auto x0 = FUN(F, NCMT, TIME, AMT, RATE, II, EVID, CMT, ADDL, SS, THETA,    BIOVAR,    tlag_v);                 \
      auto x1 = FUN(F, NCMT, TIME, AMT, RATE, II, EVID, CMT, ADDL, SS, THETA[0], BIOVAR,    tlag_v);                 \
      auto x2 = FUN(F, NCMT, TIME, AMT, RATE, II, EVID, CMT, ADDL, SS, THETA[0], BIOVAR[0], tlag_v);                 \
      auto x3 = FUN(F, NCMT, TIME, AMT, RATE, II, EVID, CMT, ADDL, SS, THETA[0], BIOVAR[0], tlag_v[0]);              \
      auto x4 = FUN(F, NCMT, TIME, AMT, RATE, II, EVID, CMT, ADDL, SS, THETA[0], BIOVAR,    tlag_v[0]);              \
      auto x5 = FUN(F, NCMT, TIME, AMT, RATE, II, EVID, CMT, ADDL, SS, THETA,    BIOVAR[0], tlag_v);                 \
      auto x6 = FUN(F, NCMT, TIME, AMT, RATE, II, EVID, CMT, ADDL, SS, THETA,    BIOVAR[0], tlag_v[0]);              \
      auto x7 = FUN(F, NCMT, TIME, AMT, RATE, II, EVID, CMT, ADDL, SS, THETA,    BIOVAR,    tlag_v[0]);              \
      torsten::test::test_grad(tlag_v[0],  x0, x1, EPS_VAL, EPS_GRAD);                                               \
      torsten::test::test_grad(tlag_v[0],  x0, x2, EPS_VAL, EPS_GRAD);                                               \
      torsten::test::test_grad(tlag_v[0],  x0, x3, EPS_VAL, EPS_GRAD);                                               \
      torsten::test::test_grad(tlag_v[0],  x0, x4, EPS_VAL, EPS_GRAD);                                               \
      torsten::test::test_grad(tlag_v[0],  x0, x5, EPS_VAL, EPS_GRAD);                                               \
      torsten::test::test_grad(tlag_v[0],  x0, x6, EPS_VAL, EPS_GRAD);                                               \
      torsten::test::test_grad(tlag_v[0],  x0, x7, EPS_VAL, EPS_GRAD);                                               \
  }                                                                                                                  \
  {                                                                                                                  \
      auto x0 = FUN(F, NCMT, TIME, AMT, RATE, II, EVID, CMT, ADDL, SS, theta_v,    biovar_v,    TLAG);               \
      auto x1 = FUN(F, NCMT, TIME, AMT, RATE, II, EVID, CMT, ADDL, SS, theta_v[0], biovar_v,    TLAG);               \
      auto x2 = FUN(F, NCMT, TIME, AMT, RATE, II, EVID, CMT, ADDL, SS, theta_v[0], biovar_v[0], TLAG);               \
      auto x3 = FUN(F, NCMT, TIME, AMT, RATE, II, EVID, CMT, ADDL, SS, theta_v[0], biovar_v[0], TLAG[0]);            \
      auto x4 = FUN(F, NCMT, TIME, AMT, RATE, II, EVID, CMT, ADDL, SS, theta_v[0], biovar_v,    TLAG[0]);            \
      auto x5 = FUN(F, NCMT, TIME, AMT, RATE, II, EVID, CMT, ADDL, SS, theta_v,    biovar_v[0], TLAG);               \
      auto x6 = FUN(F, NCMT, TIME, AMT, RATE, II, EVID, CMT, ADDL, SS, theta_v,    biovar_v[0], TLAG[0]);            \
      auto x7 = FUN(F, NCMT, TIME, AMT, RATE, II, EVID, CMT, ADDL, SS, theta_v,    biovar_v,    TLAG[0]);            \
      torsten::test::test_grad(theta_v[0],  x0, x1, EPS_VAL, EPS_GRAD);                                              \
      torsten::test::test_grad(theta_v[0],  x0, x2, EPS_VAL, EPS_GRAD);                                              \
      torsten::test::test_grad(theta_v[0],  x0, x3, EPS_VAL, EPS_GRAD);                                              \
      torsten::test::test_grad(theta_v[0],  x0, x4, EPS_VAL, EPS_GRAD);                                              \
      torsten::test::test_grad(theta_v[0],  x0, x5, EPS_VAL, EPS_GRAD);                                              \
      torsten::test::test_grad(theta_v[0],  x0, x6, EPS_VAL, EPS_GRAD);                                              \
      torsten::test::test_grad(theta_v[0],  x0, x7, EPS_VAL, EPS_GRAD);                                              \
      torsten::test::test_grad(biovar_v[0],  x0, x1, EPS_VAL, EPS_GRAD);                                             \
      torsten::test::test_grad(biovar_v[0],  x0, x2, EPS_VAL, EPS_GRAD);                                             \
      torsten::test::test_grad(biovar_v[0],  x0, x3, EPS_VAL, EPS_GRAD);                                             \
      torsten::test::test_grad(biovar_v[0],  x0, x4, EPS_VAL, EPS_GRAD);                                             \
      torsten::test::test_grad(biovar_v[0],  x0, x5, EPS_VAL, EPS_GRAD);                                             \
      torsten::test::test_grad(biovar_v[0],  x0, x6, EPS_VAL, EPS_GRAD);                                             \
      torsten::test::test_grad(biovar_v[0],  x0, x7, EPS_VAL, EPS_GRAD);                                             \
  }                                                                                                                  \
  {                                                                                                                  \
      auto x0 = FUN(F, NCMT, TIME, AMT, RATE, II, EVID, CMT, ADDL, SS, theta_v,    BIOVAR,    tlag_v);               \
      auto x1 = FUN(F, NCMT, TIME, AMT, RATE, II, EVID, CMT, ADDL, SS, theta_v[0], BIOVAR,    tlag_v);               \
      auto x2 = FUN(F, NCMT, TIME, AMT, RATE, II, EVID, CMT, ADDL, SS, theta_v[0], BIOVAR[0], tlag_v);               \
      auto x3 = FUN(F, NCMT, TIME, AMT, RATE, II, EVID, CMT, ADDL, SS, theta_v[0], BIOVAR[0], tlag_v[0]);            \
      auto x4 = FUN(F, NCMT, TIME, AMT, RATE, II, EVID, CMT, ADDL, SS, theta_v[0], BIOVAR,    tlag_v[0]);            \
      auto x5 = FUN(F, NCMT, TIME, AMT, RATE, II, EVID, CMT, ADDL, SS, theta_v,    BIOVAR[0], tlag_v);               \
      auto x6 = FUN(F, NCMT, TIME, AMT, RATE, II, EVID, CMT, ADDL, SS, theta_v,    BIOVAR[0], tlag_v[0]);            \
      auto x7 = FUN(F, NCMT, TIME, AMT, RATE, II, EVID, CMT, ADDL, SS, theta_v,    BIOVAR,    tlag_v[0]);            \
      torsten::test::test_grad(theta_v[0],  x0, x1, EPS_VAL, EPS_GRAD);                                              \
      torsten::test::test_grad(theta_v[0],  x0, x2, EPS_VAL, EPS_GRAD);                                              \
      torsten::test::test_grad(theta_v[0],  x0, x3, EPS_VAL, EPS_GRAD);                                              \
      torsten::test::test_grad(theta_v[0],  x0, x4, EPS_VAL, EPS_GRAD);                                              \
      torsten::test::test_grad(theta_v[0],  x0, x5, EPS_VAL, EPS_GRAD);                                              \
      torsten::test::test_grad(theta_v[0],  x0, x6, EPS_VAL, EPS_GRAD);                                              \
      torsten::test::test_grad(theta_v[0],  x0, x7, EPS_VAL, EPS_GRAD);                                              \
      torsten::test::test_grad(tlag_v[0],  x0, x1, EPS_VAL, EPS_GRAD);                                               \
      torsten::test::test_grad(tlag_v[0],  x0, x2, EPS_VAL, EPS_GRAD);                                               \
      torsten::test::test_grad(tlag_v[0],  x0, x3, EPS_VAL, EPS_GRAD);                                               \
      torsten::test::test_grad(tlag_v[0],  x0, x4, EPS_VAL, EPS_GRAD);                                               \
      torsten::test::test_grad(tlag_v[0],  x0, x5, EPS_VAL, EPS_GRAD);                                               \
      torsten::test::test_grad(tlag_v[0],  x0, x6, EPS_VAL, EPS_GRAD);                                               \
      torsten::test::test_grad(tlag_v[0],  x0, x7, EPS_VAL, EPS_GRAD);                                               \
  }                                                                                                                  \
  {                                                                                                                  \
      auto x0 = FUN(F, NCMT, TIME, AMT, RATE, II, EVID, CMT, ADDL, SS, THETA,    biovar_v,    tlag_v);               \
      auto x1 = FUN(F, NCMT, TIME, AMT, RATE, II, EVID, CMT, ADDL, SS, THETA[0], biovar_v,    tlag_v);               \
      auto x2 = FUN(F, NCMT, TIME, AMT, RATE, II, EVID, CMT, ADDL, SS, THETA[0], biovar_v[0], tlag_v);               \
      auto x3 = FUN(F, NCMT, TIME, AMT, RATE, II, EVID, CMT, ADDL, SS, THETA[0], biovar_v[0], tlag_v[0]);            \
      auto x4 = FUN(F, NCMT, TIME, AMT, RATE, II, EVID, CMT, ADDL, SS, THETA[0], biovar_v,    tlag_v[0]);            \
      auto x5 = FUN(F, NCMT, TIME, AMT, RATE, II, EVID, CMT, ADDL, SS, THETA,    biovar_v[0], tlag_v);               \
      auto x6 = FUN(F, NCMT, TIME, AMT, RATE, II, EVID, CMT, ADDL, SS, THETA,    biovar_v[0], tlag_v[0]);            \
      auto x7 = FUN(F, NCMT, TIME, AMT, RATE, II, EVID, CMT, ADDL, SS, THETA,    biovar_v,    tlag_v[0]);            \
      torsten::test::test_grad(biovar_v[0],  x0, x1, EPS_VAL, EPS_GRAD);                                             \
      torsten::test::test_grad(biovar_v[0],  x0, x2, EPS_VAL, EPS_GRAD);                                             \
      torsten::test::test_grad(biovar_v[0],  x0, x3, EPS_VAL, EPS_GRAD);                                             \
      torsten::test::test_grad(biovar_v[0],  x0, x4, EPS_VAL, EPS_GRAD);                                             \
      torsten::test::test_grad(biovar_v[0],  x0, x5, EPS_VAL, EPS_GRAD);                                             \
      torsten::test::test_grad(biovar_v[0],  x0, x6, EPS_VAL, EPS_GRAD);                                             \
      torsten::test::test_grad(biovar_v[0],  x0, x7, EPS_VAL, EPS_GRAD);                                             \
      torsten::test::test_grad(tlag_v[0],  x0, x1, EPS_VAL, EPS_GRAD);                                               \
      torsten::test::test_grad(tlag_v[0],  x0, x2, EPS_VAL, EPS_GRAD);                                               \
      torsten::test::test_grad(tlag_v[0],  x0, x3, EPS_VAL, EPS_GRAD);                                               \
      torsten::test::test_grad(tlag_v[0],  x0, x4, EPS_VAL, EPS_GRAD);                                               \
      torsten::test::test_grad(tlag_v[0],  x0, x5, EPS_VAL, EPS_GRAD);                                               \
      torsten::test::test_grad(tlag_v[0],  x0, x6, EPS_VAL, EPS_GRAD);                                               \
      torsten::test::test_grad(tlag_v[0],  x0, x7, EPS_VAL, EPS_GRAD);                                               \
  }                                                                                                                  \
  {                                                                                                                  \
      auto x0 = FUN(F, NCMT, TIME, AMT, RATE, II, EVID, CMT, ADDL, SS, theta_v,    biovar_v,    tlag_v);             \
      auto x1 = FUN(F, NCMT, TIME, AMT, RATE, II, EVID, CMT, ADDL, SS, theta_v[0], biovar_v,    tlag_v);             \
      auto x2 = FUN(F, NCMT, TIME, AMT, RATE, II, EVID, CMT, ADDL, SS, theta_v[0], biovar_v[0], tlag_v);             \
      auto x3 = FUN(F, NCMT, TIME, AMT, RATE, II, EVID, CMT, ADDL, SS, theta_v[0], biovar_v[0], tlag_v[0]);          \
      auto x4 = FUN(F, NCMT, TIME, AMT, RATE, II, EVID, CMT, ADDL, SS, theta_v[0], biovar_v,    tlag_v[0]);          \
      auto x5 = FUN(F, NCMT, TIME, AMT, RATE, II, EVID, CMT, ADDL, SS, theta_v,    biovar_v[0], tlag_v);             \
      auto x6 = FUN(F, NCMT, TIME, AMT, RATE, II, EVID, CMT, ADDL, SS, theta_v,    biovar_v[0], tlag_v[0]);          \
      auto x7 = FUN(F, NCMT, TIME, AMT, RATE, II, EVID, CMT, ADDL, SS, theta_v,    biovar_v,    tlag_v[0]);          \
      torsten::test::test_grad(theta_v[0],  x0, x1, EPS_VAL, EPS_GRAD);                                              \
      torsten::test::test_grad(theta_v[0],  x0, x2, EPS_VAL, EPS_GRAD);                                              \
      torsten::test::test_grad(theta_v[0],  x0, x3, EPS_VAL, EPS_GRAD);                                              \
      torsten::test::test_grad(theta_v[0],  x0, x4, EPS_VAL, EPS_GRAD);                                              \
      torsten::test::test_grad(theta_v[0],  x0, x5, EPS_VAL, EPS_GRAD);                                              \
      torsten::test::test_grad(theta_v[0],  x0, x6, EPS_VAL, EPS_GRAD);                                              \
      torsten::test::test_grad(theta_v[0],  x0, x7, EPS_VAL, EPS_GRAD);                                              \
      torsten::test::test_grad(biovar_v[0],  x0, x1, EPS_VAL, EPS_GRAD);                                             \
      torsten::test::test_grad(biovar_v[0],  x0, x2, EPS_VAL, EPS_GRAD);                                             \
      torsten::test::test_grad(biovar_v[0],  x0, x3, EPS_VAL, EPS_GRAD);                                             \
      torsten::test::test_grad(biovar_v[0],  x0, x4, EPS_VAL, EPS_GRAD);                                             \
      torsten::test::test_grad(biovar_v[0],  x0, x5, EPS_VAL, EPS_GRAD);                                             \
      torsten::test::test_grad(biovar_v[0],  x0, x6, EPS_VAL, EPS_GRAD);                                             \
      torsten::test::test_grad(biovar_v[0],  x0, x7, EPS_VAL, EPS_GRAD);                                             \
      torsten::test::test_grad(tlag_v[0],  x0, x1, EPS_VAL, EPS_GRAD);                                               \
      torsten::test::test_grad(tlag_v[0],  x0, x2, EPS_VAL, EPS_GRAD);                                               \
      torsten::test::test_grad(tlag_v[0],  x0, x3, EPS_VAL, EPS_GRAD);                                               \
      torsten::test::test_grad(tlag_v[0],  x0, x4, EPS_VAL, EPS_GRAD);                                               \
      torsten::test::test_grad(tlag_v[0],  x0, x5, EPS_VAL, EPS_GRAD);                                               \
      torsten::test::test_grad(tlag_v[0],  x0, x6, EPS_VAL, EPS_GRAD);                                               \
      torsten::test::test_grad(tlag_v[0],  x0, x7, EPS_VAL, EPS_GRAD);                                               \
  }                                                                                                                  \
  }


#define TORSTEN_CPT_GRAD_THETA_TEST(FUN, TIME, AMT, RATE, II, EVID, CMT, ADDL, SS, THETA, BIOVAR, TLAG, \
                                    H, EPS_VAL, EPS_RTOL, EPS_ATOL)                                     \
  {                                                                                                     \
    auto f1 = [&] (std::vector<std::vector<double> >& x) {                                              \
      return FUN(TIME, AMT, RATE, II, EVID, CMT, ADDL, SS, x, BIOVAR, TLAG);                            \
    };                                                                                                  \
    auto f2 = [&] (std::vector<std::vector<stan::math::var> >& x) {                                     \
      return FUN(TIME, AMT, RATE, II, EVID, CMT, ADDL, SS, x, BIOVAR, TLAG);                            \
    };                                                                                                  \
    torsten::test::test_grad(f1, f2, THETA, H, EPS_VAL, EPS_RTOL, EPS_ATOL);                            \
  }

#define TORSTEN_LIN_GRAD_THETA_TEST(FUN, TIME, AMT, RATE, II, EVID, CMT, ADDL, SS, THETA, BIOVAR, TLAG, \
                                    H, EPS_VAL, EPS_RTOL, EPS_ATOL)                                     \
  {                                                                                                     \
    auto f1 = [&] (std::vector<Eigen::MatrixXd>& x) {                                                   \
      return FUN(TIME, AMT, RATE, II, EVID, CMT, ADDL, SS, x, BIOVAR, TLAG);                            \
    };                                                                                                  \
    auto f2 = [&] (std::vector<stan::math::matrix_v>& x) {                                              \
      return FUN(TIME, AMT, RATE, II, EVID, CMT, ADDL, SS, x, BIOVAR, TLAG);                            \
    };                                                                                                  \
    torsten::test::test_grad(f1, f2, THETA, H, EPS_VAL, EPS_RTOL, EPS_ATOL);                            \
  }

#define TORSTEN_CPT_GRAD_BIOVAR_TEST(FUN, TIME, AMT, RATE, II, EVID, CMT, ADDL, SS, THETA, BIOVAR, TLAG, \
                                    H, EPS_VAL, EPS_RTOL, EPS_ATOL)                                      \
  {                                                                                                      \
    auto f1 = [&] (std::vector<std::vector<double> >& x) {                                               \
      return FUN(TIME, AMT, RATE, II, EVID, CMT, ADDL, SS, THETA, x, TLAG);                              \
    };                                                                                                   \
    auto f2 = [&] (std::vector<std::vector<stan::math::var> >& x) {                                      \
      return FUN(TIME, AMT, RATE, II, EVID, CMT, ADDL, SS, THETA, x, TLAG);                              \
    };                                                                                                   \
    torsten::test::test_grad(f1, f2, BIOVAR, H, EPS_VAL, EPS_RTOL, EPS_ATOL);                            \
  }

#define TORSTEN_CPT_GRAD_TLAG_TEST(FUN, TIME, AMT, RATE, II, EVID, CMT, ADDL, SS, THETA, BIOVAR, TLAG,   \
                                    H, EPS_VAL, EPS_RTOL, EPS_ATOL)                                      \
  {                                                                                                      \
    auto f1 = [&] (std::vector<std::vector<double> >& x) {                                               \
      return FUN(TIME, AMT, RATE, II, EVID, CMT, ADDL, SS, THETA, BIOVAR, x);                            \
    };                                                                                                   \
    auto f2 = [&] (std::vector<std::vector<stan::math::var> >& x) {                                      \
      return FUN(TIME, AMT, RATE, II, EVID, CMT, ADDL, SS, THETA, BIOVAR, x);                            \
    };                                                                                                   \
    torsten::test::test_grad(f1, f2, TLAG, H, EPS_VAL, EPS_RTOL, EPS_ATOL);                              \
  }

#define TORSTEN_CPT_GRAD_RATE_TEST(FUN, TIME, AMT, RATE, II, EVID, CMT, ADDL, SS, THETA, BIOVAR, TLAG,   \
                                    H, EPS_VAL, EPS_RTOL, EPS_ATOL)                                      \
  {                                                                                                      \
    auto f1 = [&] (std::vector<double>& x) {                                                             \
      return FUN(TIME, AMT, x, II, EVID, CMT, ADDL, SS, THETA, BIOVAR, TLAG);                            \
    };                                                                                                   \
    auto f2 = [&] (std::vector<stan::math::var>& x) {                                                    \
      return FUN(TIME, AMT, x, II, EVID, CMT, ADDL, SS, THETA, BIOVAR, TLAG);                            \
    };                                                                                                   \
    torsten::test::test_grad(f1, f2, RATE, H, EPS_VAL, EPS_RTOL, EPS_ATOL);                              \
  }

#define TORSTEN_CPT_ODE_GRAD_TEST(FUN_CPT, FUN_ODE, F,                                                   \
                                    NCMT, TIME, AMT, RATE, II, EVID, CMT, ADDL, SS, THETA, BIOVAR, TLAG, \
                                    EPS_VAL, EPS_GRAD)                                                   \
  {                                                                                                      \
    {                                                                                                    \
      auto theta_v = torsten::to_var(THETA);                                                             \
      auto x1 = FUN_CPT(TIME, AMT, RATE, II, EVID, CMT, ADDL, SS, theta_v, BIOVAR, TLAG);                \
      auto x2 = FUN_ODE(F, NCMT, TIME, AMT, RATE, II, EVID, CMT, ADDL, SS, theta_v, BIOVAR, TLAG);       \
      for (size_t i = 0; i < theta_v.size(); ++i) {                                                      \
        torsten::test::test_grad(theta_v[i], x1, x2, EPS_VAL, EPS_GRAD);                                 \
      }                                                                                                  \
    }                                                                                                    \
    {                                                                                                    \
      auto biovar_v = torsten::to_var(BIOVAR);                                                           \
      auto x1 = FUN_CPT(TIME, AMT, RATE, II, EVID, CMT, ADDL, SS, THETA, biovar_v, TLAG);                \
      auto x2 = FUN_ODE(F, NCMT, TIME, AMT, RATE, II, EVID, CMT, ADDL, SS, THETA, biovar_v, TLAG);       \
      for (size_t i = 0; i < biovar_v.size(); ++i) {                                                     \
        torsten::test::test_grad(biovar_v[i], x1, x2, EPS_VAL, EPS_GRAD);                                \
      }                                                                                                  \
    }                                                                                                    \
    {                                                                                                    \
      auto tlag_v = torsten::to_var(TLAG);                                                               \
      auto x1 = FUN_CPT(TIME, AMT, RATE, II, EVID, CMT, ADDL, SS, THETA, BIOVAR, tlag_v);                \
      auto x2 = FUN_ODE(F, NCMT, TIME, AMT, RATE, II, EVID, CMT, ADDL, SS, THETA, BIOVAR, tlag_v);       \
      for (size_t i = 0; i < tlag_v.size(); ++i) {                                                       \
        torsten::test::test_grad(tlag_v[i], x1, x2, EPS_VAL, EPS_GRAD);                                  \
      }                                                                                                  \
    }                                                                                                    \
  }

#define TORSTEN_ODE_GRAD_THETA_TEST(FUN, F, NCMT, TIME, AMT, RATE, II, EVID, CMT, ADDL, SS, THETA, BIOVAR, TLAG,  \
                                   RTOL, ATOL, NSTEP,                                                             \
                                   H, EPS_VAL, EPS_RTOL, EPS_ATOL)                                                \
  {                                                                                                               \
    auto f1 = [&] (std::vector<std::vector<double> >& x) {                                                        \
      return FUN(F, NCMT, TIME, AMT, RATE, II, EVID, CMT, ADDL, SS, x, BIOVAR, TLAG, 0, RTOL, ATOL, NSTEP);       \
    };                                                                                                            \
    auto f2 = [&] (std::vector<std::vector<stan::math::var> >& x) {                                               \
      return FUN(F, NCMT, TIME, AMT, RATE, II, EVID, CMT, ADDL, SS, x, BIOVAR, TLAG, 0, RTOL, ATOL, NSTEP);       \
    };                                                                                                            \
    torsten::test::test_grad(f1, f2, THETA, H, EPS_VAL, EPS_RTOL, EPS_ATOL);                                      \
  }

#define TORSTEN_ODE_GRAD_BIOVAR_TEST(FUN, F, NCMT, TIME, AMT, RATE, II, EVID, CMT, ADDL, SS, THETA, BIOVAR, TLAG, \
                                   RTOL, ATOL, NSTEP,                                                             \
                                   H, EPS_VAL, EPS_RTOL, EPS_ATOL)                                                \
  {                                                                                                               \
    auto f1 = [&] (std::vector<std::vector<double> >& x) {                                                        \
      return FUN(F, NCMT, TIME, AMT, RATE, II, EVID, CMT, ADDL, SS, THETA, x, TLAG, 0, RTOL, ATOL, NSTEP);        \
    };                                                                                                            \
    auto f2 = [&] (std::vector<std::vector<stan::math::var> >& x) {                                               \
      return FUN(F, NCMT, TIME, AMT, RATE, II, EVID, CMT, ADDL, SS, THETA, x, TLAG, 0, RTOL, ATOL, NSTEP);        \
    };                                                                                                            \
    torsten::test::test_grad(f1, f2, BIOVAR, H, EPS_VAL, EPS_RTOL, EPS_ATOL);                                     \
  }


#define TORSTEN_ODE_GRAD_TLAG_TEST(FUN, F, NCMT, TIME, AMT, RATE, II, EVID, CMT, ADDL, SS, THETA, BIOVAR, TLAG,   \
                                   RTOL, ATOL, NSTEP,                                                             \
                                   H, EPS_VAL, EPS_RTOL, EPS_ATOL)                                                \
  {                                                                                                               \
    auto f1 = [&] (std::vector<std::vector<double> >& x) {                                                        \
      return FUN(F, NCMT, TIME, AMT, RATE, II, EVID, CMT, ADDL, SS, THETA, BIOVAR, x, 0, RTOL, ATOL, NSTEP);      \
    };                                                                                                            \
    auto f2 = [&] (std::vector<std::vector<stan::math::var> >& x) {                                               \
      return FUN(F, NCMT, TIME, AMT, RATE, II, EVID, CMT, ADDL, SS, THETA, BIOVAR, x, 0, RTOL, ATOL, NSTEP);      \
    };                                                                                                            \
    torsten::test::test_grad(f1, f2, TLAG, H, EPS_VAL, EPS_RTOL, EPS_ATOL);                                       \
  }

#define TORSTEN_ODE_GRAD_RATE_TEST(FUN, F, NCMT, TIME, AMT, RATE, II, EVID, CMT, ADDL, SS, THETA, BIOVAR, TLAG,   \
                                   RTOL, ATOL, NSTEP,                                                             \
                                    H, EPS_VAL, EPS_RTOL, EPS_ATOL)                                               \
  {                                                                                                               \
    auto f1 = [&] (std::vector<double>& x) {                                                                      \
      return FUN(F, NCMT, TIME, AMT, x, II, EVID, CMT, ADDL, SS, THETA, BIOVAR, TLAG, 0, RTOL, ATOL, NSTEP);      \
    };                                                                                                            \
    auto f2 = [&] (std::vector<stan::math::var>& x) {                                                             \
      return FUN(F, NCMT, TIME, AMT, x, II, EVID, CMT, ADDL, SS, THETA, BIOVAR, TLAG, 0, RTOL, ATOL, NSTEP);      \
    };                                                                                                            \
    torsten::test::test_grad(f1, f2, RATE, H, EPS_VAL, EPS_RTOL, EPS_ATOL);                                       \
  }

#endif
