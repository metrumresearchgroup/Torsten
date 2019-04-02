#ifndef STAN_MATH_TORSTEN_MPI_PRECOMPUTED_GRADIENTS_HPP
#define STAN_MATH_TORSTEN_MPI_PRECOMPUTED_GRADIENTS_HPP

#include <stan/math/torsten/model_solve_d.hpp>
#include <stan/math/prim/mat/fun/typedefs.hpp>

namespace torsten {
  namespace mpi {

    /*
     * Generate a matrix with @c var entries that have given
     * value and gradients. The value and gradients are provided
     * through @c MatrixXd, with each row consisting of multiple
     * records in the format (value, grad1, grad2, grad3...).
     * @param d input matrix data with each row consisting of
     * multiple records, each record of format (value, grad1, grad2...)
     * @param n number of record in each row
     * @param nrec length of each record, so that @c nrec*n=d.cols()
     * @return @c var matrix with given value and gradients
     */
    // inline stan::math::matrix_v
    // precomputed_gradients(const Eigen::MatrixXd& d, const std::vector<stan::math::var>& operands) {
    //   const int nrec = operands.size() + 1;
    //   if (d.cols() % nrec != 0) {
    //     std::stringstream msg;
    //     static const char* expr("n * operands.size()");
    //     msg << "; column number expression = " << expr;
    //     static const char* caller("torten::mpi::precomputed_gradients");
    //     std::string msg_str(msg.str());
    //     stan::math::invalid_argument(caller, "d", d.cols(), "must have n * operands.size() columns, but is ",
    //                                  msg_str.c_str());
    //   }

    //   const int n = d.cols() / nrec;
    //   Eigen::Matrix<stan::math::var, -1, -1> res(d.rows(), n);
    //   std::vector<double> g(nrec - 1);
    //   for (int i = 0 ; i < d.rows(); ++i) {
    //     for (int j = 0; j < n; ++j) {
    //       for (int l = 0 ; l < nrec - 1; ++l) g[l] = d(i, j * nrec + l + 1);
    //       res(i, j) = precomputed_gradients(d(i, j * nrec), operands, g);
    //     }
    //   }
    //   return res;
    // }

    /*
     * Generate a Eigen::Vector with @c var entries that have given
     * value and gradients. The value and gradients are provided
     * through @c VectorXd consisting of multiple
     * records in the format (value, grad1, grad2, grad3...).
     * @param d input vector data with consisting of
     * multiple records, each record of format (value, grad1, grad2...)
     * @return @c var matrix with given value and gradients
     */
    inline stan::math::vector_v
    precomputed_gradients(const Eigen::VectorXd& d, const std::vector<stan::math::var>& operands) {
      const int nrec = operands.size() + 1;
      if (d.size() % nrec != 0) {
        std::stringstream msg;
        static const char* expr("n * operands.size()");
        msg << "; column number expression = " << expr;
        static const char* caller("torten::mpi::precomputed_gradients");
        std::string msg_str(msg.str());
        stan::math::invalid_argument(caller, "d", d.size(), "must have size = n * operands.size(), but is ",
                                     msg_str.c_str());
      }

      const int n = d.size() / nrec;
      stan::math::vector_v res(n);
      std::vector<double> g(nrec - 1);
      for (int j = 0; j < n; ++j) {
        for (int l = 0 ; l < nrec - 1; ++l) g[l] = d(j * nrec + l + 1);
        res(j) = precomputed_gradients(d(j * nrec), operands, g);
      }
      return res;
    }

    /*
     * Generate a matrix with @c var entries that have given
     * value and gradients. The value and gradients are provided
     * through @c MatrixXd, with each row consisting of multiple
     * records in the format (value, grad1, grad2, grad3...).
     * @param d input matrix data with each row consisting of
     * multiple records, each record of format (value, grad1, grad2...)
     * @param n number of record in each row
     * @param nrec length of each record, so that @c nrec*n=d.cols()
     * @return @c var matrix with given value and gradients
     */
    // inline stan::math::matrix_v
    // precomputed_gradients(const Eigen::MatrixXd& d, const stan::math::matrix_v& operands) {
    //   const int nrec = operands.cols() + 1;
    //   if (d.cols() % nrec != 0) {
    //     std::stringstream msg;
    //     static const char* expr("n * operands.cols()");
    //     msg << "; column number expression = " << expr;
    //     static const char* caller("torten::mpi::precomputed_gradients");
    //     std::string msg_str(msg.str());
    //     stan::math::invalid_argument(caller, "d", d.cols(), "must have n * operands.cols() columns, but is ",
    //                                  msg_str.c_str());
    //   }
    //   if (d.rows() != operands.rows()) {
    //     std::stringstream msg;
    //     static const char* expr("operands.rows()");
    //     msg << "; row number expression = " << expr;
    //     static const char* caller("torten::mpi::precomputed_gradients");
    //     std::string msg_str(msg.str());
    //     stan::math::invalid_argument(caller, "d", d.rows(), "must have operands.rows() rows, but is ",
    //                                  msg_str.c_str());
    //   }

    //   const int n = d.cols() / nrec;
    //   stan::math::matrix_v res(d.rows(), n);
    //   std::vector<double> g(nrec - 1);
    //   for (int i = 0 ; i < d.rows(); ++i) {
    //     std::vector<stan::math::var> v_i(operands.cols());
    //     for (int l = 0 ; l < nrec - 1; ++l) v_i[l] = operands(i, l);
    //     for (int j = 0; j < n; ++j) {
    //       for (int l = 0 ; l < nrec - 1; ++l) g[l] = d(i, j * nrec + l + 1);
    //       res(i, j) = precomputed_gradients(d(i, j * nrec), v_i, g);
    //     }
    //   }
    //   return res;
    // }
  }
}
#endif

    
