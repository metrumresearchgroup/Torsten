#ifndef STAN_MATH_TORSTEN_MPI_PRECOMPUTED_GRADIENTS_HPP
#define STAN_MATH_TORSTEN_MPI_PRECOMPUTED_GRADIENTS_HPP

#include <stan/math/rev/mat/fun/typedefs.hpp>

namespace torsten {

  /*
   * @tparam Ty @c vector type, can be either std::vector or eigen::vector
   */
  template<typename Ty = std::vector<double>,
           typename Ty_out = std::vector<stan::math::var> >
  inline Ty_out
  precomputed_gradients(const Ty& sol, const std::vector<stan::math::var>& theta) {
    const size_t n = sol.size()/(1 + theta.size());
    Ty_out y_res(n);
    std::vector<double> g(theta.size());
    for (size_t j = 0; j < n; j++) {
      for (size_t k = 0; k < theta.size(); k++) g[k] = sol[n + n * k + j];
      y_res[j] = precomputed_gradients(sol[j], theta, g);
    }
    return y_res;
  }

  template<typename Ty_out>
  inline Ty_out
  precomputed_gradients(const Eigen::MatrixXd& sol, int i_col, const std::vector<stan::math::var>& theta) {
    const size_t n = sol.rows()/(1 + theta.size());
    Ty_out y_res(n);
    std::vector<double> g(theta.size());
    for (size_t j = 0; j < n; j++) {
      for (size_t k = 0; k < theta.size(); k++) g[k] = sol(n + n * k + j, i_col);
      y_res[j] = precomputed_gradients(sol(j, i_col), theta, g);
    }
    return y_res;
  }

  inline stan::math::matrix_v
  precomputed_gradients (const Eigen::MatrixXd& sol, const std::vector<stan::math::var>& theta) {
    const size_t n = sol.rows()/(1 + theta.size());
    stan::math::matrix_v y_res(n, sol.cols());
    for (int i = 0; i < sol.cols(); i++) {
      y_res.col(i) = precomputed_gradients<stan::math::vector_v>(sol, i, theta);
    }
    return y_res;
  }

  namespace mpi {
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
    precomputed_gradients(const Eigen::VectorXd& sol, const std::vector<stan::math::var>& theta) {
      const size_t n = sol.size()/(1 + theta.size());
      stan::math::vector_v y_res(n);
      std::vector<double> g(theta.size());
      for (size_t j = 0; j < n; j++) {
        for (size_t k = 0; k < theta.size(); k++) g[k] = sol(n + n * k + j);
        y_res(j) = precomputed_gradients(sol(j), theta, g);
      }
      return y_res;
    }
  }
}
#endif

    
