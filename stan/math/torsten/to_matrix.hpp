#ifndef STAN_MATH_TORSTEN_TO_MATRIX_HPP
#define STAN_MATH_TORSTEN_TO_MATRIX_HPP

#include <stan/math/prim/mat/fun/to_matrix.hpp>
#include <nvector/nvector_serial.h>

namespace torsten {

  /*
   * convert a vector of vectors to col-majored matrix, with
   * each vector to each col of the matrix. This is
   * equivalent to have tranposing the output of stan's @c to_matrix.
   */
  template <typename T>
  inline Eigen::Matrix<T, -1, -1>
  to_matrix(const std::vector<std::vector<T> >& vv) {
    Eigen::Matrix<T, -1, -1> res(vv[0].size(), vv.size());
    for (size_t i = 0; i < vv.size(); ++i) {
      for (size_t j = 0; j < vv[i].size(); ++j) {
        res(j, i) = vv[i][j];
      }
    }
    return res;
  }
  
  /**
   * copy NV_Vector* array to Eigen::MatrixXd
   *
   * @param[in] nv N_Vector* array.
   * @param[in] nv_size length of nv.
   * @return Eigen::MatrixXd.
   */
  inline Eigen::MatrixXd to_matrix(const N_Vector* nv,
                                   const size_t& nv_size) {
    size_t m = nv_size;
    size_t n = NV_LENGTH_S(nv[0]);
    stan::math::matrix_d res(n, m);
    for (size_t j = 0; j < m; ++j) {
      auto nvp = N_VGetArrayPointer(nv[j]);
      for (size_t i = 0; i < n; ++i) {
        res(i, j) = nvp[i];
      }
    }
    return res;
  }

}


#endif
