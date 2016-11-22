#ifndef STAN_MATH_TORSTEN_PKMODEL_EXTRACTVECTOR_HPP
#define STAN_MATH_TORSTEN_PKMODEL_EXTRACTVECTOR_HPP

#include <iostream>
#include <string>
#include <Eigen/Dense>

/**
 * Extracts the nth column or row from a matrix and
 * returns it as an std::vector. Indexing starts at 0, meaning the
 * first row or column is found for n = 0.
 *
 * @param[in]: matrix from which the row or column will be extracted
 * @param[in]: n the row or column number
 * @param[in]: str a string that specifies whether the function extracts
 *             a column or a row.
 * @return: output dvector
 */
template <typename T>
vector<T> ExtractVector(Matrix<T, Dynamic, 1> matrix, int n, string str) {
  using std::vector;
  using std::string;
  using namespace Eigen;

  assert((str == "row")||(str == "col"));
  int length;
  vector<T> ExtVec;

  if(str == "col") {  //extract a column.
    length = matrix.rows();
    assert(n < matrix.cols());
    ExtVec.resize(length);
    for(int i = 0; i < length; i++) ExtVec[i] = matrix(i, n);
  }

  else {  // extract a row
    length = matrix.cols();
    assert(n < matrix.rows());
    ExtVec.resize(length);
    for(int i = 0; i < length; i++) ExtVec[i] = matrix(n, i);
  }

	return ExtVec;
}

#endif
