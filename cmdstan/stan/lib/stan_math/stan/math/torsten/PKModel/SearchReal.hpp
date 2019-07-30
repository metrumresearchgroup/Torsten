#ifndef STAN_MATH_TORSTEN_PKMODEL_SEARCHREAL_HPP
#define STAN_MATH_TORSTEN_PKMODEL_SEARCHREAL_HPP

#include <cassert>
#include <iostream>
#include <vector>
#include <algorithm>

namespace torsten {

int min(int a, int b);  // forward declare

/**
 * Returns the position of the largest value smaller or greater than srchNum
 * Assumes that v is sorted.
 * The numltm is used to set an upper limit for the index
 * the function can return. If numltm is greater than the size of v,
 * there is no upper limit and the function searches the entire vector.
 *
 * @tparam: type of scalar in input vector
 * @param[in]: v searched vector
 * @param[in]: numltm maximum index
 * @param[in]: srchNum searched Number
 * @return: index of largest value <= srchNum
 *
 */
template<typename T0, typename T1>
inline int SearchReal(const std::vector<T0>& v, int numltm, T1 srchNum) {
  int first = 0, last, mid, real_limit;

  assert(numltm >= 0);
  // limit cannot exceed size of searched vector
  real_limit = min(numltm, v.size());
  last = real_limit - 1;

  if (srchNum < v[first]) mid = first;
  else
    if (srchNum >= v[last]) {
      mid = last + 1;
  } else {
    while (first <= last) {
      mid = (first + last) / 2;
      if (srchNum < v[mid]) last = mid - 1;
      else if (srchNum > v[mid]) first = mid + 1;
      else
        first = last + 1;
    }

    while (srchNum >= v[mid]) mid += 1;
  }
  return mid;
}

inline int min(int a, int b) {
  if (a < b) return a;
  else
    return b;
}

}

#endif
