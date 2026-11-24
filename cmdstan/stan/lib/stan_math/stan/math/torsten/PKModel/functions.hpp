// Define a few simple functions. I guess I didn't deem them worthy of
// having their own file, but this might not be the best approach
// FIX ME: deprecate the print functions.
// FIX ME: each function should have its own header file
// FIX ME: using statements should be inside the scope of functions
// FIX ME: using functor for Marker instead of function with global variable

#ifndef STAN_MATH_TORSTEN_PKMODEL_FUNCTIONS_HPP
#define STAN_MATH_TORSTEN_PKMODEL_FUNCTIONS_HPP

#include <stan/math/torsten/PKModel/SearchReal.hpp>
#include <vector>

namespace torsten {

/**
 * Define some simple and handy functions.
 *
 * Marker(): prints 'marker', followed by the number of
 * the marker. Useful for tracking lines at which errors occur.
 *
 * find_time(): determines whether a value is present inside
 * a vector. Used to determine if times in the Events Schedule
 * are present in the ModelParameterHistory object.
 */

////////// FOR DEVELOPMENT PURPOSES //////////////
// int marker_count = 0;  // define global variable
// inline void Marker() {
//   std::cout << "MARKER "
//             << marker_count
//             << std::endl;
//   marker_count++;
// }

////////// REQUIRED FOR TORSTEN //////////////
template<typename T>
bool find_time(std::vector<T> v, T time) {
  bool found = false;
  int size = v.size();
  int k = SearchReal(v, size, time);  // find the index of the largest
                                      // element in v smaller than time.
  if ((k <= size) && (time == v[k-1])) found = true;
  return found;
}

}

#endif
