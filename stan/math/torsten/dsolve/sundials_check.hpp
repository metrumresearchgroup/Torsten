#ifndef TORSTEN_DSOLVE_SUNDIALS_CHECK_HPP
#define TORSTEN_DSOLVE_SUNDIALS_CHECK_HPP

#include <sstream>

#define CHECK_SUNDIALS_CALL(call) sundials_check(call, #call)

/**
 * check sundials return flag & throw runtime error
 *
 * @param[in] flag routine return flag
 * @param[in] func routine name
 */
inline void sundials_check(int flag, const char* func) {
  if (flag < 0) {
    std::ostringstream ss;
    ss << func << " failed with error flag " << flag;
    throw std::runtime_error(ss.str());
  }
}


#endif
