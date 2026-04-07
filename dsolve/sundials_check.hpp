#ifndef TORSTEN_DSOLVE_SUNDIALS_CHECK_HPP
#define TORSTEN_DSOLVE_SUNDIALS_CHECK_HPP

#include <sstream>
#include <kinsol/kinsol.h>
#include <cvodes/cvodes.h>

#define CHECK_SUNDIALS_CALL(call) sundials_check(call, #call)
#define CHECK_KINSOL_CALL(call) kinsol_check(call, #call)

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
    throw std::domain_error(ss.str());
  }
}

/**
 * check sundials::kinsol return flag & throw runtime error
 *
 * @param[in] flag routine return flag
 * @param[in] func routine name
 */
inline void kinsol_check(int flag, const char* func) {
  if (flag < 0) {
    std::ostringstream ss;
    // std::string flag_name = KINGetLinReturnFlagName(-11);
    ss << func << " failed with error flag " << flag << ": " << KINGetReturnFlagName(flag);
    throw std::domain_error(ss.str());
  }
}


#endif
