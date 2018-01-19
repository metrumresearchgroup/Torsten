#ifndef STAN_MATH_TORSTEN_PKMODEL_FUNCTORS_CHECK_MTI_HPP
#define STAN_MATH_TORSTEN_PKMODEL_FUNCTORS_CHECK_MTI_HPP

#include <stan/math/torsten/PKModel/Pred/unpromote.hpp>
#include <string>

namespace torsten {

/**
 * Checks arguments in a multiple truncated infusion 
 * case (arises when rate != 0 and ii > 0). Checks:
 * (1) amt > 0
 * (2) delta < ii, where delta = amt / rate
 * 
 * @tparam T0 scalar type of amt
 * @tparam T1 scalar type of delta
 * @tparam T2 scalar type of ii
 * @param amt dosing amount administered during the infusion
 * @param delta duration of the infusion
 * @param ii inter-dose interval
 * @return void
 */
template <typename T0, typename T1, typename T2>
void check_mti(const T0& amt,
               const T1& delta,
               const T2& ii,
               const char* function) {
  using stan::math::invalid_argument;

  if (!(unpromote(amt) > 0)) {
    invalid_argument(function, "Amount (amt)", amt, "is ",
                     " but must be stricly positive when ii > 0!");
  }

  if (unpromote(delta) > unpromote(ii)) {
    std::string msg = " but must be smaller than the interdose interval (ii): "  // NOLINT
    + boost::lexical_cast<std::string>(ii) + "!";
    const char* msg2 = msg.c_str();
    stan::math::invalid_argument(function,
                                 "Infusion time (F * amt / rate)", delta,
                                 "is ", msg2);
  }
}

}
#endif
