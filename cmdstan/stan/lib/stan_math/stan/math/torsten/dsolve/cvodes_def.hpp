#ifndef STAN_MATH_TORSTEN_DSOLVE_CVODES_DEF_HPP
#define STAN_MATH_TORSTEN_DSOLVE_CVODES_DEF_HPP

#include <cvodes/cvodes.h>
#include <cvodes/cvodes_direct.h>
#include <nvector/nvector_serial.h>
#include <stan/math/torsten/cvodes_sens_method.hpp>

#ifndef TORSTEN_CV_ISM
#define TORSTEN_CV_ISM CV_STAGGERED
#endif

#ifndef TORSTEN_CV_SENS
#define TORSTEN_CV_SENS AD
#endif

namespace torsten {
  namespace dsolve {

    /**
     * define CVODES constants
     * 
     */
    template<PMXCvodesSensMethod Sm, int Lmm, int ism>
    struct cvodes_def {
      static constexpr PMXCvodesSensMethod cv_sm = Sm; 
      static constexpr int cv_lmm = Lmm; 
      static constexpr int cv_ism = ism; 
    };
  }
}

#endif
