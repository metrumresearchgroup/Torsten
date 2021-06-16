#ifndef STAN_MATH_TORSTEN_DSOLVE_CVODES_SERVICE_HPP
#define STAN_MATH_TORSTEN_DSOLVE_CVODES_SERVICE_HPP

#include <stan/math/rev/core/recover_memory.hpp>
#include <stan/math/rev/meta/is_var.hpp>
#include <stan/math/torsten/dsolve/sundials_check.hpp>
#include <stan/math/torsten/dsolve/ode_func_type.hpp>
#include <stan/math/torsten/dsolve/pmx_ode_system.hpp>
#include <cvodes/cvodes.h>
#include <nvector/nvector_serial.h>
#include <sunmatrix/sunmatrix_dense.h>
#include <sunlinsol/sunlinsol_dense.h>
#include <ostream>
#include <vector>
#include <algorithm>

namespace torsten {
  namespace dsolve {

    /* For each type of Ode(with different rhs functor F and
     * senstivity parameters), we allocate mem and workspace for
     * cvodes. This service manages the
     * allocation/deallocation, so ODE systems only request
     * service by injection.
     */
    template <typename Ode, int lmm_type>
    struct PMXOdeService {
      int ns;
      N_Vector nv_y;
      N_Vector* nv_ys;
      void* mem;
      SUNMatrix A;
      SUNLinearSolver LS;
      std::vector<std::complex<double> > yy_cplx;
      std::vector<std::complex<double> > theta_cplx;
      std::vector<std::complex<double> > fval_cplx;
      bool sens_inited;

      /**
       * Construct CVODES ODE mem & workspace
       *
       * @param[in] n ODE system size
       * @param[in] m length of parameter theta
       * @param[in] f ODE RHS function
       */
      PMXOdeService(int n, int m, int ns0, Ode& ode) :
        ns(ns0),
        nv_y(N_VNew_Serial(n)),
        nv_ys(nullptr),
        mem(CVodeCreate(lmm_type)),
        A(SUNDenseMatrix(n, n)),
        LS(SUNLinSol_Dense(nv_y, A)),
        yy_cplx(n),
        theta_cplx(m),
        fval_cplx(n),
        sens_inited(false)
      {
        const double t0 = 0.0;
        N_VConst(RCONST(0.0), nv_y);

        /*
         * allocate sensitivity array if need fwd sens calculation
         */ 
        if (Ode::use_fwd_sens) {
          nv_ys = N_VCloneVectorArray(ns, nv_y);
        }

        /*
         * initialize cvodes system and attach linear solver
         */ 
        CHECK_SUNDIALS_CALL(CVodeInit(mem, Ode::cvodes_rhs, t0, nv_y));
        CHECK_SUNDIALS_CALL(CVodeSetUserData(mem, static_cast<void*>(&ode)));
        CHECK_SUNDIALS_CALL(CVodeSetLinearSolver(mem, LS, A));

        // default CV_STAGGERED will be overwritten later during reinit
        if (Ode::use_fwd_sens) {
          CHECK_SUNDIALS_CALL(CVodeSensInit(mem, ns, CV_STAGGERED, Ode::cvodes_sens_rhs, nv_ys)); 
        }
      }

      ~PMXOdeService() {
        SUNLinSolFree(LS);
        SUNMatDestroy(A);
        CVodeFree(&mem);
        if (Ode::use_fwd_sens) {
          CVodeSensFree(mem);
        }
        if (Ode::use_fwd_sens) {
          N_VDestroyVectorArray(nv_ys, ns);
        }
        N_VDestroy(nv_y);
      }
    };
  }
}

#endif
