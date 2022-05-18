#ifndef STAN_MATH_TORSTEN_DSOLVE_CVODES_SERVICE_HPP
#define STAN_MATH_TORSTEN_DSOLVE_CVODES_SERVICE_HPP

#include <stan/math/rev/core/recover_memory.hpp>
#include <stan/math/rev/meta/is_var.hpp>
#include <stan/math/prim/meta/require_helpers.hpp>
#include <stan/math/torsten/dsolve/sundials_check.hpp>
#include <stan/math/torsten/dsolve/ode_func_type.hpp>
#include <stan/math/torsten/dsolve/pmx_ode_system.hpp>
#include <cvodes/cvodes.h>
// #include <arkode/arkode.h>
// #include <arkode/arkode_erkstep.h>
#include <sundials/sundials_context.h>
#include <nvector/nvector_serial.h>
#include <sunmatrix/sunmatrix_dense.h>
#include <sunlinsol/sunlinsol_dense.h>
#include <type_traits>
#include <ostream>
#include <vector>
#include <algorithm>

namespace torsten {
  namespace dsolve {

    /** For each type of Ode(with different rhs functor F and
     * senstivity parameters), we allocate mem and workspace for
     * cvodes. This service manages the
     * allocation/deallocation, so ODE systems only request
     * service by injection.
     * @tparam ode ode type
     * @tparam lmm_type CVODES solver type (BDF & ADAMS)
     * @tparam butcher_tab AKRODE Butcher table
     */
    template <typename Ode, int lmm_type, int butcher_tab = 0, typename = void>
    struct PMXOdeService {
      sundials::Context sundials_context_;
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
        sundials_context_(),
        ns(ns0),
        nv_y(N_VNew_Serial(n, sundials_context_)),
        nv_ys(nullptr),
        mem(CVodeCreate(lmm_type, sundials_context_)),
        A(SUNDenseMatrix(n, n, sundials_context_)),
        LS(SUNLinSol_Dense(nv_y, A, sundials_context_)),
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

    /** 
     * ARKODE's Butcher table: 0-99 are for ERK.
     * 
     */
    template<int butcher_tab>
    struct is_erk_tab {
      static constexpr bool value = butcher_tab >= 0 && butcher_tab < 100;
    };

    // template <typename Ode, int butcher_tab>
    // struct PMXOdeService<Ode, 0, butcher_tab,
    //                      stan::require_t<is_erk_tab<butcher_tab> > > {
    //   sundials::Context sundials_context_;
    //   int ns;
    //   N_Vector nv_y;
    //   void* mem;

    //   /**
    //    * Construct CVODES ODE mem & workspace
    //    *
    //    * @param[in] n ODE system size
    //    * @param[in] m length of parameter theta
    //    * @param[in] f ODE RHS function
    //    */
    //   PMXOdeService(int n, int m, int ns0, Ode& ode) :
    //     sundials_context_(),
    //     ns(ns0),
    //     nv_y(N_VNew_Serial(ode.system_size, sundials_context_)),
    //     mem(ERKStepCreate(Ode::cvodes_rhs, 0.0, nv_y))
    //   {
    //     N_VConst(RCONST(0.0), nv_y);
    //     CHECK_SUNDIALS_CALL(ERKStepSetUserData(mem, static_cast<void*>(&ode)));
    //   }

    //   ~PMXOdeService() {
    //     ERKStepFree(&mem);    // Free integrator memory
    //     N_VDestroy(nv_y);        // Free y vector
    //   }
    // };
  }
}

#endif
