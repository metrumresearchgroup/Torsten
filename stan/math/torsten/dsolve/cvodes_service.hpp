#ifndef STAN_MATH_TORSTEN_DSOLVE_CVODES_SERVICE_HPP
#define STAN_MATH_TORSTEN_DSOLVE_CVODES_SERVICE_HPP

#include <stan/math/rev/scal/meta/is_var.hpp>
#include <stan/math/torsten/dsolve/sundials_check.hpp>
#include <cvodes/cvodes.h>
#include <cvodes/cvodes_direct.h>
#include <sunmatrix/sunmatrix_dense.h>
#include <sunlinsol/sunlinsol_dense.h>
#include <nvector/nvector_serial.h>
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
    template <typename Ode>
    struct cvodes_service {
      const size_t N;
      const size_t M;
      const size_t ns;
      N_Vector nv_y;
      N_Vector* nv_ys;
      std::vector<double> y;
      std::vector<double> fval;
      void* mem;
      SUNMatrix A;
      SUNLinearSolver LS;
      std::vector<std::complex<double> > yy_cplx;
      std::vector<std::complex<double> > theta_cplx;
      std::vector<std::complex<double> > fval_cplx;

      /**
       * Construct CVODES ODE mem & workspace
       *
       * @param[in] n ODE system size
       * @param[in] m length of parameter theta
       * @param[in] f ODE RHS function
       */
      cvodes_service(int n, int m, CVRhsFn f) :
        N(n),
        M(m),
        ns((Ode::is_var_y0 ? n : 0) + (Ode::is_var_par ? m : 0)),
        nv_y(N_VNew_Serial(n)),
        nv_ys(nullptr),
        y(n),
        fval(n),
        mem(CVodeCreate(Ode::lmm_type, CV_NEWTON)),
        yy_cplx(n),
        theta_cplx(m),
        fval_cplx(n)
      {
        const double t0 = 0.0;
        for (int i = 0; i < n; ++i)  N_VConst(RCONST(0.0), nv_y);

        /*
         * allocate sensitivity array if need fwd sens calculation
         */ 
        if (Ode::need_fwd_sens) {
          nv_ys = N_VCloneVectorArray(ns, nv_y);
          for (size_t i = 0; i < ns; ++i) N_VConst(RCONST(0.0), nv_ys[i]);
        }

        /*
         * initialize cvodes system and allocate linear solver mem
         */ 
        CHECK_SUNDIALS_CALL(CVodeInit(mem, Ode::rhs(), t0, nv_y));
        A = SUNDenseMatrix(n, n);
        LS = SUNDenseLinearSolver(nv_y, A);
        CHECK_SUNDIALS_CALL(CVDlsSetLinearSolver(mem, LS, A));
      }

      ~cvodes_service() {
        SUNLinSolFree(LS);
        SUNMatDestroy(A);
        N_VDestroy(nv_y);
        CVodeFree(&mem);
      }
    };

  }
}

#endif
