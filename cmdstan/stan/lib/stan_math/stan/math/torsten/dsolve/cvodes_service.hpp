#ifndef STAN_MATH_TORSTEN_DSOLVE_CVODES_SERVICE_HPP
#define STAN_MATH_TORSTEN_DSOLVE_CVODES_SERVICE_HPP

#include <stan/math/rev/scal/meta/is_var.hpp>
#include <stan/math/torsten/dsolve/sundials_check.hpp>
#include <stan/math/torsten/dsolve/cvodes_rhs.hpp>
#include <stan/math/torsten/dsolve/ode_forms.hpp>
#include <cvodes/cvodes.h>
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
    template <typename Ode, enum PMXOdeForms = OdeForm<Ode>::value>
    struct PMXOdeService;

    template <typename Ode>
    struct PMXOdeService<Ode, Cvodes> {
      const size_t N;
      const size_t M;
      const size_t ns;
      const size_t size;
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
      PMXOdeService(int n, int m) :
        N(n),
        M(m),
        ns((Ode::is_var_y0 ? n : 0) + (Ode::is_var_par ? m : 0)),
        size(N + N * ns),
        nv_y(N_VNew_Serial(n)),
        nv_ys(nullptr),
        y(n),
        fval(n),
        mem(CVodeCreate(Ode::lmm_type)),
        A(SUNDenseMatrix(n, n)),
        LS(SUNLinSol_Dense(nv_y, A)),
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
         * initialize cvodes system and attach linear solver
         */ 
        CHECK_SUNDIALS_CALL(CVodeInit(mem, cvodes_rhs<Ode>(), t0, nv_y));
        CHECK_SUNDIALS_CALL(CVDlsSetLinearSolver(mem, LS, A));
      }

      ~PMXOdeService() {
        SUNLinSolFree(LS);
        SUNMatDestroy(A);
        CVodeFree(&mem);
        N_VDestroyVectorArray(nv_ys, ns);
        N_VDestroy(nv_y);
      }
    };

    template <typename Ode>
    struct PMXOdeService<Ode, Odeint> {
      const size_t N;
      const size_t M;
      const size_t ns;
      const size_t size;
      std::vector<double> y;

      /**
       * Construct Boost Odeint workspace
       */
      PMXOdeService(int n, int m) :
        N(n),
        M(m),
        ns((Ode::is_var_y0 ? n : 0) + (Ode::is_var_par ? m : 0)),
        size(n + n * ns),
        y(size, 0.0)
      {}

      ~PMXOdeService() {}
    };

  }
}

#endif
