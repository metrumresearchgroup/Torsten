#ifndef STAN_MATH_TORSTEN_DSOLVE_CVODES_SERVICE_HPP
#define STAN_MATH_TORSTEN_DSOLVE_CVODES_SERVICE_HPP

#include <stan/math/rev/meta/is_var.hpp>
#include <stan/math/torsten/dsolve/sundials_check.hpp>
#include <stan/math/torsten/dsolve/cvodes_rhs.hpp>
#include <stan/math/torsten/dsolve/cvodes_sens_rhs.hpp>
#include <stan/math/torsten/dsolve/ode_func_type.hpp>
#include <stan/math/torsten/dsolve/ode_forms.hpp>
#include <cvodes/cvodes.h>
#include <sunmatrix/sunmatrix_dense.h>
#include <sunlinsol/sunlinsol_dense.h>
#include <ostream>
#include <vector>
#include <algorithm>

namespace torsten {
  namespace dsolve {

    template<typename Ode>
    struct cvodes_user_data {
      using F = typename ode_func<Ode>::type;

      const F* f;
      const size_t N;
      const size_t M;
      const size_t ns;
      const size_t fwd_ode_dim;
      std::vector<double> y;
      std::vector<double> fval;
      std::vector<double> theta_d;
      const std::vector<double>* px_r;
      const std::vector<int>* px_i;
      std::ostream* msgs;

      /** 
       * constructor
       * 
       * @param n dim of original ODE
       * @param m # of theta params
       * 
       */
      cvodes_user_data(int n, int m) :
        f(nullptr), N(n), M(m),
        ns((Ode::is_var_y0 ? n : 0) + (Ode::is_var_par ? m : 0)),
        fwd_ode_dim(N + N * ns),
        y(n), fval(n), theta_d(m),
        px_r(nullptr), px_i(nullptr), msgs(nullptr)
      {}

      /**
       * evaluate RHS function using current state, store
       * the result in @c N_Vector.
       */
      inline void eval_rhs(double t, N_Vector& nv_y, N_Vector& ydot) {
        for (size_t i = 0; i < N; ++i) {
         y[i] = NV_Ith_S(nv_y, i); 
        }
        fval = (*f)(t, y, theta_d, *px_r, *px_i, msgs);
        for (size_t i = 0; i < N; ++i) {
          NV_Ith_S(ydot, i) = fval[i];
        }
      }

      /**
       * Calculate sensitivity rhs using CVODES vectors. The
       * internal workspace is allocated by @c PMXOdeService.
       */
      void eval_sens_rhs(int ns, double t, N_Vector nv_y, N_Vector ydot,
                         N_Vector* ys, N_Vector* ysdot,
                         N_Vector temp1, N_Vector temp2) {
        using stan::math::var;

        for (int i = 0; i < N; ++i) y[i] = NV_Ith_S(nv_y, i);

        // initialize ysdot
        for (int i = 0; i < ns; ++i) N_VConst(0.0, ysdot[i]);

        try {
          stan::math::start_nested();

          std::vector<var> yv_work(NV_DATA_S(nv_y), NV_DATA_S(nv_y) + N);
          std::vector<var> theta_work(theta_d.begin(), theta_d.end());
          std::vector<var> fyv_work(Ode::is_var_par ?
                                    (*f)(t, yv_work, theta_work, *px_r, *px_i, msgs) :
                                    (*f)(t, yv_work, theta_d, *px_r, *px_i, msgs));

          stan::math::check_size_match("PMXOdeSystem", "dz_dt", fyv_work.size(), "states", N);

          for (int j = 0; j < N; ++j) {
            stan::math::set_zero_all_adjoints_nested();
            fyv_work[j].grad();

            // df/dy*s_i term, for i = 1...ns
            for (int i = 0; i < ns; ++i) {
              auto ysp = N_VGetArrayPointer(ys[i]);
              auto nvp = N_VGetArrayPointer(ysdot[i]);
              for (int k = 0; k < N; ++k) nvp[j] += yv_work[k].adj() * ysp[k];
            }

            // df/dp_i term, for i = n...n+m-1
            if (Ode::is_var_par) {
              for (int i = 0; i < M; ++i) {
                auto nvp = N_VGetArrayPointer(ysdot[ns - M + i]);
                nvp[j] += theta_work[i].adj();
              }
            }
          }
        } catch (const std::exception& e) {
          stan::math::recover_memory_nested();
          throw;
        }
        stan::math::recover_memory_nested();
      }

      /**
       * return a closure for CVODES residual callback using a
       * non-capture lambda.
       *
       * @tparam Ode type of Ode
       * @return RHS function for Cvodes
       */
      inline CVDlsJacFn cvodes_jac() {
        return [](realtype t, N_Vector y, N_Vector fy, SUNMatrix J, void* user_data, // NOLINT
                  N_Vector tmp1, N_Vector tmp2, N_Vector tmp3) -> int {
          cvodes_user_data<Ode>* ode = static_cast<cvodes_user_data<Ode>*>(user_data);
          ode -> eval_jac(t, y, fy, J);
          return 0;
        };
      }

      /**
       * evaluate Jacobian matrix using current state, store
       * the result in @c SUNMatrix J.
       *
       * @param t current time
       * @param y current y
       * @param fy current f(y)
       * @param J Jacobian matrix J(i,j) = df_i/dy_j
       */
      inline void eval_jac(double t, N_Vector& nv_y, N_Vector& fy, SUNMatrix& J) {
        using stan::math::var;

        try {
          stan::math::start_nested();

          std::vector<var> yv_work(NV_DATA_S(nv_y), NV_DATA_S(nv_y) + N);
          std::vector<var> fyv_work((*f)(t, yv_work, theta_d, *px_r, *px_i, msgs));

          for (int i = 0; i < N; ++i) {
            stan::math::set_zero_all_adjoints_nested();
            fyv_work[i].grad();
            for (int j = 0; j < N; ++j) {
              SM_ELEMENT_D(J, i, j) = yv_work[j].adj();
            }
          }
        } catch (const std::exception& e) {
          stan::math::recover_memory_nested();
          throw;
        }
        stan::math::recover_memory_nested();
      }
    };

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
      cvodes_user_data<Ode> user_data;
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
      PMXOdeService(int n, int m) :
        user_data(n, m),
        nv_y(N_VNew_Serial(n)),
        nv_ys(nullptr),
        mem(CVodeCreate(Ode::lmm_type)),
        A(SUNDenseMatrix(n, n)),
        LS(SUNLinSol_Dense(nv_y, A)),
        yy_cplx(n),
        theta_cplx(m),
        fval_cplx(n),
        sens_inited(false)
      {
        const double t0 = 0.0;
        for (int i = 0; i < n; ++i)  N_VConst(RCONST(0.0), nv_y);

        /*
         * allocate sensitivity array if need fwd sens calculation
         */ 
        if (Ode::need_fwd_sens) {
          nv_ys = N_VCloneVectorArray(user_data.ns, nv_y);
          for (size_t i = 0; i < user_data.ns; ++i) N_VConst(RCONST(0.0), nv_ys[i]);
        }

        /*
         * initialize cvodes system and attach linear solver
         */ 
        CHECK_SUNDIALS_CALL(CVodeInit(mem, cvodes_rhs, t0, nv_y));
        CHECK_SUNDIALS_CALL(CVodeSetUserData(mem, static_cast<void*>(&user_data)));
        CHECK_SUNDIALS_CALL(CVDlsSetLinearSolver(mem, LS, A));

        if (Ode::need_fwd_sens) {
          CHECK_SUNDIALS_CALL(CVodeSensInit(mem, user_data.ns, Ode::ism_type, cvodes_sens_rhs, nv_ys)); 
        }
      }

      ~PMXOdeService() {
        SUNLinSolFree(LS);
        SUNMatDestroy(A);
        CVodeFree(&mem);
        if (Ode::need_fwd_sens) {
          CVodeSensFree(mem);
        }
        N_VDestroyVectorArray(nv_ys, user_data.ns);
        N_VDestroy(nv_y);
      }

      static int cvodes_rhs(double t, N_Vector y, N_Vector ydot, void* user_data) {
          cvodes_user_data<Ode>* ode = static_cast<cvodes_user_data<Ode>*>(user_data);
          ode -> eval_rhs(t, y, ydot);
          return 0;
      }

      // sens
      template <bool needs_sens>
      struct cvodes_sens_rhs_impl {
        static int f(int ns, double t, N_Vector y, N_Vector ydot,
                     N_Vector* ys, N_Vector* ysdot, void* user_data,
                     N_Vector temp1, N_Vector temp2) {
            cvodes_user_data<Ode>* ode = static_cast<cvodes_user_data<Ode>*>(user_data);
            ode -> eval_sens_rhs(ns, t, y, ydot, ys, ysdot, temp1, temp2);
            return 0;
          }
      };

      template <>
      struct cvodes_sens_rhs_impl<false> {
        static int f(int ns, double t, N_Vector y, N_Vector ydot,
                     N_Vector* ys, N_Vector* ysdot, void* user_data,
                     N_Vector temp1, N_Vector temp2) {
          return 0;
        }
      };

      /**
       * return a closure for CVODES sensitivity RHS callback using a
       * non-capture lambda.
       *
       * @tparam Ode type of Ode
       * @return RHS function for Cvodes
       */
      static int cvodes_sens_rhs(int ns, double t, N_Vector y, N_Vector ydot,
                                 N_Vector* ys, N_Vector* ysdot, void* user_data,
                                 N_Vector temp1, N_Vector temp2) {
        return cvodes_sens_rhs_impl<Ode::need_fwd_sens>::f(ns, t, y, ydot, ys, ysdot, user_data, temp1, temp2);
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
