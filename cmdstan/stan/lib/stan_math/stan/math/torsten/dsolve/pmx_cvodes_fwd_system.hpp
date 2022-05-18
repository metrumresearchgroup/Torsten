#ifndef STAN_MATH_TORSTEN_DSOLVE_CVODES_FWD_SYSTEM_HPP
#define STAN_MATH_TORSTEN_DSOLVE_CVODES_FWD_SYSTEM_HPP

#include <stan/math/prim/err/check_size_match.hpp>
#include <stan/math/torsten/dsolve/pmx_cvodes_system.hpp>
#include <stan/math/torsten/pmx_csda.hpp>
#include <stan/math/torsten/cvodes_sens_method.hpp>
#include <stan/math/torsten/dsolve/cvodes_def.hpp>

namespace torsten {
  namespace dsolve {
    /**
     * CVODES ODE system with forward sensitivity calculation
     *
     * @tparam F type of functor for ODE residual.
     * @tparam Tts type of time
     * @tparam Ty0 type of initial unknown values.
     * @tparam Tpar type of parameters.
     * @tparam Lmm method of integration(CV_ADAMS or CV_BDF)
     * @tparam Sm method of sensitivity calculatioin, choose among @c PMXCvodesSensMethod.
     */
    template <typename F, typename Tts, typename Ty0, typename Tpar, typename cv_def>  // NOLINT
    class PMXCvodesFwdSystem;

    /**
     * when use CSDA to calculate sensitivity, the
     * user-supplied RHS function of the ODE system must be
     * able to operate on complex numbers, most of current Stan
     * math functions do not support this.
     */
    template <typename F, typename Tts, typename Ty0, typename Tpar, int Lmm, int ism>
    class PMXCvodesFwdSystem<F, Tts, Ty0, Tpar, cvodes_def<CSDA, Lmm, ism>> :
      public PMXCvodesSystem<F, Tts, Ty0, Tpar, cvodes_def<CSDA, Lmm, ism>> {
    public:
      using Ode = PMXCvodesFwdSystem<F, Tts, Ty0, Tpar, cvodes_def<CSDA, Lmm, ism>>;
      static constexpr bool need_fwd_sens = Ode::is_var_y0 || Ode::is_var_par;
      PMXOdeService<Ode>& serv;
    private:
      N_Vector* nv_ys_;
      std::vector<std::complex<double> >& yy_cplx_;
      std::vector<std::complex<double> >& theta_cplx_;
      std::vector<std::complex<double> >& fval_cplx_;
    public:
      /**
       * Construct CVODES ODE system from initial condition and parameters
       *
       * @param[in] f ODE residual functor
       * @param[in] y0 initial condition
       * @param[in] theta parameters of the base ODE
       * @param[in] x_r continuous data vector for the ODE
       * @param[in] x_i integer data vector for the ODE
       * @param[in] msgs stream to which messages are printed
       */
      PMXCvodesFwdSystem(PMXOdeService<Ode>& serv0,
                         const F& f,
                         double t0,
                         const std::vector<Tts>& ts,
                         const std::vector<Ty0>& y0,
                         const std::vector<Tpar>& theta,
                         const std::vector<double>& x_r,
                         const std::vector<int>& x_i,
                         std::ostream* msgs) :
        PMXCvodesSystem<F, Tts, Ty0, Tpar, cvodes_def<CSDA, Lmm, ism>>(serv0, f, t0, ts, y0, theta, x_r, x_i, msgs),
        nv_ys_(serv0.nv_ys),
        serv(serv0),
        yy_cplx_(serv0.yy_cplx),
        theta_cplx_(serv0.theta_cplx),
        fval_cplx_(serv0.fval_cplx)
      {}

      /**
       * Dummy destructor. Deallocation of CVODES memory is done
       * in @c PMXOdeService.
       */
      ~PMXCvodesFwdSystem() {
      }

      // /**
      //  * return N_Vector pointer array of sensitivity
      //  */
      // N_Vector* nv_ys() { return nv_ys_; }

      /**
       * convert to void pointer for CVODES callbacks
       */
      void* to_user_data() {  // prepare to inject ODE info
        return static_cast<void*>(this);
      }

      // void reset_sens_mem() {
      //   serv.reset_sens_mem();
      // }

      /**
       * Calculate sensitivity rhs using CVODES vectors. The
       * internal workspace is allocated by @c PMXOdeService.
       * We use CSDA to compute senstivity, so we need to
       * generate complex version of parameters.
       */
      void eval_sens_rhs(int ns, double t, N_Vector y, N_Vector ydot,
                         N_Vector* ys, N_Vector* ysdot,
                         N_Vector temp1, N_Vector temp2) {
        using std::complex;
        using cplx = complex<double>;
        using B = PMXCvodesSystem<F, Tts, Ty0, Tpar, cvodes_def<CSDA, Lmm, ism>>;
        const int n = B::N_;
        const double h = 1.E-20;
        for (int i = 0; i < ns; ++i) {
          for (int j = 0; j < n; ++j) {
            yy_cplx_[j] = cplx(NV_Ith_S(y, j), h * NV_Ith_S(ys[i], j));
          }

          /* if y0 is the only parameter, use tangent linear
           * model(TLM). Otherwise use full forward sensitivity model.
           * Note that when both y0 and theta are
           * parameters, the first n vector of ys are for y0 sensitivity.
           */
          if (B::is_var_y0 && i < n) {
            fval_cplx_ =
              B::f_(t, yy_cplx_, B::theta_dbl_, B::x_r_, B::x_i_, B::msgs_);
          } else {
            std::transform(B::theta_dbl_.begin(),
                           B::theta_dbl_.end(),
                           theta_cplx_.begin(),
                           [](double r) -> cplx {return cplx(r, 0.0); });
            theta_cplx_.at(i - B::is_var_y0 * n) += cplx(0.0, h);
            fval_cplx_ =
              B::f_(t, yy_cplx_, theta_cplx_, B::x_r_, B::x_i_, B::msgs_);
          }

          std::transform(fval_cplx_.begin(),
                         fval_cplx_.end(),
                         B::fval_.begin(),
                         [&h](cplx x) -> double { return std::imag(x)/h; });
          for (int j = 0; j < n; ++j) NV_Ith_S(ysdot[i], j) = B::fval_[j];
        }
      }
    };

    /**
     * use autodiff to calculate sensitivity
     */
    template <typename F, typename Tts, typename Ty0, typename Tpar, int Lmm, int ism>
    class PMXCvodesFwdSystem<F, Tts, Ty0, Tpar, cvodes_def<AD, Lmm, ism>> :
      public PMXCvodesSystem<F, Tts, Ty0, Tpar, cvodes_def<AD, Lmm, ism>> {
    public:
      using Ode = PMXCvodesFwdSystem<F, Tts, Ty0, Tpar, cvodes_def<AD, Lmm, ism>>;
      static constexpr bool need_fwd_sens = Ode::is_var_y0 || Ode::is_var_par;
      PMXOdeService<Ode>& serv;
    public:
      /**
       * Construct CVODES ODE system from initial condition and parameters
       *
       * @param[in] f ODE residual functor
       * @param[in] y0 initial condition
       * @param[in] theta parameters of the base ODE
       * @param[in] x_r continuous data vector for the ODE
       * @param[in] x_i integer data vector for the ODE
       * @param[in] msgs stream to which messages are printed
       */
      PMXCvodesFwdSystem(PMXOdeService<Ode>& serv0,
                           const F& f,
                           double t0,
                           const std::vector<Tts>& ts,
                           const std::vector<Ty0>& y0,
                           const std::vector<Tpar>& theta,
                           const std::vector<double>& x_r,
                           const std::vector<int>& x_i,
                           std::ostream* msgs) :
        PMXCvodesSystem<F, Tts, Ty0, Tpar, cvodes_def<AD, Lmm, ism>>(serv0, f, t0, ts, y0, theta, x_r, x_i, msgs),
        serv(serv0)
      {}

      /**
       * Dummy destructor. Deallocation of CVODES memory is done
       * in @c PMXOdeService.
       */
      ~PMXCvodesFwdSystem() {}
    };


    template<typename F_type, typename t_type, typename initial_type, typename param_type>
    using PMXOdeSystem = PMXCvodesFwdSystem<F_type, t_type, initial_type, param_type,
                                                        cvodes_def<TORSTEN_CV_SENS, CV_ADAMS, TORSTEN_CV_ISM>>;

    template<typename F_type, typename t_type, typename initial_type, typename param_type>
    using PMXOdeSystem = PMXCvodesFwdSystem<F_type, t_type, initial_type, param_type,
                                                        cvodes_def<TORSTEN_CV_SENS, CV_BDF, TORSTEN_CV_ISM>>;

  }  // namespace dsolve
}  // namespace torsten

#endif
