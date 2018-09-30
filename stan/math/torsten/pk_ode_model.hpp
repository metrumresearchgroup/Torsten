#ifndef STAN_MATH_TORSTEN_REFACTOR_ODE_MODEL_HPP
#define STAN_MATH_TORSTEN_REFACTOR_ODE_MODEL_HPP

#include <stan/math/torsten/torsten_def.hpp>
#include <stan/math/torsten/pk_ode_integrator.hpp>
#include <stan/math/torsten/pk_ss_system.hpp>
#include <stan/math/torsten/PKModel/Pred/unpromote.hpp>

namespace refactor {

  using boost::math::tools::promote_args;
  using Eigen::Matrix;
  using Eigen::Dynamic;
  using torsten::PkOdeIntegrator;
  using torsten::PkOdeIntegratorId;

  template<typename F, typename T_rate>
  struct PKOdeFunctorRateAdaptor;

  /*
   * Adaptor for ODE functor when rate is data. In this
   * case rate should be passed in by @c x_r. 
   */
  template<typename F>
  struct PKOdeFunctorRateAdaptor<F, double> {
    const F f;
    const int dummy = -1;
    explicit PKOdeFunctorRateAdaptor(const F& f0) : f(f0) {}
    PKOdeFunctorRateAdaptor(const F& f0, const int i) : f(f0), dummy(i) {}

    /*
     * @c algebra_solver requires a default constructor for
     * the functor type passed as its 1st argument. In this
     * case the ODE functor @c F has default constructor, so
     * we add a default value of @c dummy, as it is not
     * relevant when @c rate is data anyway.
     */
    PKOdeFunctorRateAdaptor() {}

    /*
     * Evaluate ODE functor and add @c rate data afterwards.
     */
    template <typename T0, typename T1, typename T2>
    inline std::vector<typename stan::return_type<T1, T2>::type>
    operator()(const T0& t,
               const std::vector<T1>& y,
               const std::vector<T2>& theta,
               const std::vector<double>& x_r,
               const std::vector<int>& x_i,
               std::ostream* msgs) const {
      std::vector<typename stan::return_type<T1, T2>::type> res;
      res = f(t, y, theta, x_r, x_i, msgs);
      for (size_t i = 0; i < y.size(); i++) res.at(i) += x_r.at(i);
      return res;
    }
  };

  /*
   * Adaptor for ODE functor when rate is @c var. In this
   * case rate should be passed in by @c theta. When
   * constructed, this type stores an index @c index_rate of @c rate
   * within @c theta. Note that we only allow @c T_rate to
   * be @c var when @c T_par is @c var, so in this case @c
   * theta is always a @c var vector.
   */
  template<typename F>
  struct PKOdeFunctorRateAdaptor<F, stan::math::var> {
    const F f;
    const int index_rate = -1;
    PKOdeFunctorRateAdaptor(const F& f0, int i) : f(f0), index_rate(i) {}

    /*
     * @c algebra_solver requires a default constructor. In
     * this case, the index of @c rate in @c theta is
     * calculated by assumption that @c theta is laid out as
     * @c [theta, rate]
     */
    PKOdeFunctorRateAdaptor() {}

    /*
     * Evaluate ODE functor, with @c theta contains original
     * parameter vector followed by @c rate params. So the
     * @c theta passed in must be modfified to reflect this
     * data arrangement.
     */
    template <typename T0, typename T1, typename T2>
    inline std::vector<typename stan::return_type<T1, T2>::type>
    operator()(const T0& t,
               const std::vector<T1>& y,
               const std::vector<T2>& theta,
               const std::vector<double>& x_r,
               const std::vector<int>& x_i,
               std::ostream* msgs) const {
      std::vector<typename stan::return_type<T1, T2>::type> res;
      res = f(t, y, theta, x_r, x_i, msgs);
      if (index_rate == -1) {
        for (size_t i = 0; i < y.size(); i++) res.at(i) += theta.at(i + theta.size() - res.size()); // NOLINT
      } else {
        for (size_t i = 0; i < y.size(); i++) res.at(i) += theta.at(i + index_rate); // NOLINT
      }

      return res;
    }
  };


  /**
   * ODE-based PKPD models.
   *
   * @tparam T_time t type
   * @tparam T_init initial condition type
   * @tparam T_rate dosing rate type
   * @tparam T_par PK parameters type
   * @tparam F ODE functor
   * @tparam Ti ODE additional parameter type, usually the ODE size
   */
  template<typename T_time, typename T_init, typename T_rate, typename T_par, typename F, typename Ti> // NOLINT
  class PKODEModel {
    const T_time &t0_;
    const Eigen::Matrix<T_init, 1, Eigen::Dynamic>& y0_;
    const std::vector<T_rate> &rate_;
    const std::vector<T_par> &par_;
    const F &f_;
    const PKOdeFunctorRateAdaptor<F, T_rate> f1;
    const int ncmt_;
  public:
    using scalar_type = typename promote_args<T_time, T_rate, T_par, T_init>::type; // NOLINT
    using aug_par_type = typename promote_args<T_rate, T_par, T_init>::type;
    using init_type   = T_init;
    using time_type   = T_time;
    using par_type    = T_par;
    using rate_type   = T_rate;
    using f_type      = F;

    /*
     * FIX ME: we need to get rid of @c ModelParameters in
     * @c Pred2
     */
    template<typename T0, typename T1, typename T2, typename T3>
    static std::vector<T1>
    get_param(const torsten::ModelParameters<T0, T1, T2, T3>& p) {
      return p.get_RealParameters();
    }

    /**
     * Constructor
     * FIXME need to remove parameter as this is for linode only.
     *
     * @tparam T_mp parameters class
     * @tparam Ts parameter types
     * @param t0 initial time
     * @param y0 initial condition
     * @param rate dosing rate
     * @param par model parameters
     * @param parameter ModelParameter type
     * @param f ODE functor
     * @param ncmt the ODE size.
     */
    template<template<typename...> class T_mp, typename... Ts>
    PKODEModel(const T_time& t0,
               const Eigen::Matrix<T_init, 1, Eigen::Dynamic>& y0,
               const std::vector<T_rate> &rate,
               const std::vector<T_par> &par,
               const T_mp<Ts...> &parameter,
               const F& f,
               const Ti &ncmt) :
      t0_(t0), y0_(y0), rate_(rate), par_(par), f_(f), f1(f_, par_.size()), ncmt_(ncmt) // NOLINT
    {}

    /**
     * Constructor
     *
     * @param t0 initial time
     * @param y0 initial condition
     * @param rate dosing rate
     * @param par model parameters
     * @param f ODE functor
     * @param ncmt the ODE size.
     */
    PKODEModel(const T_time& t0,
               const Eigen::Matrix<T_init, 1, Eigen::Dynamic>& y0,
               const std::vector<T_rate> &rate,
               const std::vector<T_par> &par,
               const F& f,
               const Ti &ncmt) :
      t0_(t0), y0_(y0), rate_(rate), par_(par), f_(f), f1(f_, par_.size()), ncmt_(ncmt) // NOLINT
    {}

    /**
     * Constructor from any other model type as long as it
     * provides enough information to build an ODE model.
     *
     * @tparam T_model model type that can be used as ODE model.
     * @tparam Ts type parameters for @c T_model
     * @param m model that provides ODE information.
     */
    template<template<typename...> class T_model, typename... Ts>    
    PKODEModel(const T_model<Ts...>& m) :
      t0_(m.t0()),
      y0_(m.y0()),
      rate_(m.rate()),
      par_(m.par()),
      f_(m.f()),
      f1(f_, par_.size()),
      ncmt_(m.ncmt())
    {}
    
    /**
     * @return initial time
     */
    const T_time              & t0()       const { return t0_; }
    /**
     * @return initial condition
     */
    const PKRec<T_init>    & y0()       const { return y0_; }
    /**
     * @return dosing rate
     */
    const std::vector<T_rate> & rate()     const { return rate_; }
    /**
     * @return model parameters
     */
    const std::vector<T_par>  & par()      const { return par_; }
    /**
     * @return RHS functor
     */
    const F                   & f ()      const { return f_; }
    /**
     * @return ODE size
     */
    const int                 & ncmt ()    const { return ncmt_; }

  private:
    /*
     * We overload @c integrate so that we can pass @c rate
     * with different types, due to limit of c++ of partial
     * spec of member functions. The first version is for
     * rate being data.
     */
    template<PkOdeIntegratorId It>
    Eigen::Matrix<scalar_type, Eigen::Dynamic, 1>
    integrate(const std::vector<double> &rate,
              const T_time& dt,
              const double rtol,
              const double atol,
              const int max_num_steps,
              std::ostream* msgs) const {
      using stan::math::value_of;

      const double t0 = value_of(t0_);
      std::vector<T_time> ts{t0_ + dt};
      Eigen::Matrix<scalar_type, Eigen::Dynamic, 1> res;
      if (ts[0] == t0_) {
        res = y0_;
      } else {
        auto y = stan::math::to_array_1d(y0_);
        std::vector<int> x_i;
        std::vector<std::vector<scalar_type> > res_v =
          PkOdeIntegrator<It>(rtol, atol, max_num_steps, msgs)(f1, y, t0, ts, par_, rate, x_i); // NOLINT
        res = stan::math::to_vector(res_v[0]);
      }
      return res;
    }

    /*
     * For steady-state. Same as the first version but with given @c init.
     */
    template<PkOdeIntegratorId It>
    Eigen::Matrix<double, Eigen::Dynamic, 1>
    integrate(const std::vector<double> &rate,
              const Eigen::Matrix<double, 1, Eigen::Dynamic>& y0,
              const double& dt,
              const double rtol,
              const double atol,
              const int max_num_steps,
              std::ostream* msgs) const {
      using stan::math::value_of;

      const double t0 = value_of(t0_) - dt;
      std::vector<double> ts{value_of(t0_)};
      Eigen::Matrix<double, Eigen::Dynamic, 1> res;
      if (ts[0] == t0) {
        res = y0;
      } else {
        auto y = stan::math::to_array_1d(y0);
        PKOdeFunctorRateAdaptor<F, double> f(f_);
        std::vector<int> x_i;
        const std::vector<double> pars{value_of(par_)};
        std::vector<std::vector<double> > res_v =
          PkOdeIntegrator<It>(rtol, atol, max_num_steps, msgs)(f, y, t0, ts, pars, rate, x_i); // NOLINT
        res = stan::math::to_vector(res_v[0]);
      }
      return res;
    }

    /*
     * When @c rate is a parameter, we append it to @c theta.
     */
    template<PkOdeIntegratorId It>
    Eigen::Matrix<scalar_type, Eigen::Dynamic, 1>
    integrate(const std::vector<stan::math::var> &rate,
              const T_time& dt,
              const double rtol,
              const double atol,
              const int max_num_steps,
              std::ostream* msgs) const {
      using stan::math::var;
      using stan::math::value_of;

      const double t0 = value_of(t0_);
      std::vector<T_time> ts{t0_ + dt};
      Eigen::Matrix<scalar_type, Eigen::Dynamic, 1> res;
      if (ts[0] == t0_) {
        res = y0_;
      } else {
        auto y = stan::math::to_array_1d(y0_);
        std::vector<double> x_r;
        std::vector<int> x_i;
        std::vector<stan::math::var> theta(par_.size() + rate.size());
        for (size_t i = 0; i < par_.size(); ++i) theta[i] = par_[i];
        for (size_t i = 0; i < rate.size(); ++i) theta[i + par_.size()] = rate[i];
        std::vector<std::vector<scalar_type> > res_v =
          PkOdeIntegrator<It>(rtol, atol, max_num_steps, msgs)(f1, y, t0, ts, theta, x_r, x_i); // NOLINT
        res = stan::math::to_vector(res_v[0]);
      }
      return res;
    }

  public:
    /*
     * solve ODE system. The different cases when @c rate is
     * @c var or data are handled by private methods @c integrate.
     *
     * @parm[in] dt next time point when result is to be *
     *           solved. The actual time point will be @c t0+dt
     * @parm[in] rtol relative tolerance for ODE solver
     * @parm[in] atol absolute tolerance for ODE solver
     * @parm[in] max_num_steps max number of steps between @c t0 and @c t0+dt.
     * @parm[in] msgs output stream.
     * @return an Eigen matrix with each row being the
     *         solution at certain time step. Hence the returned
     *         matrix is of dim (numer of time steps) x (siez of ODE system).
     */
    template<PkOdeIntegratorId It>
    Eigen::Matrix<scalar_type, Eigen::Dynamic, 1>
    solve(const T_time& dt,
          const double rtol,
          const double atol,
          const int max_num_steps,
          std::ostream* msgs) const {
      return integrate<It>(rate_, dt, rtol, atol, max_num_steps, msgs);
    }

    /**
     * Solve steady state ODE-based PKPD model with @c amt
     * and @c rate both being data.
     *
     * @tparam It ODE numerical integrator ID.
     * @tparam T_ii dosing interval type
     * @param amt dosing amount
     * @param rate dosing rate
     * @param ii dosing interval
     * @param cmt compartment where the dosing occurs
     * @param rtol relative tolerance of ODE integrator.
     * @param atol absolute tolerance of ODE integrator.
     * @param max_num_steps max number of steps of ODE integrator.
     * @param msgs @c ostream output of ODE integrator.
     * @return col vector of the steady state ODE solution
     */
    template<PkOdeIntegratorId It, typename T_ii>
    Eigen::Matrix<typename promote_args<T_ii, par_type>::type, Eigen::Dynamic, 1> // NOLINT
    solve(const double& amt,
          const double& rate,
          const T_ii& ii,
          const int& cmt,
          const double rtol,
          const double atol,
          const int max_num_steps,
          std::ostream* msgs) const {
      using Eigen::Matrix;
      using Eigen::Dynamic;
      using Eigen::VectorXd;
      using std::vector;
      using stan::math::algebra_solver;
      using stan::math::to_vector;
      using stan::math::value_of;
      using refactor::PKRec;

      typedef typename promote_args<T_ii, T_par>::type scalar;

      const int& ncmt = ncmt_;

      Matrix<scalar, Dynamic, 1> pred;

      // Arguments for ODE integrator (and initial guess)
      double ii_dbl = torsten::unpromote(ii);
      Matrix<double, 1, Dynamic> init_dbl(ncmt);
      for (int i = 0; i < ncmt; i++) init_dbl(i) = 0;
      vector<double> x_r(ncmt, 0);
      vector<int> x_i(0);

      // Arguments for algebraic solver
      Matrix<double, Dynamic, 1> y;
      const double alge_rtol = 1e-10;  // default
      const double f_tol = 1e-4;  // empirical
      const long int alge_max_steps = 1e3;  // default

      using F_ss = PKOdeFunctorRateAdaptor<F, double>;
      PkOdeIntegrator<It> integrator(rtol, atol, max_num_steps, msgs);
      torsten::SSFunctor<It, double, double, F_ss, void> system(f1, ii_dbl, cmt, integrator); // NOLINT

      if (rate == 0) {  // bolus dose
        // compute initial guess
        init_dbl(cmt - 1) = amt;
        y = integrate<It>(x_r, init_dbl, ii_dbl, rtol, atol, max_num_steps, msgs); // NOLINT
        x_r.push_back(amt);
        pred = algebra_solver(system, y,
                              to_vector(par_),
                              x_r, x_i,
                              0, alge_rtol, f_tol, alge_max_steps);
        // DEV - what tuning parameters should we use for the algebra solver?
        // DEV - update initial guess or tuning parameters if result not good?
      }  else if (ii > 0) {  // multiple truncated infusions
        x_r[cmt - 1] = rate;
        // compute initial guess
        y = integrate<It>(x_r, init_dbl, ii_dbl, rtol, atol, max_num_steps, msgs); // NOLINT
        x_r.push_back(amt);
        pred = algebra_solver(system, y,
                              to_vector(par_),
                              x_r, x_i,
                              0, alge_rtol, 1e-3, alge_max_steps);  // FIX ME
        // use ftol
      } else {  // constant infusion
        x_r[cmt - 1] = rate;
        y = integrate<It>(x_r, init_dbl, 100.0, rtol, atol, max_num_steps, msgs); // NOLINT
        x_r.push_back(amt);
        pred = algebra_solver(system, y,
                              to_vector(par_),
                              x_r, x_i,
                              0, alge_rtol, f_tol, alge_max_steps);
      }
      return pred;
    }

    /**
     * Solve steady state ODE-based PKPD model with @c amt
     * being param and @c rate both being data.
     *
     * @tparam It ODE numerical integrator ID.
     * @tparam T_amt dosing amount type.
     * @tparam T_ii dosing interval type
     * @param amt dosing amount
     * @param rate dosing rate
     * @param ii dosing interval
     * @param cmt compartment where the dosing occurs
     * @param rtol relative tolerance of ODE integrator.
     * @param atol absolute tolerance of ODE integrator.
     * @param max_num_steps max number of steps of ODE integrator.
     * @param msgs @c ostream output of ODE integrator.
     * @return col vector of the steady state ODE solution
     */
    template<PkOdeIntegratorId It, typename T_amt, typename T_ii>
    Eigen::Matrix<typename promote_args<T_amt, T_ii, T_par>::type, Eigen::Dynamic, 1> // NOLINT
    solve(const T_amt& amt,
          const double& rate,
          const T_ii& ii,
          const int& cmt,
          const double rtol,
          const double atol,
          const int max_num_steps,
          std::ostream* msgs) const {
      using Eigen::Matrix;
      using Eigen::Dynamic;
      using Eigen::VectorXd;
      using std::vector;
      using stan::math::algebra_solver;
      using stan::math::to_vector;
      using stan::math::invalid_argument;

      typedef typename promote_args<T_ii, T_amt, T_time, T_par>::type scalar;

      const int& ncmt = ncmt_;
      const std::vector<T_par>& pars = par_;

      Matrix<scalar, Dynamic, 1> pred;

      // Arguments for the ODE integrator
      double ii_dbl = torsten::unpromote(ii);
      Matrix<double, 1, Dynamic> init_dbl(ncmt);
      for (int i = 0; i < ncmt; i++) init_dbl(i) = 0;
      vector<double> x_r(ncmt, 0);
      vector<int> x_i(0);

      // Arguments for algebraic solver
      Matrix<double, Dynamic, 1> y;
      const double alge_rtol = 1e-10;  // default
      const double f_tol = 5e-4;  // empirical (note: differs from other function)
      const long int alge_max_steps = 1e4;  // default  // NOLINT

      // construct algebraic function
      torsten::general_functor<F> f1(f_);
      using F_ss = PKOdeFunctorRateAdaptor<F, double>;
      PkOdeIntegrator<It> integrator(rtol, atol, max_num_steps, msgs);
      torsten::SSFunctor<It, T_amt, double, F_ss, void> system(F_ss(f_), ii_dbl, cmt, integrator); // NOLINT

      int npar = pars.size();
      Matrix<scalar, Dynamic, 1> parms(npar + 1);
      for (int i = 0; i < npar; i++) parms(i) = pars[i];
      parms(npar) = amt;

      if (rate == 0) {  // bolus dose
        // compute initial guess
        init_dbl(cmt - 1) = torsten::unpromote(amt);
        y = integrate<It>(x_r, init_dbl, ii_dbl, rtol, atol, max_num_steps, msgs); // NOLINT
        pred = algebra_solver(system, y, parms, x_r, x_i, 0, alge_rtol, f_tol, alge_max_steps); // NOLINT
      }  else if (ii > 0) {  // multiple truncated infusions
        // compute initial guess
        x_r[cmt - 1] = rate;
        y = integrate<It>(x_r, init_dbl, ii_dbl, rtol, atol, max_num_steps, msgs); // NOLINT
        pred = algebra_solver(system, y, parms, x_r, x_i, 0, alge_rtol, 1e-3, alge_max_steps); // NOLINT
      } else {  // constant infusion
        x_r[cmt - 1] = rate;
        y = integrate<It>(x_r, init_dbl, 100.0, rtol, atol, max_num_steps, msgs); // NOLINT
        pred = algebra_solver(system, y, parms, x_r, x_i, 0, alge_rtol, f_tol, alge_max_steps); // NOLINT
      }
      return pred;
    }
  };
}

#endif
