#ifndef STAN_MATH_TORSTEN_REFACTOR_ODE_MODEL_HPP
#define STAN_MATH_TORSTEN_REFACTOR_ODE_MODEL_HPP

#include <stan/math/rev/mat/functor/algebra_solver.hpp>
#include <stan/math/torsten/model_solve_d.hpp>
#include <stan/math/torsten/pk_nvars.hpp>
#include <stan/math/torsten/pk_ss_system.hpp>
#include <stan/math/torsten/PKModel/Pred/unpromote.hpp>

namespace refactor {

  using boost::math::tools::promote_args;
  using Eigen::Matrix;
  using Eigen::Dynamic;
  using torsten::PMXOdeIntegrator;
  using torsten::PMXOdeIntegratorId;

  /*
   * when solving ODE with @c var rate, we append it to
   * parameter vector. Note that spurious @c var
   * parameters will be generated if the original parameters
   * are data.
   */
  template<typename T>
  inline std::vector<stan::math::var> parameters_with_rate(const std::vector<T> &par,
                                                           const std::vector<stan::math::var> &rate)
  {
    std::vector<stan::math::var> theta(par.size() + rate.size());
    for (size_t i = 0; i < par.size();  ++i) theta[i] = par[i];
    for (size_t i = 0; i < rate.size(); ++i) theta[i + par.size()] = rate[i];
    return theta;
  }

  template<typename T>
  inline std::vector<stan::math::var> parameters_with_rate(const std::vector<T> &par,
                                                           const std::vector<double> &rate)
  {
    std::vector<stan::math::var> res;
    return res;
  }

  template<typename F, typename T_rate>
  struct PMXOdeFunctorRateAdaptor;

  /*
   * Adaptor for ODE functor when rate is data. In this
   * case rate should be passed in by @c x_r. 
   */
  template<typename F>
  struct PMXOdeFunctorRateAdaptor<F, double> {
    const F f;
    const int dummy = -1;
    explicit PMXOdeFunctorRateAdaptor(const F& f0) : f(f0) {}
    PMXOdeFunctorRateAdaptor(const F& f0, const int i) : f(f0), dummy(i) {}

    /*
     * @c algebra_solver requires a default constructor for
     * the functor type passed as its 1st argument. In this
     * case the ODE functor @c F has default constructor, so
     * we add a default value of @c dummy, as it is not
     * relevant when @c rate is data anyway.
     */
    PMXOdeFunctorRateAdaptor() : f() {}

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
  struct PMXOdeFunctorRateAdaptor<F, stan::math::var> {
    const F f;
    const int index_rate = -1;
    PMXOdeFunctorRateAdaptor(const F& f0, int i) : f(f0), index_rate(i) {}

    /*
     * @c algebra_solver requires a default constructor. In
     * this case, the index of @c rate in @c theta is
     * calculated by assumption that @c theta is laid out as
     * @c [theta, rate]
     */
    PMXOdeFunctorRateAdaptor() {}

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
  template<typename T_time, typename T_init, typename T_rate, typename T_par, typename F> // NOLINT
  class PKODEModel {
    const T_time &t0_;
    const Eigen::Matrix<T_init, 1, Eigen::Dynamic>& y0_;
    const std::vector<T_rate> &rate_;
    const std::vector<T_par> &par_;
    const F &f_;
    const PMXOdeFunctorRateAdaptor<F, T_rate> f1;
    const int ncmt_;
  public:
    using scalar_type = typename promote_args<T_time, T_rate, T_par, T_init>::type; // NOLINT
    using aug_par_type = typename promote_args<T_rate, T_par, T_init>::type;
    using init_type   = T_init;
    using time_type   = T_time;
    using par_type    = T_par;
    using rate_type   = T_rate;
    using f_type      = F;

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
               const F& f) :
      t0_(t0), y0_(y0), rate_(rate), par_(par), f_(f), f1(f_, par_.size()), ncmt_(y0.size()) // NOLINT
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
               const F& f) :
      t0_(t0), y0_(y0), rate_(rate), par_(par), f_(f), f1(f_, par_.size()), ncmt_(y0.size()) // NOLINT
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
    
    /*
     * The number of @c var that will be in the ODE
     * integrator. Note that not all @c var will be
     * explicitly in an integrator's call signature.
     * Since we only step one time-step, if @c t0_ is @c var
     * it only adds one(because @c ts will be of size one). Also note that regardless if @c
     * par_ is @c var or not, when @c rate_ is @c var we
     * generate a new @c theta of @c var vector to pass to ODE integrator that
     * contains both @c par_ and @c rate_.
     */

    /*
     * calculate number of @c vars for transient dosing.
     */
    static int nvars(int ncmt, int npar) {
      using stan::is_var;
      int n = 0;
      if (is_var<T_time>::value) n++; // t0
      if (is_var<T_init>::value) n += ncmt;
      if (is_var<T_rate>::value) {
        n += ncmt + npar;
      } else if (is_var<T_par>::value) {
        n += npar;
      }
      return n;
    }

    /*
     * calculate number of @c vars for steady-state dosing.
     */
    template<typename T_a, typename T_r, typename T_ii>
    static int nvars(int npar) {
      using stan::is_var;
      int n = 0;
      if (is_var<T_a>::value) n++; // amt
      if (is_var<T_r>::value) n++; // rate
      if (is_var<T_ii>::value) n++; // ii
      if (is_var<T_par>::value) n += npar;
      return n;
    }

    template<typename T0, typename T1, typename T2, typename T3>
    static int nvars(const T0& t0,
                     const Eigen::Matrix<T1, 1, Eigen::Dynamic>& y0,
                     const std::vector<T2> &rate,
                     const std::vector<T3> &par) {
      using stan::is_var;
      int res = 0;
      if (is_var<T0>::value) res += 1;
      if (is_var<T1>::value) res += y0.size();
      if (is_var<T2>::value) {
        res += rate.size() + par.size();
      } else if (is_var<T3>::value) {
        res += par.size();
      }
      return res;
    }

    /*
     * generate @c var vector that consisting of all the
     * parameters used passed to an ODE integrator. The output
     * will be arranged in order (y0, par, rate, dt)
     */
    template<typename T0, typename T1, typename T2, typename T3>
    static std::vector<stan::math::var> vars(const T0& t1,
                                             const Eigen::Matrix<T1, 1, Eigen::Dynamic>& y0,
                                             const std::vector<T2> &rate,
                                             const std::vector<T3> &par) {
      using stan::is_var;      
      using stan::math::var;
      std::vector<stan::math::var> res(nvars(t1, y0, rate, par));
      int i = 0;
      if (is_var<T1>::value) {
        for (int j = 0; j < y0.size(); ++j) {
          res[i] = y0(j);
          i++;
        }
      }
      if (is_var<T2>::value) {
        std::vector<stan::math::var> theta(parameters_with_rate(par, rate));
        for (size_t j = 0; j < theta.size(); ++j) {
          res[i] = theta[j];
          i++;
        }
      } else if (is_var<T3>::value) {
        for (size_t j = 0; j < par.size(); ++j) {
          res[i] = par[j];
          i++;
        }
      }

      if (is_var<T0>::value) {
        res[i] = t1;
      }
      return res;
    }

    /*
     * Calculate the size of the entire system, with ODe
     * solutions and their gradients w.r.t. the parameters.
     */
    template<typename T0, typename T1, typename T2, typename T3>
    static int n_sys(const T0& t0,
                     const Eigen::Matrix<T1, 1, Eigen::Dynamic>& y0,
                     const std::vector<T2> &rate,
                     const std::vector<T3> &par) {
      using stan::is_var;
      int n = nvars(t0, y0, rate, par);
      return y0.size() * (n + 1);
    }

    /*
     * calculate number of vars with constructed data
     */
    template<typename T0>
    int nvars(const T0& t1) {
      return nvars(t1, y0_, rate_, par_);
    }

    /*
     * return the number @c var that will be the parameters
     * of the stead-state dosing event's solution
     */
    template<typename T_a, typename T_r, typename T_ii>
    int nvars(const T_a& a, const T_r& r, const T_ii& ii) {
      return torsten::pk_nvars(a, r, ii, par_);
    }

    /*
     * return @c vars that will be in integrator calls
     */
    template<typename T0>
    std::vector<stan::math::var> vars(const T0 t1) {
      return vars(t1, y0_, rate_, par_);
    }

    /*
     * return @c vars that will be steady-state
     * solution. For SS solution @c rate_ or @ y0_ will not
     * be in the solution.
     */
    template<typename T_a, typename T_r, typename T_ii>
    std::vector<stan::math::var> vars(const T_a& a, const T_r& r, const T_ii& ii) {
      return torsten::dsolve::pk_vars(a, r, ii, par_);
    }

    /*
     * calculate size of the entire system
     */
    int n_sys() const {
      return n_sys(t0_, y0_, rate_, par_);
    }

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
     * When time is data, the time step is trivial
     */
    std::vector<double> time_step(const double t_next) const {
      return {t_next};
    }

    /*
     * When time is param, the time step needs to
     * incorporate information of initial time.
     */
    std::vector<stan::math::var> time_step(const stan::math::var& t_next) const {
      return {stan::math::value_of(t0_) + t_next - t0_};
    }

    /*
     * We overload @c integrate so that we can pass @c rate
     * with different types, due to limit of c++ of partial
     * spec of member functions. The first version is for
     * rate being data.
     */
    template<PMXOdeIntegratorId It>
    Eigen::Matrix<scalar_type, Eigen::Dynamic, 1>
    integrate(const std::vector<double> &rate,
              const T_time& t_next,
              const PMXOdeIntegrator<It>& integrator) const {
      using stan::math::value_of;

      const double t0 = value_of(t0_);
      std::vector<T_time> ts(time_step(t_next));
      Eigen::Matrix<scalar_type, Eigen::Dynamic, 1> res;
      if (ts[0] == t0_) {
        res = y0_;
      } else {
        auto y = stan::math::to_array_1d(y0_);
        std::vector<int> x_i;
        std::vector<std::vector<scalar_type> > res_v =
          integrator(f1, y, t0, ts, par_, rate, x_i); // NOLINT
        res = stan::math::to_vector(res_v[0]);
      }
      return res;
    }

    /*
     * For steady-state. Same as the first version but with given @c init.
     */
    template<PMXOdeIntegratorId It>
    Eigen::Matrix<double, Eigen::Dynamic, 1>
    integrate(const std::vector<double> &rate,
              const Eigen::Matrix<double, 1, Eigen::Dynamic>& y0,
              const double& dt,
              const PMXOdeIntegrator<It>& integrator) const {
      using stan::math::value_of;

      const double t0 = value_of(t0_) - dt;
      std::vector<double> ts{value_of(t0_)};
      Eigen::Matrix<double, Eigen::Dynamic, 1> res;
      if (ts[0] == t0) {
        res = y0;
      } else {
        auto y = stan::math::to_array_1d(y0);
        PMXOdeFunctorRateAdaptor<F, double> f(f_);
        std::vector<int> x_i;
        const std::vector<double> pars{value_of(par_)};
        std::vector<std::vector<double> > res_v =
          integrator(f, y, t0, ts, pars, rate, x_i); // NOLINT
        res = stan::math::to_vector(res_v[0]);
      }
      return res;
    }

    /*
     * When @c rate is a parameter, we append it to @c theta.
     */
    template<PMXOdeIntegratorId It>
    Eigen::Matrix<scalar_type, Eigen::Dynamic, 1>
    integrate(const std::vector<stan::math::var> &rate,
              const T_time& t_next,
              const PMXOdeIntegrator<It>& integrator) const {
      using stan::math::var;
      using stan::math::value_of;

      const double t0 = value_of(t0_);
      std::vector<T_time> ts(time_step(t_next));
      Eigen::Matrix<scalar_type, Eigen::Dynamic, 1> res;
      if (ts[0] == t0_) {
        res = y0_;
      } else {
        auto y = stan::math::to_array_1d(y0_);
        std::vector<double> x_r;
        std::vector<int> x_i;
        std::vector<stan::math::var> theta(parameters_with_rate(par_, rate));
        std::vector<std::vector<scalar_type> > res_v =
          integrator(f1, y, t0, ts, theta, x_r, x_i); // NOLINT
        res = stan::math::to_vector(res_v[0]);
      }
      return res;
    }

    /*
     * Solve the ODE but return the results in data
     * consisting of solution values and gradients.
     */
    template<PMXOdeIntegratorId It,
             typename std::enable_if_t<It == torsten::PkBdf || It == torsten::PkAdams || It == torsten::PkRk45>* = nullptr>
    Eigen::VectorXd integrate_d(const std::vector<stan::math::var> &rate,
                                const T_time& t_next,
                                const PMXOdeIntegrator<It>& integrator) const {
      using stan::math::var;
      using stan::math::value_of;

      const double t0 = value_of(t0_);
      std::vector<T_time> ts(time_step(t_next));
      Eigen::VectorXd res(n_sys());
      if (ts[0] == t0_) {
        res = stan::math::value_of(y0_);
        res.resize(n_sys());
      } else {
        auto y = stan::math::to_array_1d(y0_);
        std::vector<double> x_r;
        std::vector<int> x_i;
        std::vector<stan::math::var> theta(parameters_with_rate(par_, rate));
        res = integrator.solve_d(f1, y, t0, ts, theta, x_r, x_i).col(0);
      }
      return res;
    }

    /*
     * Solve the ODE but return the results in data
     * consisting of solution values and gradients.
     */
    template<PMXOdeIntegratorId It,
             typename std::enable_if_t<It == torsten::PkBdf || It == torsten::PkAdams || It == torsten::PkRk45>* = nullptr>
    Eigen::VectorXd integrate_d(const std::vector<double> &rate,
                                const T_time& t_next,
                                const PMXOdeIntegrator<It>& integrator) const {
      using stan::math::value_of;

      Eigen::VectorXd res(n_sys());
      const double t0 = value_of(t0_);
      std::vector<T_time> ts(time_step(t_next));
      if (ts[0] == t0_) {

      } else {
        auto y = stan::math::to_array_1d(y0_);
        std::vector<int> x_i;
        res = integrator.solve_d(f1, y, t0, ts, par_, rate, x_i).col(0);
      }
      return res;
    }

  public:
    /*
     * solve ODE system. The different cases when @c rate is
     * @c var or data are handled by private methods @c integrate.
     *
     * @parm[in] t_next next time point when result is to be
     *           solved. The actual time point will be @c t0+dt
     * @parm[in] rtol relative tolerance for ODE solver
     * @parm[in] atol absolute tolerance for ODE solver
     * @parm[in] max_num_steps max number of steps between @c t0 and @c t0+dt.
     * @parm[in] msgs output stream.
     * @return an Eigen matrix with each row being the
     *         solution at certain time step. Hence the returned
     *         matrix is of dim (numer of time steps) x (siez of ODE system).
     */
    template<PMXOdeIntegratorId It>
    Eigen::Matrix<scalar_type, Eigen::Dynamic, 1>
    solve(const T_time& t_next,
          const PMXOdeIntegrator<It>& integrator) const {
      return integrate(rate_, t_next, integrator);
    }

    /*
     * Solve the ODE but return the results in form of data.
     * The default behavior is defined in function template
     * @c model_solve_d(), using autodiff to recalculate gradients.
     * Some integrators, such as @c PkBdf, have their own implementation
     * that can return data directly, so we skip them.
     */
    template<PMXOdeIntegratorId It,
             typename std::enable_if_t<It == torsten::StanBdf || It == torsten::StanAdams || It == torsten::StanRk45>* = nullptr>
    Eigen::VectorXd solve_d(const T_time& t_next,
                            const PMXOdeIntegrator<It>& integrator) const
    {
      static const char* caller = "PKODEModel::solve";
      stan::math::check_greater(caller, "time step", t_next, t0_);

      return torsten::model_solve_d(*this, t_next, integrator);
    }

    /*
     * @c PkBdf can return results in form of data directly,
     * thanks to @c pk_cvodes_integrator implementation.
     */
    template<PMXOdeIntegratorId It,
             typename std::enable_if_t<It == torsten::PkBdf || It == torsten::PkAdams || It == torsten::PkRk45>* = nullptr>
    Eigen::VectorXd solve_d(const T_time& t_next,
                            const PMXOdeIntegrator<It>& integrator) const {
      static const char* caller = "PMXOdeModel::solve_d";
      stan::math::check_greater(caller, "next time", t_next, t0_);

      return integrate_d(rate_, t_next, integrator);
    }

    /**
     * Solve steady state ODE-based PKPD model with @c amt
     * and @c rate both being data. Specifically for
     * constant rate we solve ODEs of the form
     *
     * x′(t)=F(x,t,θ)+ Rate,
     *
     * and solution x satisfies
     *
     * F(x,t,θ) + Rate = 0.
     *
     * For the periodic steady-state resulting from a
     * periodic input D into one or more cpts at a constant
     * interval τ, the root finder is combined with the ODE
     * solver to solve
     *
     * x(t0 + τ) − x(t0) = 0.
     *
     * Denote x(t0) by x0, it is equivalent to solve x0 in
     *
     * x(t0 + τ) = x0, where in the interval (t0, t0 + τ)
     * x(t) satisfies the ODE with x0 as init condition.
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
    template<PMXOdeIntegratorId It, typename T_ii>
    Eigen::Matrix<typename promote_args<T_ii, par_type>::type, Eigen::Dynamic, 1> // NOLINT
    solve(const double& amt,
          const double& rate,
          const T_ii& ii,
          const int& cmt,
          const PMXOdeIntegrator<It>& integrator) const {
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

      using F_ss = PMXOdeFunctorRateAdaptor<F, double>;
      torsten::SSFunctor<It, double, double, F_ss, void> system(f1, ii_dbl, cmt, integrator); // NOLINT

      if (rate == 0) {          // multiple bolus dose: AMT>0, RATE=0, SS=1, and II>0.
        // compute initial guess
        init_dbl(cmt - 1) = amt;
        y = integrate(x_r, init_dbl, ii_dbl, integrator); // NOLINT
        x_r.push_back(amt);
        pred = algebra_solver(system, y,
                              to_vector(par_),
                              x_r, x_i,
                              0, alge_rtol, f_tol, alge_max_steps);
        // DEV - what tuning parameters should we use for the algebra solver?
        // DEV - update initial guess or tuning parameters if result not good?
      }  else if (ii > 0) {     // multiple infusion: AMT>0, RATE>0, SS=1, and II>0.
        x_r[cmt - 1] = rate;
        // compute initial guess
        y = integrate(x_r, init_dbl, ii_dbl, integrator); // NOLINT
        x_r.push_back(amt);
        pred = algebra_solver(system, y,
                              to_vector(par_),
                              x_r, x_i,
                              0, alge_rtol, 1e-3, alge_max_steps);  // FIX ME
        // use ftol
      } else {                  // constant infusion: AMT=0, RATE>0, SS=1, and II=0.
        x_r[cmt - 1] = rate;
        y = integrate(x_r, init_dbl, 100.0, integrator); // NOLINT
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
    template<PMXOdeIntegratorId It, typename T_amt, typename T_ii>
    Eigen::Matrix<typename promote_args<T_amt, T_ii, T_par>::type, Eigen::Dynamic, 1> // NOLINT
    solve(const T_amt& amt,
          const double& rate,
          const T_ii& ii,
          const int& cmt,
          const PMXOdeIntegrator<It>& integrator) const {
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
      using F_ss = PMXOdeFunctorRateAdaptor<F, double>;
      torsten::SSFunctor<It, T_amt, double, F_ss, void> system(F_ss(f_), ii_dbl, cmt, integrator); // NOLINT

      int npar = pars.size();
      Matrix<scalar, Dynamic, 1> parms(npar + 1);
      for (int i = 0; i < npar; i++) parms(i) = pars[i];
      parms(npar) = amt;

      if (rate == 0) {  // bolus dose
        // compute initial guess
        init_dbl(cmt - 1) = torsten::unpromote(amt);
        y = integrate(x_r, init_dbl, ii_dbl, integrator); // NOLINT
        pred = algebra_solver(system, y, parms, x_r, x_i, 0, alge_rtol, f_tol, alge_max_steps); // NOLINT
      }  else if (ii > 0) {  // multiple truncated infusions
        // compute initial guess
        x_r[cmt - 1] = rate;
        y = integrate(x_r, init_dbl, ii_dbl, integrator); // NOLINT
        pred = algebra_solver(system, y, parms, x_r, x_i, 0, alge_rtol, 1e-3, alge_max_steps); // NOLINT
      } else {  // constant infusion
        x_r[cmt - 1] = rate;
        y = integrate(x_r, init_dbl, 100.0, integrator); // NOLINT
        pred = algebra_solver(system, y, parms, x_r, x_i, 0, alge_rtol, f_tol, alge_max_steps); // NOLINT
      }
      return pred;
    }

    /*
     * return steady state solution in form of data, use
     * default behavior, namely take gradients using autodiff.
     */
    template<PMXOdeIntegratorId It, typename T_amt, typename T_ii>
    Eigen::VectorXd solve_d(const T_amt& amt,
                            const double& rate,
                            const T_ii& ii,
                            const int& cmt,
                            const PMXOdeIntegrator<It>& integrator) const {
      return torsten::model_solve_d(*this, amt, rate, ii, cmt, integrator);
    }

  };

}

#endif
