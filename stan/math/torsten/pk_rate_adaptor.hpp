#ifndef STAN_MATH_TORSTEN_RATE_ADAPTOR_HPP
#define STAN_MATH_TORSTEN_RATE_ADAPTOR_HPP

#include <stan/math/rev/core.hpp>
#include <stan/math/fwd/core.hpp>
#include <vector>
#include <iostream>

namespace refactor {

  using boost::math::tools::promote_args;
  using refactor::PKODEModel;

  template<typename F, typename T_rate, typename T_par>
  struct PKOdeFunctorRateAdaptor;

  /*
   * Adaptor for ODE functor when rate is data. In this
   * case rate should be passed in by @c x_r.
   */
  template<typename F, typename T_par>  
  struct PKOdeFunctorRateAdaptor<F, double, T_par> {
    const F& f;
    PKOdeFunctorRateAdaptor(const F& f0) : f(f0) {}

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
  struct PKOdeFunctorRateAdaptor<F, stan::math::var, stan::math::var> {
    const F& f;
    const int index_rate;
    PKOdeFunctorRateAdaptor(const F& f0, int i) : f(f0), index_rate(i) {}

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
      for (size_t i = 0; i < y.size(); i++) res.at(i) += theta.at(i + index_rate);
      return res;
    }
  };

  /*
   * Adaptor for a model in which @c rate is data.
   */
  template<typename T_model>
  class PKODERateAdaptor {
  private:
    const PKOdeFunctorRateAdaptor<torsten::f_t<T_model>,
                                  torsten::rate_t<T_model>,
                                  torsten::par_t<T_model> > f_;
    std::vector<typename stan::return_type<
                  torsten::rate_t<T_model>, torsten::par_t<T_model>>::type> theta_;
    const PKODEModel<torsten::time_t<T_model>,
                     torsten::init_t<T_model>,
                     torsten::rate_t<T_model>,
                     torsten::par_t<T_model>,
                     PKOdeFunctorRateAdaptor<torsten::f_t<T_model>,
                                             torsten::rate_t<T_model>,
                                             torsten::par_t<T_model>>, int> model_;


  public:
    template<typename T, typename std::enable_if_t<
                           !torsten::has_var_rate<T>::value>* = nullptr>
    PKODERateAdaptor(const T & pkmodel) :
      f_(pkmodel.f()),
      model_(pkmodel.t0(),
             pkmodel.y0(),
             pkmodel.rate(),
             pkmodel.par(), f_,
             pkmodel.ncmt())
    {}

    template<typename T, typename std::enable_if_t<
                           torsten::has_var_rate<T>::value>* = nullptr>
    PKODERateAdaptor(const T & pkmodel) :
      f_(pkmodel.f(), pkmodel.par().size()),
      theta_(pkmodel.par().size() + pkmodel.rate().size()),
      model_(pkmodel.t0(),
             pkmodel.y0(),
             pkmodel.rate(),
             theta_, f_,
             pkmodel.ncmt())
    {
      auto par = pkmodel.par();
      auto rate = pkmodel.rate();
      for (size_t i = 0; i < par.size(); ++i) theta_[i] = par[i]; 
      for (size_t i = 0; i < rate.size(); ++i) theta_[i + par.size()] = rate[i];
    }

    /* a rate adaptor can be used as an ODE model, as long
     * as it has the following methods.
     */
    const auto  & t0()       const { return model_.t0(); }
    const auto  & y0()       const { return model_.y0(); }
    const auto  & rate()     const { return model_.rate(); }
    const auto  & par()      const { return model_.par(); }
    const auto  & f ()       const { return model_.f(); }
    const int   & ncmt ()    const { return model_.ncmt(); }
   };

  // /*
  //  * Adaptor for a model that has @c var type rate
  //  */
  // template<typename T_model,
  //          typename std::enable_if_t<
  //            stan::math::is_var<typename T_model::rate_type>::value>* = nullptr>
  // class PKODERateAdaptor {
  //   using T_time = typename T_model::torsten::time_t<T_model>;
  //   using T_init = typename T_model::init_type;
  //   using rate_type = typename T_model::rate_type;
  //   using par_type  = typename T_model::par_type;
  //   using torsten::f_t<T_model>      = typename T_model::torsten::f_t<T_model>;
  //   using par_type = typename promote_args<par_type, rate_type>::type;
  //   using FA = torsten::ode_rate_var_functor<torsten::general_functor<torsten::f_t<T_model>> >;

  //   const std::vector<double> dummy_;
  //   const FA f_;
  //   const std::vector<par_type> theta_;
  //   const PKODEModel<T_time, T_init, double, par_type, torsten::f_t<T_model>A, int> model_;
    
  //   std::vector<par_type>
  //   concat_par_rate(const std::vector<par_type> & par,
  //                   const std::vector<rate_type> & rate) {
  //     const size_t n = par.size();
  //     std::vector<par_type> theta(n + rate.size());
  //     for (size_t i = 0; i < n; i++) theta[i] = par[i];
  //     for (size_t i = 0; i < rate.size(); i++) theta[n + i] = rate[i];
  //     return theta;
  //   }

  // public:

  //   PKODERateAdaptor(const T_model & pkmodel) :
  //     dummy_{0.0},
  //     f_(torsten::general_functor<F>(pkmodel.rhs_fun())),
  //     theta_(concat_par_rate(pkmodel.par(), pkmodel.rate())),
  //     model_(pkmodel.t0(),
  //            pkmodel.y0(),
  //            dummy_,
  //            theta_, f_,
  //            pkmodel.ncmt())
  //   {}

  //   const PKODEModel<T_time, T_init, double, par_type, FA, int>& model() {
  //     return model_;
  //   }
  // };
      
  // template<typename T_time, typename T_init, typename par_type, typename F>
  // class PKODERateAdaptor<PKODEModel<T_time, T_init, double, par_type, F, int>> {
  //   using FA = torsten::ode_rate_dbl_functor<torsten::general_functor<F> >;
  //   const torsten::ode_rate_dbl_functor<torsten::general_functor<F> > f_;
  //   const PKODEModel<T_time, T_init, double, par_type, FA, int> model_;
  // public:
  //   PKODERateAdaptor(const PKODEModel<T_time,
  //                    T_init, double, par_type, F, int> & pkmodel) :
  //     f_(torsten::general_functor<F>(pkmodel.rhs_fun())),
  //     model_(pkmodel.t0(),
  //            pkmodel.y0(),
  //            pkmodel.rate(),
  //            pkmodel.par(), f_,
  //            pkmodel.ncmt())
  //   {}

  //   const PKODEModel<T_time, T_init, double, par_type, FA, int>& model() {
  //     return model_;
  //   }
  // };
}

#endif
