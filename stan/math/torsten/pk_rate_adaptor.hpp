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
