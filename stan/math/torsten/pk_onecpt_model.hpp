#ifndef STAN_MATH_TORSTEN_ONECPT_MODEL_HPP
#define STAN_MATH_TORSTEN_ONECPT_MODEL_HPP

#include <stan/math/torsten/torsten_def.hpp>

namespace refactor {

  using boost::math::tools::promote_args;

  /**
   * standard one compartment PK ODE functor
   */
  struct PKOneCptODE {
  /**
   * standard one compartment PK ODE RHS function
   * @tparam T0 t type
   * @tparam T1 initial condition type
   * @tparam T2 parameter type
   * @tparam T3 real data/rate type
   * @param t type
   * @param x initial condition type
   * @param parms parameters
   * @param rate dosing rate
   * @param dummy dummy
   */
    template <typename T0, typename T1, typename T2, typename T3>
    inline
    std::vector<typename stan::return_type<T0, T1, T2, T3>::type>
    operator()(const T0& t,
               const std::vector<T1>& x,
               const std::vector<T2>& parms,
               const std::vector<T3>& rate,
               const std::vector<int>& dummy,
               std::ostream* pstream__) const {
      typedef typename stan::return_type<T0, T1, T2, T3>::type scalar;

      scalar CL = parms[0], V1 = parms[1], ka = parms[2], k10 = CL / V1;
      std::vector<scalar> y(2, 0);

      y[0] = -ka * x[0];
      y[1] = ka * x[0] - k10 * x[1];

      return y;
    }
  };

  /**
   * One-compartment PK model. The static memebers provide
   * universal information, i.e. nb. of compartments,
   * nb. of parameters, and the RHS functor.
   *
   * @tparam T_time t type
   * @tparam T_init initial condition type
   * @tparam T_rate dosing rate type
   * @tparam T_par PK parameters type
   */
  template<typename T_time, typename T_init, typename T_rate, typename T_par>
  class PKOneCptModel {
    const T_time &t0_;
    const Eigen::Matrix<T_init, 1, Eigen::Dynamic>& y0_;
    const std::vector<T_rate> &rate_;
    const T_par &CL_;
    const T_par &V2_;
    const T_par &ka_;
    const T_par k10_;
    const std::vector<T_par> alpha_;

  public:
    static constexpr int Ncmt = 2;
    static constexpr int Npar = 3;
    static constexpr PKOneCptODE f_ = PKOneCptODE();

    using scalar_type = typename promote_args<T_time, T_rate, T_par, T_init>::type;
    using aug_par_type = typename promote_args<T_rate, T_par, T_init>::type;
    using init_type   = T_init;
    using time_type   = T_time;
    using par_type    = T_par;
    using rate_type   = T_rate;
    using f_type      = PKOneCptODE;

  /**
   * One-compartment PK model constructor
   *
   * @param t0 initial time
   * @param y0 initial condition
   * @param rate dosing rate
   * @param par model parameters
   * @param CL clearance
   * @param V2 central cpt vol
   * @param ka absorption
   */
    PKOneCptModel(const T_time& t0,
                  const Eigen::Matrix<T_init, 1, Eigen::Dynamic>& y0,
                  const std::vector<T_rate> &rate,
                  const T_par& CL,
                  const T_par& V2,
                  const T_par& ka) :
      t0_(t0),
      y0_(y0),
      rate_(rate),
      CL_(CL),
      V2_(V2),
      ka_(ka),
      k10_(CL / V2),
      alpha_{k10_, ka_}
    {}

  /**
   * One-compartment PK model constructor
   * FIXME need to remove parameter as this is for linode only.
   *
   * @tparam T_mp parameters class
   * @tparam Ts parameter types
   * @param t0 initial time
   * @param y0 initial condition
   * @param rate dosing rate
   * @param par model parameters
   * @param parameter ModelParameter type
   */
    template<template<typename...> class T_mp, typename... Ts>
    PKOneCptModel(const T_time& t0,
                  const Eigen::Matrix<T_init, 1, Eigen::Dynamic>& y0,
                  const std::vector<T_rate> &rate,
                  const std::vector<T_par> & par,
                  const T_mp<Ts...> &parameter) :
      PKOneCptModel(t0, y0, rate, par[0], par[1], par[2])
    {}

  /**
   * One-compartment PK model constructor
   *
   * @param t0 initial time
   * @param y0 initial condition
   * @param rate dosing rate
   * @param par model parameters
   */
    PKOneCptModel(const T_time& t0,
                  const Eigen::Matrix<T_init, 1, Eigen::Dynamic>& y0,
                  const std::vector<T_rate> &rate,
                  const std::vector<T_par> & par) :
      PKOneCptModel(t0, y0, rate, par[0], par[1], par[2])
    {}

  /**
   * One-compartment PK model get methods
   */
    const T_time              & t0()      const { return t0_;    }
    const PKRec<T_init>       & y0()      const { return y0_;    }
    const std::vector<T_rate> & rate()    const { return rate_;  }
    const T_par               & CL()      const { return CL_;    }
    const T_par               & V2()      const { return V2_;    }
    const T_par               & ka()      const { return ka_;    }
    const T_par               & k10()     const { return k10_;   }
    const std::vector<T_par>  & alpha()   const { return alpha_; }
    const std::vector<T_par>    par()     const { return {CL_, V2_, ka_}; }
    const PKOneCptODE         & f()       const { return f_;     }
    const int                 & ncmt ()   const { return Ncmt;   }

  };

  template<typename T_time, typename T_init, typename T_rate, typename T_par>
  constexpr int PKOneCptModel<T_time, T_init, T_rate, T_par>::Ncmt;

  template<typename T_time, typename T_init, typename T_rate, typename T_par>
  constexpr int PKOneCptModel<T_time, T_init, T_rate, T_par>::Npar;

  template<typename T_time, typename T_init, typename T_rate, typename T_par>
  constexpr PKOneCptODE PKOneCptModel<T_time, T_init, T_rate, T_par>::f_;


  /**
   * Barebone definition of a one-cpt model requires only a
   * vector of parameters. Sometimes this is enough.
   *
   * @tparam T_par parameter type
   */
  template<typename T_par>
  class PKOneCptModelParameters {
    const T_par &CL_;
    const T_par &V2_;
    const T_par &ka_;
    const T_par k10_;
    const std::vector<T_par> alpha_;    
  public:
    static constexpr int Ncmt = 2;
    static constexpr int Npar = 3;
    using scalar_type = T_par;
    using par_type    = T_par;
  /**
   * Constructor
   * @tparam par parameter vector of size three.
   */
    PKOneCptModelParameters(const std::vector<T_par> & par) :
      CL_(par[0]),
      V2_(par[1]),
      ka_(par[2]),
      k10_(CL_ / V2_),
      alpha_{k10_, ka_}
    {}
  /**
   * Get methods
   */
    const T_par               & CL()      const { return CL_;    }
    const T_par               & V2()      const { return V2_;    }
    const T_par               & ka()      const { return ka_;    }
    const T_par               & k10()     const { return k10_;   }
    const std::vector<T_par>  & alpha()   const { return alpha_; }
    const int                 & ncmt ()   const { return Ncmt;   }
  }; 
  template<typename T_par>
  constexpr int PKOneCptModelParameters<T_par>::Ncmt;
  template<typename T_par>
  constexpr int PKOneCptModelParameters<T_par>::Npar;

}

#endif
