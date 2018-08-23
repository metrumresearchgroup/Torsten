#ifndef STAN_MATH_TORSTEN_TWOCPT_MODEL_HPP
#define STAN_MATH_TORSTEN_TWOCPT_MODEL_HPP

#include <stan/math/torsten/torsten_def.hpp>
#include <stan/math/prim/scal/err/check_positive.hpp>
#include <stan/math/prim/scal/err/check_finite.hpp>

namespace refactor {

  using boost::math::tools::promote_args;

  /**
   * standard two compartment PK ODE functor
   */
  struct PKTwoCptODE {
  /**
   * standard two compartment PK ODE RHS function
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
               const std::vector<int>& dummy, std::ostream* pstream__) const {
      typedef typename stan::return_type<T0, T1, T2, T3>::type scalar;

      scalar
        CL = parms.at(0),
        Q = parms.at(1),
        V1 = parms.at(2),
        V2 = parms.at(3),
        ka = parms.at(4),
        k10 = CL / V1,
        k12 = Q / V1,
        k21 = Q / V2;

      std::vector<scalar> y(3, 0);
      y.at(0) = -ka * x.at(0);
      y.at(1) = ka * x.at(0) - (k10 + k12) * x.at(1) + k21 * x.at(2);
      y.at(2) = k12 * x.at(1) - k21 * x.at(2);

      return y;
    }
  };

  /**
   * two-compartment PK model. The static memebers provide
   * universal information, i.e. nb. of compartments,
   * nb. of parameters, and the RHS functor.
   *
   * @tparam T_time t type
   * @tparam T_init initial condition type
   * @tparam T_rate dosing rate type
   * @tparam T_par PK parameters type
   */
  template<typename T_time, typename T_init, typename T_rate, typename T_par>
  class PKTwoCptModel {
    const T_time &t0_;
    const Eigen::Matrix<T_init, 1, Eigen::Dynamic>& y0_;
    const std::vector<T_rate> &rate_;
    const T_par &CL_;
    const T_par &Q_;
    const T_par &V2_;
    const T_par &V3_;
    const T_par &ka_;
    const T_par k10_;
    const T_par k12_;
    const T_par k21_;
    const T_par ksum_;
    const std::vector<T_par> alpha_;
    const std::vector<T_par> par_;

  public:
    static constexpr int Ncmt = 3;
    static constexpr int Npar = 5;
    static constexpr PKTwoCptODE f_ = PKTwoCptODE();

  /**
   * Two-compartment PK model constructor
   *
   * @param t0 initial time
   * @param y0 initial condition
   * @param rate dosing rate
   * @param par model parameters
   * @param CL clearance
   * @param Q distributed amt
   * @param V2 central cpt vol
   * @param V3 peri cpt vol
   * @param ka absorption
   */
    PKTwoCptModel(const T_time& t0,
                  const Eigen::Matrix<T_init, 1, Eigen::Dynamic>& y0,
                  const std::vector<T_rate> &rate,
                  const T_par& CL,
                  const T_par& Q,
                  const T_par& V2,
                  const T_par& V3,
                  const T_par& ka) :
      t0_(t0),
      y0_(y0),
      rate_(rate),
      CL_(CL),
      Q_(Q),
      V2_(V2),
      V3_(V3),
      ka_(ka),
      k10_(CL_ / V2_),
      k12_(Q_ / V2_),
      k21_(Q_ / V3_),
      ksum_(k10_ + k12_ + k21_),
      alpha_{0.5 * (ksum_ + sqrt(ksum_ * ksum_ - 4 * k10_ * k21_)),
        0.5 * (ksum_ - sqrt(ksum_ * ksum_ - 4 * k10_ * k21_)),
        ka_},
      par_{CL_, Q_, V2_, V3_, ka_}
    {
      using stan::math::check_positive;
      const char* fun = "PKTwoCptModel";
      check_positive(fun, "CL", CL_);
      check_positive(fun, "Q", Q_);
      check_positive(fun, "V2", V2_);
      check_positive(fun, "V3", V3_);
      check_positive(fun, "ka", ka_);
    }

  /**
   * two-Compartment PK model constructor
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
    PKTwoCptModel(const T_time& t0,
                  const Eigen::Matrix<T_init, 1, Eigen::Dynamic>& y0,
                  const std::vector<T_rate> &rate,
                  const std::vector<T_par> & par,
                  const T_mp<Ts...> &parameter) :
      PKTwoCptModel(t0, y0, rate, par[0], par[1], par[2], par[3], par[4])
    {}

  /**
   * two-compartment PK model constructor
   *
   * @param t0 initial time
   * @param y0 initial condition
   * @param rate dosing rate
   * @param par model parameters
   */
    PKTwoCptModel(const T_time& t0,
                  const Eigen::Matrix<T_init, 1, Eigen::Dynamic>& y0,
                  const std::vector<T_rate> &rate,
                  const std::vector<T_par> & par) :
      PKTwoCptModel(t0, y0, rate, par[0], par[1], par[2], par[3], par[4])
    {}

  /**
   * two-compartment PK model get methods
   */
    const T_time              & t0()      const { return t0_;    }
    const PKRec<T_init>    & y0()         const { return y0_;    }
    const std::vector<T_rate> & rate()    const { return rate_;  }
    const T_par               & CL()      const { return CL_;    }
    const T_par               & Q()       const { return Q_;     }
    const T_par               & V2()      const { return V2_;    }
    const T_par               & V3()      const { return V3_;    }
    const T_par               & ka()      const { return ka_;    }
    const T_par               & k10()     const { return k10_;   }
    const T_par               & k12()     const { return k12_;   }
    const T_par               & k21()     const { return k21_;   }
    const std::vector<T_par>  & par()     const { return par_;   }
    const std::vector<T_par>  & alpha()   const { return alpha_; }
    const PKTwoCptODE         & f()       const { return f_;     }
    const int                 & ncmt ()   const { return Ncmt;   }

  };

  template<typename T_time, typename T_init, typename T_rate, typename T_par>
  constexpr int PKTwoCptModel<T_time, T_init, T_rate, T_par>::Ncmt;

  template<typename T_time, typename T_init, typename T_rate, typename T_par>
  constexpr int PKTwoCptModel<T_time, T_init, T_rate, T_par>::Npar;

  template<typename T_time, typename T_init, typename T_rate, typename T_par>
  constexpr PKTwoCptODE PKTwoCptModel<T_time, T_init, T_rate, T_par>::f_;

}

#endif
