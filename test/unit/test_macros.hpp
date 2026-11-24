#ifndef STAN_MATH_TORSTEN_TEST_MACROS
#define STAN_MATH_TORSTEN_TEST_MACROS

#include <stan/math/torsten/test/unit/test_util.hpp>

/*
 * Macro to test overloaded torsten functions when @c theta,
 * @c biovar and @c tlag that can be constatnt or time-dependent.
 */
#define TORSTEN_CPT_PARAM_OVERLOAD_TEST(FUN, TIME, AMT, RATE, II, EVID, CMT, ADDL, SS, THETA, BIOVAR, TLAG, \
                                    EPS_VAL, EPS_GRAD)                                                      \
  {                                                                                                         \
  std::vector<std::vector<stan::math::var> > theta_v(1, stan::math::to_var(THETA[0]));                      \
  std::vector<std::vector<stan::math::var> > biovar_v(1, stan::math::to_var(BIOVAR[0]));                    \
  std::vector<std::vector<stan::math::var> > tlag_v(1, stan::math::to_var(tlag[0]));                        \
  {                                                                                                         \
      auto x0 = FUN(TIME, AMT, RATE, II, EVID, CMT, ADDL, SS, theta_v,    BIOVAR,    TLAG);                 \
      auto x1 = FUN(TIME, AMT, RATE, II, EVID, CMT, ADDL, SS, theta_v[0], BIOVAR,    TLAG);                 \
      auto x2 = FUN(TIME, AMT, RATE, II, EVID, CMT, ADDL, SS, theta_v[0], BIOVAR[0], TLAG);                 \
      auto x3 = FUN(TIME, AMT, RATE, II, EVID, CMT, ADDL, SS, theta_v[0], BIOVAR[0], TLAG[0]);              \
      auto x4 = FUN(TIME, AMT, RATE, II, EVID, CMT, ADDL, SS, theta_v[0], BIOVAR,    TLAG[0]);              \
      auto x5 = FUN(TIME, AMT, RATE, II, EVID, CMT, ADDL, SS, theta_v,    BIOVAR[0], TLAG);                 \
      auto x6 = FUN(TIME, AMT, RATE, II, EVID, CMT, ADDL, SS, theta_v,    BIOVAR[0], TLAG[0]);              \
      auto x7 = FUN(TIME, AMT, RATE, II, EVID, CMT, ADDL, SS, theta_v,    BIOVAR,    TLAG[0]);              \
      torsten::test::test_grad(theta_v[0],  x0, x1, EPS_VAL, EPS_GRAD);                                     \
      torsten::test::test_grad(theta_v[0],  x0, x2, EPS_VAL, EPS_GRAD);                                     \
      torsten::test::test_grad(theta_v[0],  x0, x3, EPS_VAL, EPS_GRAD);                                     \
      torsten::test::test_grad(theta_v[0],  x0, x4, EPS_VAL, EPS_GRAD);                                     \
      torsten::test::test_grad(theta_v[0],  x0, x5, EPS_VAL, EPS_GRAD);                                     \
      torsten::test::test_grad(theta_v[0],  x0, x6, EPS_VAL, EPS_GRAD);                                     \
      torsten::test::test_grad(theta_v[0],  x0, x7, EPS_VAL, EPS_GRAD);                                     \
  }                                                                                                         \
  {                                                                                                         \
      auto x0 = FUN(TIME, AMT, RATE, II, EVID, CMT, ADDL, SS, THETA,    biovar_v,    TLAG);                 \
      auto x1 = FUN(TIME, AMT, RATE, II, EVID, CMT, ADDL, SS, THETA[0], biovar_v,    TLAG);                 \
      auto x2 = FUN(TIME, AMT, RATE, II, EVID, CMT, ADDL, SS, THETA[0], biovar_v[0], TLAG);                 \
      auto x3 = FUN(TIME, AMT, RATE, II, EVID, CMT, ADDL, SS, THETA[0], biovar_v[0], TLAG[0]);              \
      auto x4 = FUN(TIME, AMT, RATE, II, EVID, CMT, ADDL, SS, THETA[0], biovar_v,    TLAG[0]);              \
      auto x5 = FUN(TIME, AMT, RATE, II, EVID, CMT, ADDL, SS, THETA,    biovar_v[0], TLAG);                 \
      auto x6 = FUN(TIME, AMT, RATE, II, EVID, CMT, ADDL, SS, THETA,    biovar_v[0], TLAG[0]);              \
      auto x7 = FUN(TIME, AMT, RATE, II, EVID, CMT, ADDL, SS, THETA,    biovar_v,    TLAG[0]);              \
      torsten::test::test_grad(biovar_v[0],  x0, x1, EPS_VAL, EPS_GRAD);                                    \
      torsten::test::test_grad(biovar_v[0],  x0, x2, EPS_VAL, EPS_GRAD);                                    \
      torsten::test::test_grad(biovar_v[0],  x0, x3, EPS_VAL, EPS_GRAD);                                    \
      torsten::test::test_grad(biovar_v[0],  x0, x4, EPS_VAL, EPS_GRAD);                                    \
      torsten::test::test_grad(biovar_v[0],  x0, x5, EPS_VAL, EPS_GRAD);                                    \
      torsten::test::test_grad(biovar_v[0],  x0, x6, EPS_VAL, EPS_GRAD);                                    \
      torsten::test::test_grad(biovar_v[0],  x0, x7, EPS_VAL, EPS_GRAD);                                    \
  }                                                                                                         \
  {                                                                                                         \
      auto x0 = FUN(TIME, AMT, RATE, II, EVID, CMT, ADDL, SS, THETA,    BIOVAR,    tlag_v);                 \
      auto x1 = FUN(TIME, AMT, RATE, II, EVID, CMT, ADDL, SS, THETA[0], BIOVAR,    tlag_v);                 \
      auto x2 = FUN(TIME, AMT, RATE, II, EVID, CMT, ADDL, SS, THETA[0], BIOVAR[0], tlag_v);                 \
      auto x3 = FUN(TIME, AMT, RATE, II, EVID, CMT, ADDL, SS, THETA[0], BIOVAR[0], tlag_v[0]);              \
      auto x4 = FUN(TIME, AMT, RATE, II, EVID, CMT, ADDL, SS, THETA[0], BIOVAR,    tlag_v[0]);              \
      auto x5 = FUN(TIME, AMT, RATE, II, EVID, CMT, ADDL, SS, THETA,    BIOVAR[0], tlag_v);                 \
      auto x6 = FUN(TIME, AMT, RATE, II, EVID, CMT, ADDL, SS, THETA,    BIOVAR[0], tlag_v[0]);              \
      auto x7 = FUN(TIME, AMT, RATE, II, EVID, CMT, ADDL, SS, THETA,    BIOVAR,    tlag_v[0]);              \
      torsten::test::test_grad(tlag_v[0],  x0, x1, EPS_VAL, EPS_GRAD);                                      \
      torsten::test::test_grad(tlag_v[0],  x0, x2, EPS_VAL, EPS_GRAD);                                      \
      torsten::test::test_grad(tlag_v[0],  x0, x3, EPS_VAL, EPS_GRAD);                                      \
      torsten::test::test_grad(tlag_v[0],  x0, x4, EPS_VAL, EPS_GRAD);                                      \
      torsten::test::test_grad(tlag_v[0],  x0, x5, EPS_VAL, EPS_GRAD);                                      \
      torsten::test::test_grad(tlag_v[0],  x0, x6, EPS_VAL, EPS_GRAD);                                      \
      torsten::test::test_grad(tlag_v[0],  x0, x7, EPS_VAL, EPS_GRAD);                                      \
  }                                                                                                         \
  {                                                                                                         \
      auto x0 = FUN(TIME, AMT, RATE, II, EVID, CMT, ADDL, SS, theta_v,    biovar_v,    TLAG);               \
      auto x1 = FUN(TIME, AMT, RATE, II, EVID, CMT, ADDL, SS, theta_v[0], biovar_v,    TLAG);               \
      auto x2 = FUN(TIME, AMT, RATE, II, EVID, CMT, ADDL, SS, theta_v[0], biovar_v[0], TLAG);               \
      auto x3 = FUN(TIME, AMT, RATE, II, EVID, CMT, ADDL, SS, theta_v[0], biovar_v[0], TLAG[0]);            \
      auto x4 = FUN(TIME, AMT, RATE, II, EVID, CMT, ADDL, SS, theta_v[0], biovar_v,    TLAG[0]);            \
      auto x5 = FUN(TIME, AMT, RATE, II, EVID, CMT, ADDL, SS, theta_v,    biovar_v[0], TLAG);               \
      auto x6 = FUN(TIME, AMT, RATE, II, EVID, CMT, ADDL, SS, theta_v,    biovar_v[0], TLAG[0]);            \
      auto x7 = FUN(TIME, AMT, RATE, II, EVID, CMT, ADDL, SS, theta_v,    biovar_v,    TLAG[0]);            \
      torsten::test::test_grad(theta_v[0],  x0, x1, EPS_VAL, EPS_GRAD);                                     \
      torsten::test::test_grad(theta_v[0],  x0, x2, EPS_VAL, EPS_GRAD);                                     \
      torsten::test::test_grad(theta_v[0],  x0, x3, EPS_VAL, EPS_GRAD);                                     \
      torsten::test::test_grad(theta_v[0],  x0, x4, EPS_VAL, EPS_GRAD);                                     \
      torsten::test::test_grad(theta_v[0],  x0, x5, EPS_VAL, EPS_GRAD);                                     \
      torsten::test::test_grad(theta_v[0],  x0, x6, EPS_VAL, EPS_GRAD);                                     \
      torsten::test::test_grad(theta_v[0],  x0, x7, EPS_VAL, EPS_GRAD);                                     \
      torsten::test::test_grad(biovar_v[0],  x0, x1, EPS_VAL, EPS_GRAD);                                    \
      torsten::test::test_grad(biovar_v[0],  x0, x2, EPS_VAL, EPS_GRAD);                                    \
      torsten::test::test_grad(biovar_v[0],  x0, x3, EPS_VAL, EPS_GRAD);                                    \
      torsten::test::test_grad(biovar_v[0],  x0, x4, EPS_VAL, EPS_GRAD);                                    \
      torsten::test::test_grad(biovar_v[0],  x0, x5, EPS_VAL, EPS_GRAD);                                    \
      torsten::test::test_grad(biovar_v[0],  x0, x6, EPS_VAL, EPS_GRAD);                                    \
      torsten::test::test_grad(biovar_v[0],  x0, x7, EPS_VAL, EPS_GRAD);                                    \
  }                                                                                                         \
  {                                                                                                         \
      auto x0 = FUN(TIME, AMT, RATE, II, EVID, CMT, ADDL, SS, theta_v,    BIOVAR,    tlag_v);               \
      auto x1 = FUN(TIME, AMT, RATE, II, EVID, CMT, ADDL, SS, theta_v[0], BIOVAR,    tlag_v);               \
      auto x2 = FUN(TIME, AMT, RATE, II, EVID, CMT, ADDL, SS, theta_v[0], BIOVAR[0], tlag_v);               \
      auto x3 = FUN(TIME, AMT, RATE, II, EVID, CMT, ADDL, SS, theta_v[0], BIOVAR[0], tlag_v[0]);            \
      auto x4 = FUN(TIME, AMT, RATE, II, EVID, CMT, ADDL, SS, theta_v[0], BIOVAR,    tlag_v[0]);            \
      auto x5 = FUN(TIME, AMT, RATE, II, EVID, CMT, ADDL, SS, theta_v,    BIOVAR[0], tlag_v);               \
      auto x6 = FUN(TIME, AMT, RATE, II, EVID, CMT, ADDL, SS, theta_v,    BIOVAR[0], tlag_v[0]);            \
      auto x7 = FUN(TIME, AMT, RATE, II, EVID, CMT, ADDL, SS, theta_v,    BIOVAR,    tlag_v[0]);            \
      torsten::test::test_grad(theta_v[0],  x0, x1, EPS_VAL, EPS_GRAD);                                     \
      torsten::test::test_grad(theta_v[0],  x0, x2, EPS_VAL, EPS_GRAD);                                     \
      torsten::test::test_grad(theta_v[0],  x0, x3, EPS_VAL, EPS_GRAD);                                     \
      torsten::test::test_grad(theta_v[0],  x0, x4, EPS_VAL, EPS_GRAD);                                     \
      torsten::test::test_grad(theta_v[0],  x0, x5, EPS_VAL, EPS_GRAD);                                     \
      torsten::test::test_grad(theta_v[0],  x0, x6, EPS_VAL, EPS_GRAD);                                     \
      torsten::test::test_grad(theta_v[0],  x0, x7, EPS_VAL, EPS_GRAD);                                     \
      torsten::test::test_grad(tlag_v[0],  x0, x1, EPS_VAL, EPS_GRAD);                                      \
      torsten::test::test_grad(tlag_v[0],  x0, x2, EPS_VAL, EPS_GRAD);                                      \
      torsten::test::test_grad(tlag_v[0],  x0, x3, EPS_VAL, EPS_GRAD);                                      \
      torsten::test::test_grad(tlag_v[0],  x0, x4, EPS_VAL, EPS_GRAD);                                      \
      torsten::test::test_grad(tlag_v[0],  x0, x5, EPS_VAL, EPS_GRAD);                                      \
      torsten::test::test_grad(tlag_v[0],  x0, x6, EPS_VAL, EPS_GRAD);                                      \
      torsten::test::test_grad(tlag_v[0],  x0, x7, EPS_VAL, EPS_GRAD);                                      \
  }                                                                                                         \
  {                                                                                                         \
      auto x0 = FUN(TIME, AMT, RATE, II, EVID, CMT, ADDL, SS, THETA,    biovar_v,    tlag_v);               \
      auto x1 = FUN(TIME, AMT, RATE, II, EVID, CMT, ADDL, SS, THETA[0], biovar_v,    tlag_v);               \
      auto x2 = FUN(TIME, AMT, RATE, II, EVID, CMT, ADDL, SS, THETA[0], biovar_v[0], tlag_v);               \
      auto x3 = FUN(TIME, AMT, RATE, II, EVID, CMT, ADDL, SS, THETA[0], biovar_v[0], tlag_v[0]);            \
      auto x4 = FUN(TIME, AMT, RATE, II, EVID, CMT, ADDL, SS, THETA[0], biovar_v,    tlag_v[0]);            \
      auto x5 = FUN(TIME, AMT, RATE, II, EVID, CMT, ADDL, SS, THETA,    biovar_v[0], tlag_v);               \
      auto x6 = FUN(TIME, AMT, RATE, II, EVID, CMT, ADDL, SS, THETA,    biovar_v[0], tlag_v[0]);            \
      auto x7 = FUN(TIME, AMT, RATE, II, EVID, CMT, ADDL, SS, THETA,    biovar_v,    tlag_v[0]);            \
      torsten::test::test_grad(biovar_v[0],  x0, x1, EPS_VAL, EPS_GRAD);                                    \
      torsten::test::test_grad(biovar_v[0],  x0, x2, EPS_VAL, EPS_GRAD);                                    \
      torsten::test::test_grad(biovar_v[0],  x0, x3, EPS_VAL, EPS_GRAD);                                    \
      torsten::test::test_grad(biovar_v[0],  x0, x4, EPS_VAL, EPS_GRAD);                                    \
      torsten::test::test_grad(biovar_v[0],  x0, x5, EPS_VAL, EPS_GRAD);                                    \
      torsten::test::test_grad(biovar_v[0],  x0, x6, EPS_VAL, EPS_GRAD);                                    \
      torsten::test::test_grad(biovar_v[0],  x0, x7, EPS_VAL, EPS_GRAD);                                    \
      torsten::test::test_grad(tlag_v[0],  x0, x1, EPS_VAL, EPS_GRAD);                                      \
      torsten::test::test_grad(tlag_v[0],  x0, x2, EPS_VAL, EPS_GRAD);                                      \
      torsten::test::test_grad(tlag_v[0],  x0, x3, EPS_VAL, EPS_GRAD);                                      \
      torsten::test::test_grad(tlag_v[0],  x0, x4, EPS_VAL, EPS_GRAD);                                      \
      torsten::test::test_grad(tlag_v[0],  x0, x5, EPS_VAL, EPS_GRAD);                                      \
      torsten::test::test_grad(tlag_v[0],  x0, x6, EPS_VAL, EPS_GRAD);                                      \
      torsten::test::test_grad(tlag_v[0],  x0, x7, EPS_VAL, EPS_GRAD);                                      \
  }                                                                                                         \
  {                                                                                                         \
      auto x0 = FUN(TIME, AMT, RATE, II, EVID, CMT, ADDL, SS, theta_v,    biovar_v,    tlag_v);             \
      auto x1 = FUN(TIME, AMT, RATE, II, EVID, CMT, ADDL, SS, theta_v[0], biovar_v,    tlag_v);             \
      auto x2 = FUN(TIME, AMT, RATE, II, EVID, CMT, ADDL, SS, theta_v[0], biovar_v[0], tlag_v);             \
      auto x3 = FUN(TIME, AMT, RATE, II, EVID, CMT, ADDL, SS, theta_v[0], biovar_v[0], tlag_v[0]);          \
      auto x4 = FUN(TIME, AMT, RATE, II, EVID, CMT, ADDL, SS, theta_v[0], biovar_v,    tlag_v[0]);          \
      auto x5 = FUN(TIME, AMT, RATE, II, EVID, CMT, ADDL, SS, theta_v,    biovar_v[0], tlag_v);             \
      auto x6 = FUN(TIME, AMT, RATE, II, EVID, CMT, ADDL, SS, theta_v,    biovar_v[0], tlag_v[0]);          \
      auto x7 = FUN(TIME, AMT, RATE, II, EVID, CMT, ADDL, SS, theta_v,    biovar_v,    tlag_v[0]);          \
      torsten::test::test_grad(theta_v[0],  x0, x1, EPS_VAL, EPS_GRAD);                                     \
      torsten::test::test_grad(theta_v[0],  x0, x2, EPS_VAL, EPS_GRAD);                                     \
      torsten::test::test_grad(theta_v[0],  x0, x3, EPS_VAL, EPS_GRAD);                                     \
      torsten::test::test_grad(theta_v[0],  x0, x4, EPS_VAL, EPS_GRAD);                                     \
      torsten::test::test_grad(theta_v[0],  x0, x5, EPS_VAL, EPS_GRAD);                                     \
      torsten::test::test_grad(theta_v[0],  x0, x6, EPS_VAL, EPS_GRAD);                                     \
      torsten::test::test_grad(theta_v[0],  x0, x7, EPS_VAL, EPS_GRAD);                                     \
      torsten::test::test_grad(biovar_v[0],  x0, x1, EPS_VAL, EPS_GRAD);                                    \
      torsten::test::test_grad(biovar_v[0],  x0, x2, EPS_VAL, EPS_GRAD);                                    \
      torsten::test::test_grad(biovar_v[0],  x0, x3, EPS_VAL, EPS_GRAD);                                    \
      torsten::test::test_grad(biovar_v[0],  x0, x4, EPS_VAL, EPS_GRAD);                                    \
      torsten::test::test_grad(biovar_v[0],  x0, x5, EPS_VAL, EPS_GRAD);                                    \
      torsten::test::test_grad(biovar_v[0],  x0, x6, EPS_VAL, EPS_GRAD);                                    \
      torsten::test::test_grad(biovar_v[0],  x0, x7, EPS_VAL, EPS_GRAD);                                    \
      torsten::test::test_grad(tlag_v[0],  x0, x1, EPS_VAL, EPS_GRAD);                                      \
      torsten::test::test_grad(tlag_v[0],  x0, x2, EPS_VAL, EPS_GRAD);                                      \
      torsten::test::test_grad(tlag_v[0],  x0, x3, EPS_VAL, EPS_GRAD);                                      \
      torsten::test::test_grad(tlag_v[0],  x0, x4, EPS_VAL, EPS_GRAD);                                      \
      torsten::test::test_grad(tlag_v[0],  x0, x5, EPS_VAL, EPS_GRAD);                                      \
      torsten::test::test_grad(tlag_v[0],  x0, x6, EPS_VAL, EPS_GRAD);                                      \
      torsten::test::test_grad(tlag_v[0],  x0, x7, EPS_VAL, EPS_GRAD);                                      \
  }                                                                                                         \
  }

/*
 * Macro to test overloaded torsten functions when @c theta,
 * @c biovar and @c tlag that can be constatnt or time-dependent.
 */
#define TORSTEN_ODE_PARAM_OVERLOAD_TEST(FUN, F, NCMT, TIME, AMT, RATE, II, EVID, CMT, ADDL, SS, THETA, BIOVAR, TLAG, \
                                    EPS_VAL, EPS_GRAD)                                                               \
  {                                                                                                                  \
  std::vector<std::vector<stan::math::var> > theta_v(1, stan::math::to_var(THETA[0]));                               \
  std::vector<std::vector<stan::math::var> > biovar_v(1, stan::math::to_var(BIOVAR[0]));                             \
  std::vector<std::vector<stan::math::var> > tlag_v(1, stan::math::to_var(tlag[0]));                                 \
  {                                                                                                                  \
      auto x0 = FUN(F, NCMT, TIME, AMT, RATE, II, EVID, CMT, ADDL, SS, theta_v,    BIOVAR,    TLAG, nullptr);                 \
      auto x1 = FUN(F, NCMT, TIME, AMT, RATE, II, EVID, CMT, ADDL, SS, theta_v[0], BIOVAR,    TLAG, nullptr);                 \
      auto x2 = FUN(F, NCMT, TIME, AMT, RATE, II, EVID, CMT, ADDL, SS, theta_v[0], BIOVAR[0], TLAG, nullptr);                 \
      auto x3 = FUN(F, NCMT, TIME, AMT, RATE, II, EVID, CMT, ADDL, SS, theta_v[0], BIOVAR[0], TLAG[0], nullptr);              \
      auto x4 = FUN(F, NCMT, TIME, AMT, RATE, II, EVID, CMT, ADDL, SS, theta_v[0], BIOVAR,    TLAG[0], nullptr);              \
      auto x5 = FUN(F, NCMT, TIME, AMT, RATE, II, EVID, CMT, ADDL, SS, theta_v,    BIOVAR[0], TLAG, nullptr);                 \
      auto x6 = FUN(F, NCMT, TIME, AMT, RATE, II, EVID, CMT, ADDL, SS, theta_v,    BIOVAR[0], TLAG[0], nullptr);              \
      auto x7 = FUN(F, NCMT, TIME, AMT, RATE, II, EVID, CMT, ADDL, SS, theta_v,    BIOVAR,    TLAG[0], nullptr);              \
      torsten::test::test_grad(theta_v[0],  x0, x1, EPS_VAL, EPS_GRAD);                                              \
      torsten::test::test_grad(theta_v[0],  x0, x2, EPS_VAL, EPS_GRAD);                                              \
      torsten::test::test_grad(theta_v[0],  x0, x3, EPS_VAL, EPS_GRAD);                                              \
      torsten::test::test_grad(theta_v[0],  x0, x4, EPS_VAL, EPS_GRAD);                                              \
      torsten::test::test_grad(theta_v[0],  x0, x5, EPS_VAL, EPS_GRAD);                                              \
      torsten::test::test_grad(theta_v[0],  x0, x6, EPS_VAL, EPS_GRAD);                                              \
      torsten::test::test_grad(theta_v[0],  x0, x7, EPS_VAL, EPS_GRAD);                                              \
  }                                                                                                                  \
  {                                                                                                                  \
      auto x0 = FUN(F, NCMT, TIME, AMT, RATE, II, EVID, CMT, ADDL, SS, THETA,    biovar_v,    TLAG, nullptr);                 \
      auto x1 = FUN(F, NCMT, TIME, AMT, RATE, II, EVID, CMT, ADDL, SS, THETA[0], biovar_v,    TLAG, nullptr);                 \
      auto x2 = FUN(F, NCMT, TIME, AMT, RATE, II, EVID, CMT, ADDL, SS, THETA[0], biovar_v[0], TLAG, nullptr);                 \
      auto x3 = FUN(F, NCMT, TIME, AMT, RATE, II, EVID, CMT, ADDL, SS, THETA[0], biovar_v[0], TLAG[0], nullptr);              \
      auto x4 = FUN(F, NCMT, TIME, AMT, RATE, II, EVID, CMT, ADDL, SS, THETA[0], biovar_v,    TLAG[0], nullptr);              \
      auto x5 = FUN(F, NCMT, TIME, AMT, RATE, II, EVID, CMT, ADDL, SS, THETA,    biovar_v[0], TLAG, nullptr);                 \
      auto x6 = FUN(F, NCMT, TIME, AMT, RATE, II, EVID, CMT, ADDL, SS, THETA,    biovar_v[0], TLAG[0], nullptr);              \
      auto x7 = FUN(F, NCMT, TIME, AMT, RATE, II, EVID, CMT, ADDL, SS, THETA,    biovar_v,    TLAG[0], nullptr);              \
      torsten::test::test_grad(biovar_v[0],  x0, x1, EPS_VAL, EPS_GRAD);                                             \
      torsten::test::test_grad(biovar_v[0],  x0, x2, EPS_VAL, EPS_GRAD);                                             \
      torsten::test::test_grad(biovar_v[0],  x0, x3, EPS_VAL, EPS_GRAD);                                             \
      torsten::test::test_grad(biovar_v[0],  x0, x4, EPS_VAL, EPS_GRAD);                                             \
      torsten::test::test_grad(biovar_v[0],  x0, x5, EPS_VAL, EPS_GRAD);                                             \
      torsten::test::test_grad(biovar_v[0],  x0, x6, EPS_VAL, EPS_GRAD);                                             \
      torsten::test::test_grad(biovar_v[0],  x0, x7, EPS_VAL, EPS_GRAD);                                             \
  }                                                                                                                  \
  {                                                                                                                  \
      auto x0 = FUN(F, NCMT, TIME, AMT, RATE, II, EVID, CMT, ADDL, SS, THETA,    BIOVAR,    tlag_v, nullptr);                 \
      auto x1 = FUN(F, NCMT, TIME, AMT, RATE, II, EVID, CMT, ADDL, SS, THETA[0], BIOVAR,    tlag_v, nullptr);                 \
      auto x2 = FUN(F, NCMT, TIME, AMT, RATE, II, EVID, CMT, ADDL, SS, THETA[0], BIOVAR[0], tlag_v, nullptr);                 \
      auto x3 = FUN(F, NCMT, TIME, AMT, RATE, II, EVID, CMT, ADDL, SS, THETA[0], BIOVAR[0], tlag_v[0], nullptr);              \
      auto x4 = FUN(F, NCMT, TIME, AMT, RATE, II, EVID, CMT, ADDL, SS, THETA[0], BIOVAR,    tlag_v[0], nullptr);              \
      auto x5 = FUN(F, NCMT, TIME, AMT, RATE, II, EVID, CMT, ADDL, SS, THETA,    BIOVAR[0], tlag_v, nullptr);                 \
      auto x6 = FUN(F, NCMT, TIME, AMT, RATE, II, EVID, CMT, ADDL, SS, THETA,    BIOVAR[0], tlag_v[0], nullptr);              \
      auto x7 = FUN(F, NCMT, TIME, AMT, RATE, II, EVID, CMT, ADDL, SS, THETA,    BIOVAR,    tlag_v[0], nullptr);              \
      torsten::test::test_grad(tlag_v[0],  x0, x1, EPS_VAL, EPS_GRAD);                                               \
      torsten::test::test_grad(tlag_v[0],  x0, x2, EPS_VAL, EPS_GRAD);                                               \
      torsten::test::test_grad(tlag_v[0],  x0, x3, EPS_VAL, EPS_GRAD);                                               \
      torsten::test::test_grad(tlag_v[0],  x0, x4, EPS_VAL, EPS_GRAD);                                               \
      torsten::test::test_grad(tlag_v[0],  x0, x5, EPS_VAL, EPS_GRAD);                                               \
      torsten::test::test_grad(tlag_v[0],  x0, x6, EPS_VAL, EPS_GRAD);                                               \
      torsten::test::test_grad(tlag_v[0],  x0, x7, EPS_VAL, EPS_GRAD);                                               \
  }                                                                                                                  \
  {                                                                                                                  \
      auto x0 = FUN(F, NCMT, TIME, AMT, RATE, II, EVID, CMT, ADDL, SS, theta_v,    biovar_v,    TLAG, nullptr);               \
      auto x1 = FUN(F, NCMT, TIME, AMT, RATE, II, EVID, CMT, ADDL, SS, theta_v[0], biovar_v,    TLAG, nullptr);               \
      auto x2 = FUN(F, NCMT, TIME, AMT, RATE, II, EVID, CMT, ADDL, SS, theta_v[0], biovar_v[0], TLAG, nullptr);               \
      auto x3 = FUN(F, NCMT, TIME, AMT, RATE, II, EVID, CMT, ADDL, SS, theta_v[0], biovar_v[0], TLAG[0], nullptr);            \
      auto x4 = FUN(F, NCMT, TIME, AMT, RATE, II, EVID, CMT, ADDL, SS, theta_v[0], biovar_v,    TLAG[0], nullptr);            \
      auto x5 = FUN(F, NCMT, TIME, AMT, RATE, II, EVID, CMT, ADDL, SS, theta_v,    biovar_v[0], TLAG, nullptr);               \
      auto x6 = FUN(F, NCMT, TIME, AMT, RATE, II, EVID, CMT, ADDL, SS, theta_v,    biovar_v[0], TLAG[0], nullptr);            \
      auto x7 = FUN(F, NCMT, TIME, AMT, RATE, II, EVID, CMT, ADDL, SS, theta_v,    biovar_v,    TLAG[0], nullptr);            \
      torsten::test::test_grad(theta_v[0],  x0, x1, EPS_VAL, EPS_GRAD);                                              \
      torsten::test::test_grad(theta_v[0],  x0, x2, EPS_VAL, EPS_GRAD);                                              \
      torsten::test::test_grad(theta_v[0],  x0, x3, EPS_VAL, EPS_GRAD);                                              \
      torsten::test::test_grad(theta_v[0],  x0, x4, EPS_VAL, EPS_GRAD);                                              \
      torsten::test::test_grad(theta_v[0],  x0, x5, EPS_VAL, EPS_GRAD);                                              \
      torsten::test::test_grad(theta_v[0],  x0, x6, EPS_VAL, EPS_GRAD);                                              \
      torsten::test::test_grad(theta_v[0],  x0, x7, EPS_VAL, EPS_GRAD);                                              \
      torsten::test::test_grad(biovar_v[0],  x0, x1, EPS_VAL, EPS_GRAD);                                             \
      torsten::test::test_grad(biovar_v[0],  x0, x2, EPS_VAL, EPS_GRAD);                                             \
      torsten::test::test_grad(biovar_v[0],  x0, x3, EPS_VAL, EPS_GRAD);                                             \
      torsten::test::test_grad(biovar_v[0],  x0, x4, EPS_VAL, EPS_GRAD);                                             \
      torsten::test::test_grad(biovar_v[0],  x0, x5, EPS_VAL, EPS_GRAD);                                             \
      torsten::test::test_grad(biovar_v[0],  x0, x6, EPS_VAL, EPS_GRAD);                                             \
      torsten::test::test_grad(biovar_v[0],  x0, x7, EPS_VAL, EPS_GRAD);                                             \
  }                                                                                                                  \
  {                                                                                                                  \
      auto x0 = FUN(F, NCMT, TIME, AMT, RATE, II, EVID, CMT, ADDL, SS, theta_v,    BIOVAR,    tlag_v, nullptr);               \
      auto x1 = FUN(F, NCMT, TIME, AMT, RATE, II, EVID, CMT, ADDL, SS, theta_v[0], BIOVAR,    tlag_v, nullptr);               \
      auto x2 = FUN(F, NCMT, TIME, AMT, RATE, II, EVID, CMT, ADDL, SS, theta_v[0], BIOVAR[0], tlag_v, nullptr);               \
      auto x3 = FUN(F, NCMT, TIME, AMT, RATE, II, EVID, CMT, ADDL, SS, theta_v[0], BIOVAR[0], tlag_v[0], nullptr);            \
      auto x4 = FUN(F, NCMT, TIME, AMT, RATE, II, EVID, CMT, ADDL, SS, theta_v[0], BIOVAR,    tlag_v[0], nullptr);            \
      auto x5 = FUN(F, NCMT, TIME, AMT, RATE, II, EVID, CMT, ADDL, SS, theta_v,    BIOVAR[0], tlag_v, nullptr);               \
      auto x6 = FUN(F, NCMT, TIME, AMT, RATE, II, EVID, CMT, ADDL, SS, theta_v,    BIOVAR[0], tlag_v[0], nullptr);            \
      auto x7 = FUN(F, NCMT, TIME, AMT, RATE, II, EVID, CMT, ADDL, SS, theta_v,    BIOVAR,    tlag_v[0], nullptr);            \
      torsten::test::test_grad(theta_v[0],  x0, x1, EPS_VAL, EPS_GRAD);                                              \
      torsten::test::test_grad(theta_v[0],  x0, x2, EPS_VAL, EPS_GRAD);                                              \
      torsten::test::test_grad(theta_v[0],  x0, x3, EPS_VAL, EPS_GRAD);                                              \
      torsten::test::test_grad(theta_v[0],  x0, x4, EPS_VAL, EPS_GRAD);                                              \
      torsten::test::test_grad(theta_v[0],  x0, x5, EPS_VAL, EPS_GRAD);                                              \
      torsten::test::test_grad(theta_v[0],  x0, x6, EPS_VAL, EPS_GRAD);                                              \
      torsten::test::test_grad(theta_v[0],  x0, x7, EPS_VAL, EPS_GRAD);                                              \
      torsten::test::test_grad(tlag_v[0],  x0, x1, EPS_VAL, EPS_GRAD);                                               \
      torsten::test::test_grad(tlag_v[0],  x0, x2, EPS_VAL, EPS_GRAD);                                               \
      torsten::test::test_grad(tlag_v[0],  x0, x3, EPS_VAL, EPS_GRAD);                                               \
      torsten::test::test_grad(tlag_v[0],  x0, x4, EPS_VAL, EPS_GRAD);                                               \
      torsten::test::test_grad(tlag_v[0],  x0, x5, EPS_VAL, EPS_GRAD);                                               \
      torsten::test::test_grad(tlag_v[0],  x0, x6, EPS_VAL, EPS_GRAD);                                               \
      torsten::test::test_grad(tlag_v[0],  x0, x7, EPS_VAL, EPS_GRAD);                                               \
  }                                                                                                                  \
  {                                                                                                                  \
      auto x0 = FUN(F, NCMT, TIME, AMT, RATE, II, EVID, CMT, ADDL, SS, THETA,    biovar_v,    tlag_v, nullptr);               \
      auto x1 = FUN(F, NCMT, TIME, AMT, RATE, II, EVID, CMT, ADDL, SS, THETA[0], biovar_v,    tlag_v, nullptr);               \
      auto x2 = FUN(F, NCMT, TIME, AMT, RATE, II, EVID, CMT, ADDL, SS, THETA[0], biovar_v[0], tlag_v, nullptr);               \
      auto x3 = FUN(F, NCMT, TIME, AMT, RATE, II, EVID, CMT, ADDL, SS, THETA[0], biovar_v[0], tlag_v[0], nullptr);            \
      auto x4 = FUN(F, NCMT, TIME, AMT, RATE, II, EVID, CMT, ADDL, SS, THETA[0], biovar_v,    tlag_v[0], nullptr);            \
      auto x5 = FUN(F, NCMT, TIME, AMT, RATE, II, EVID, CMT, ADDL, SS, THETA,    biovar_v[0], tlag_v, nullptr);               \
      auto x6 = FUN(F, NCMT, TIME, AMT, RATE, II, EVID, CMT, ADDL, SS, THETA,    biovar_v[0], tlag_v[0], nullptr);            \
      auto x7 = FUN(F, NCMT, TIME, AMT, RATE, II, EVID, CMT, ADDL, SS, THETA,    biovar_v,    tlag_v[0], nullptr);            \
      torsten::test::test_grad(biovar_v[0],  x0, x1, EPS_VAL, EPS_GRAD);                                             \
      torsten::test::test_grad(biovar_v[0],  x0, x2, EPS_VAL, EPS_GRAD);                                             \
      torsten::test::test_grad(biovar_v[0],  x0, x3, EPS_VAL, EPS_GRAD);                                             \
      torsten::test::test_grad(biovar_v[0],  x0, x4, EPS_VAL, EPS_GRAD);                                             \
      torsten::test::test_grad(biovar_v[0],  x0, x5, EPS_VAL, EPS_GRAD);                                             \
      torsten::test::test_grad(biovar_v[0],  x0, x6, EPS_VAL, EPS_GRAD);                                             \
      torsten::test::test_grad(biovar_v[0],  x0, x7, EPS_VAL, EPS_GRAD);                                             \
      torsten::test::test_grad(tlag_v[0],  x0, x1, EPS_VAL, EPS_GRAD);                                               \
      torsten::test::test_grad(tlag_v[0],  x0, x2, EPS_VAL, EPS_GRAD);                                               \
      torsten::test::test_grad(tlag_v[0],  x0, x3, EPS_VAL, EPS_GRAD);                                               \
      torsten::test::test_grad(tlag_v[0],  x0, x4, EPS_VAL, EPS_GRAD);                                               \
      torsten::test::test_grad(tlag_v[0],  x0, x5, EPS_VAL, EPS_GRAD);                                               \
      torsten::test::test_grad(tlag_v[0],  x0, x6, EPS_VAL, EPS_GRAD);                                               \
      torsten::test::test_grad(tlag_v[0],  x0, x7, EPS_VAL, EPS_GRAD);                                               \
  }                                                                                                                  \
  {                                                                                                                  \
      auto x0 = FUN(F, NCMT, TIME, AMT, RATE, II, EVID, CMT, ADDL, SS, theta_v,    biovar_v,    tlag_v, nullptr);             \
      auto x1 = FUN(F, NCMT, TIME, AMT, RATE, II, EVID, CMT, ADDL, SS, theta_v[0], biovar_v,    tlag_v, nullptr);             \
      auto x2 = FUN(F, NCMT, TIME, AMT, RATE, II, EVID, CMT, ADDL, SS, theta_v[0], biovar_v[0], tlag_v, nullptr);             \
      auto x3 = FUN(F, NCMT, TIME, AMT, RATE, II, EVID, CMT, ADDL, SS, theta_v[0], biovar_v[0], tlag_v[0], nullptr);          \
      auto x4 = FUN(F, NCMT, TIME, AMT, RATE, II, EVID, CMT, ADDL, SS, theta_v[0], biovar_v,    tlag_v[0], nullptr);          \
      auto x5 = FUN(F, NCMT, TIME, AMT, RATE, II, EVID, CMT, ADDL, SS, theta_v,    biovar_v[0], tlag_v, nullptr);             \
      auto x6 = FUN(F, NCMT, TIME, AMT, RATE, II, EVID, CMT, ADDL, SS, theta_v,    biovar_v[0], tlag_v[0], nullptr);          \
      auto x7 = FUN(F, NCMT, TIME, AMT, RATE, II, EVID, CMT, ADDL, SS, theta_v,    biovar_v,    tlag_v[0], nullptr);          \
      torsten::test::test_grad(theta_v[0],  x0, x1, EPS_VAL, EPS_GRAD);                                              \
      torsten::test::test_grad(theta_v[0],  x0, x2, EPS_VAL, EPS_GRAD);                                              \
      torsten::test::test_grad(theta_v[0],  x0, x3, EPS_VAL, EPS_GRAD);                                              \
      torsten::test::test_grad(theta_v[0],  x0, x4, EPS_VAL, EPS_GRAD);                                              \
      torsten::test::test_grad(theta_v[0],  x0, x5, EPS_VAL, EPS_GRAD);                                              \
      torsten::test::test_grad(theta_v[0],  x0, x6, EPS_VAL, EPS_GRAD);                                              \
      torsten::test::test_grad(theta_v[0],  x0, x7, EPS_VAL, EPS_GRAD);                                              \
      torsten::test::test_grad(biovar_v[0],  x0, x1, EPS_VAL, EPS_GRAD);                                             \
      torsten::test::test_grad(biovar_v[0],  x0, x2, EPS_VAL, EPS_GRAD);                                             \
      torsten::test::test_grad(biovar_v[0],  x0, x3, EPS_VAL, EPS_GRAD);                                             \
      torsten::test::test_grad(biovar_v[0],  x0, x4, EPS_VAL, EPS_GRAD);                                             \
      torsten::test::test_grad(biovar_v[0],  x0, x5, EPS_VAL, EPS_GRAD);                                             \
      torsten::test::test_grad(biovar_v[0],  x0, x6, EPS_VAL, EPS_GRAD);                                             \
      torsten::test::test_grad(biovar_v[0],  x0, x7, EPS_VAL, EPS_GRAD);                                             \
      torsten::test::test_grad(tlag_v[0],  x0, x1, EPS_VAL, EPS_GRAD);                                               \
      torsten::test::test_grad(tlag_v[0],  x0, x2, EPS_VAL, EPS_GRAD);                                               \
      torsten::test::test_grad(tlag_v[0],  x0, x3, EPS_VAL, EPS_GRAD);                                               \
      torsten::test::test_grad(tlag_v[0],  x0, x4, EPS_VAL, EPS_GRAD);                                               \
      torsten::test::test_grad(tlag_v[0],  x0, x5, EPS_VAL, EPS_GRAD);                                               \
      torsten::test::test_grad(tlag_v[0],  x0, x6, EPS_VAL, EPS_GRAD);                                               \
      torsten::test::test_grad(tlag_v[0],  x0, x7, EPS_VAL, EPS_GRAD);                                               \
  }                                                                                                                  \
  }

// test varyadic solvers
#define TORSTEN_ODE_PARAM_VARI_OVERLOAD_TEST(FUN, F, NCMT, TIME, AMT, RATE, II, EVID, CMT, ADDL, SS, THETA, BIOVAR, TLAG,                     \
                                             EPS_VAL, EPS_GRAD)                                                                               \
  {                                                                                                                                           \
  std::vector<std::vector<stan::math::var> > theta_v(1, stan::math::to_var(THETA[0]));                                                        \
  std::vector<std::vector<stan::math::var> > biovar_v(1, stan::math::to_var(BIOVAR[0]));                                                      \
  std::vector<std::vector<stan::math::var> > tlag_v(1, stan::math::to_var(tlag[0]));                                                          \
  {                                                                                                                                           \
      auto x0 = FUN(F, NCMT, TIME, AMT, RATE, II, EVID, CMT, ADDL, SS, theta_v,    BIOVAR,    TLAG, nullptr);                                          \
      auto x1 = FUN(F, NCMT, TIME, AMT, RATE, II, EVID, CMT, ADDL, SS, theta_v[0], BIOVAR,    TLAG, 1.e-6, 1.e-6, 1e6, nullptr);                       \
      auto x2 = FUN(F, NCMT, TIME, AMT, RATE, II, EVID, CMT, ADDL, SS, theta_v[0], BIOVAR[0], TLAG, 1.e-6, 1.e-6, 1e6, nullptr);                       \
      auto x3 = FUN(F, NCMT, TIME, AMT, RATE, II, EVID, CMT, ADDL, SS, theta_v[0], BIOVAR[0], TLAG[0], 1.e-6, 1.e-6, 1e6, 1.e-6, 1.e-6, 1e2, nullptr); \
      auto x4 = FUN(F, NCMT, TIME, AMT, RATE, II, EVID, CMT, ADDL, SS, theta_v[0], BIOVAR,    TLAG[0], 1.e-6, 1.e-6, 1e6, 1.e-6, 1.e-6, 1e2, nullptr); \
      auto x5 = FUN(F, NCMT, TIME, AMT, RATE, II, EVID, CMT, ADDL, SS, theta_v,    BIOVAR[0], TLAG, 1.e-6, 1.e-6, 1e6, nullptr);                       \
      auto x6 = FUN(F, NCMT, TIME, AMT, RATE, II, EVID, CMT, ADDL, SS, theta_v,    BIOVAR[0], TLAG[0], 1.e-6, 1.e-6, 1e6, nullptr);                    \
      auto x7 = FUN(F, NCMT, TIME, AMT, RATE, II, EVID, CMT, ADDL, SS, theta_v,    BIOVAR,    TLAG[0], 1.e-6, 1.e-6, 1e6, nullptr);                    \
      torsten::test::test_grad(theta_v[0],  x0, x1, EPS_VAL, EPS_GRAD);                                                                       \
      torsten::test::test_grad(theta_v[0],  x0, x2, EPS_VAL, EPS_GRAD);                                                                       \
      torsten::test::test_grad(theta_v[0],  x0, x3, EPS_VAL, EPS_GRAD);                                                                       \
      torsten::test::test_grad(theta_v[0],  x0, x4, EPS_VAL, EPS_GRAD);                                                                       \
      torsten::test::test_grad(theta_v[0],  x0, x5, EPS_VAL, EPS_GRAD);                                                                       \
      torsten::test::test_grad(theta_v[0],  x0, x6, EPS_VAL, EPS_GRAD);                                                                       \
      torsten::test::test_grad(theta_v[0],  x0, x7, EPS_VAL, EPS_GRAD);                                                                       \
  }                                                                                                                                           \
  {                                                                                                                                           \
      auto x0 = FUN(F, NCMT, TIME, AMT, RATE, II, EVID, CMT, ADDL, SS, theta_v,    BIOVAR, nullptr);                                          \
      auto x1 = FUN(F, NCMT, TIME, AMT, RATE, II, EVID, CMT, ADDL, SS, theta_v[0], 1.e-6, 1.e-6, 1e6, nullptr);                             \
      auto x2 = FUN(F, NCMT, TIME, AMT, RATE, II, EVID, CMT, ADDL, SS, theta_v[0], 1.e-6, 1.e-6, 1e6, nullptr);                             \
      auto x3 = FUN(F, NCMT, TIME, AMT, RATE, II, EVID, CMT, ADDL, SS, theta_v[0], 1.e-6, 1.e-6, 1e6, 1.e-6, 1.e-6, 1e2, nullptr);          \
      auto x4 = FUN(F, NCMT, TIME, AMT, RATE, II, EVID, CMT, ADDL, SS, theta_v[0], 1.e-6, 1.e-6, 1e6, 1.e-6, 1.e-6, 1e2, nullptr);          \
      auto x5 = FUN(F, NCMT, TIME, AMT, RATE, II, EVID, CMT, ADDL, SS, theta_v,    1.e-6, 1.e-6, 1e6, nullptr);                             \
      auto x6 = FUN(F, NCMT, TIME, AMT, RATE, II, EVID, CMT, ADDL, SS, theta_v,    1.e-6, 1.e-6, 1e6, nullptr);                             \
      auto x7 = FUN(F, NCMT, TIME, AMT, RATE, II, EVID, CMT, ADDL, SS, theta_v,    1.e-6, 1.e-6, 1e6, nullptr);                             \
      torsten::test::test_grad(theta_v[0],  x0, x1, EPS_VAL, EPS_GRAD);                                                                       \
      torsten::test::test_grad(theta_v[0],  x0, x2, EPS_VAL, EPS_GRAD);                                                                       \
      torsten::test::test_grad(theta_v[0],  x0, x3, EPS_VAL, EPS_GRAD);                                                                       \
      torsten::test::test_grad(theta_v[0],  x0, x4, EPS_VAL, EPS_GRAD);                                                                       \
      torsten::test::test_grad(theta_v[0],  x0, x5, EPS_VAL, EPS_GRAD);                                                                       \
      torsten::test::test_grad(theta_v[0],  x0, x6, EPS_VAL, EPS_GRAD);                                                                       \
      torsten::test::test_grad(theta_v[0],  x0, x7, EPS_VAL, EPS_GRAD);                                                                       \
  }                                                                                                                                           \
  }

// test varyadic solvers for tlag
#define TORSTEN_ODE_PARAM_VARI_TLAG_OVERLOAD_TEST(FUN, F, NCMT, TIME, AMT, RATE, II, EVID, CMT, ADDL, SS, THETA, BIOVAR, TLAG,                     \
                                             EPS_VAL, EPS_GRAD)                                                                               \
  {                                                                                                                                           \
  std::vector<std::vector<stan::math::var> > theta_v(1, stan::math::to_var(THETA[0]));                                                        \
  std::vector<std::vector<stan::math::var> > biovar_v(1, stan::math::to_var(BIOVAR[0]));                                                      \
  std::vector<std::vector<stan::math::var> > tlag_v(1, stan::math::to_var(tlag[0]));                                                          \
  {                                                                                                                                           \
      auto x0 = FUN(F, NCMT, TIME, AMT, RATE, II, EVID, CMT, ADDL, SS, theta_v,    BIOVAR,    TLAG, nullptr);                                          \
      auto x1 = FUN(F, NCMT, TIME, AMT, RATE, II, EVID, CMT, ADDL, SS, theta_v[0], BIOVAR,    1.e-6, 1.e-6, 1e6, nullptr);                       \
      auto x2 = FUN(F, NCMT, TIME, AMT, RATE, II, EVID, CMT, ADDL, SS, theta_v[0], BIOVAR[0], 1.e-6, 1.e-6, 1e6, nullptr);                       \
      auto x3 = FUN(F, NCMT, TIME, AMT, RATE, II, EVID, CMT, ADDL, SS, theta_v[0], BIOVAR[0], 1.e-6, 1.e-6, 1e6, 1.e-6, 1.e-6, 1e2, nullptr); \
      auto x4 = FUN(F, NCMT, TIME, AMT, RATE, II, EVID, CMT, ADDL, SS, theta_v[0], BIOVAR,    1.e-6, 1.e-6, 1e6, 1.e-6, 1.e-6, 1e2, nullptr); \
      auto x5 = FUN(F, NCMT, TIME, AMT, RATE, II, EVID, CMT, ADDL, SS, theta_v,    BIOVAR[0], 1.e-6, 1.e-6, 1e6, nullptr);                       \
      auto x6 = FUN(F, NCMT, TIME, AMT, RATE, II, EVID, CMT, ADDL, SS, theta_v,    BIOVAR[0], 1.e-6, 1.e-6, 1e6, nullptr);                    \
      auto x7 = FUN(F, NCMT, TIME, AMT, RATE, II, EVID, CMT, ADDL, SS, theta_v,    BIOVAR,    1.e-6, 1.e-6, 1e6, nullptr);                    \
      torsten::test::test_grad(theta_v[0],  x0, x1, EPS_VAL, EPS_GRAD);                                                                       \
      torsten::test::test_grad(theta_v[0],  x0, x2, EPS_VAL, EPS_GRAD);                                                                       \
      torsten::test::test_grad(theta_v[0],  x0, x3, EPS_VAL, EPS_GRAD);                                                                       \
      torsten::test::test_grad(theta_v[0],  x0, x4, EPS_VAL, EPS_GRAD);                                                                       \
      torsten::test::test_grad(theta_v[0],  x0, x5, EPS_VAL, EPS_GRAD);                                                                       \
      torsten::test::test_grad(theta_v[0],  x0, x6, EPS_VAL, EPS_GRAD);                                                                       \
      torsten::test::test_grad(theta_v[0],  x0, x7, EPS_VAL, EPS_GRAD);                                                                       \
  }                                                                                                                                           \
  }

#define TORSTEN_CPT_GRAD_THETA_TEST(FUN, TIME, AMT, RATE, II, EVID, CMT, ADDL, SS, THETA, BIOVAR, TLAG, \
                                    H, EPS_VAL, EPS_RTOL, EPS_ATOL)                                     \
  {                                                                                                     \
    auto f1 = [&] (const std::vector<std::vector<double> >& x) {                                        \
      return FUN(TIME, AMT, RATE, II, EVID, CMT, ADDL, SS, x, BIOVAR, TLAG);                            \
    };                                                                                                  \
    auto f2 = [&] (const std::vector<std::vector<stan::math::var> >& x) {                               \
      return FUN(TIME, AMT, RATE, II, EVID, CMT, ADDL, SS, x, BIOVAR, TLAG);                            \
    };                                                                                                  \
    torsten::test::test_grad(f1, f2, THETA, H, EPS_VAL, EPS_RTOL, EPS_ATOL);                            \
  }

#define TORSTEN_LIN_GRAD_THETA_TEST(FUN, TIME, AMT, RATE, II, EVID, CMT, ADDL, SS, THETA, BIOVAR, TLAG, \
                                    H, EPS_VAL, EPS_RTOL, EPS_ATOL)                                     \
  {                                                                                                     \
    auto f1 = [&] (const std::vector<Eigen::MatrixXd>& x) {                                             \
      return FUN(TIME, AMT, RATE, II, EVID, CMT, ADDL, SS, x, BIOVAR, TLAG);                            \
    };                                                                                                  \
    auto f2 = [&] (const std::vector<stan::math::matrix_v>& x) {                                        \
      return FUN(TIME, AMT, RATE, II, EVID, CMT, ADDL, SS, x, BIOVAR, TLAG);                            \
    };                                                                                                  \
    torsten::test::test_grad(f1, f2, THETA, H, EPS_VAL, EPS_RTOL, EPS_ATOL);                            \
  }

#define TORSTEN_CPT_GRAD_BIOVAR_TEST(FUN, TIME, AMT, RATE, II, EVID, CMT, ADDL, SS, THETA, BIOVAR, TLAG, \
                                    H, EPS_VAL, EPS_RTOL, EPS_ATOL)                                      \
  {                                                                                                      \
    auto f1 = [&] (const std::vector<std::vector<double> >& x) {                                         \
      return FUN(TIME, AMT, RATE, II, EVID, CMT, ADDL, SS, THETA, x, TLAG);                              \
    };                                                                                                   \
    auto f2 = [&] (const std::vector<std::vector<stan::math::var> >& x) {                                \
      return FUN(TIME, AMT, RATE, II, EVID, CMT, ADDL, SS, THETA, x, TLAG);                              \
    };                                                                                                   \
    torsten::test::test_grad(f1, f2, BIOVAR, H, EPS_VAL, EPS_RTOL, EPS_ATOL);                            \
  }

#define TORSTEN_CPT_GRAD_TLAG_TEST(FUN, TIME, AMT, RATE, II, EVID, CMT, ADDL, SS, THETA, BIOVAR, TLAG,   \
                                    H, EPS_VAL, EPS_RTOL, EPS_ATOL)                                      \
  {                                                                                                      \
    auto f1 = [&] (const std::vector<std::vector<double> >& x) {                                         \
      return FUN(TIME, AMT, RATE, II, EVID, CMT, ADDL, SS, THETA, BIOVAR, x);                            \
    };                                                                                                   \
    auto f2 = [&] (const std::vector<std::vector<stan::math::var> >& x) {                                \
      return FUN(TIME, AMT, RATE, II, EVID, CMT, ADDL, SS, THETA, BIOVAR, x);                            \
    };                                                                                                   \
    torsten::test::test_grad(f1, f2, TLAG, H, EPS_VAL, EPS_RTOL, EPS_ATOL);                              \
  }

#define TORSTEN_CPT_GRAD_RATE_TEST(FUN, TIME, AMT, RATE, II, EVID, CMT, ADDL, SS, THETA, BIOVAR, TLAG,   \
                                    H, EPS_VAL, EPS_RTOL, EPS_ATOL)                                      \
  {                                                                                                      \
    auto f1 = [&] (const std::vector<double>& x) {                                                       \
      return FUN(TIME, AMT, x, II, EVID, CMT, ADDL, SS, THETA, BIOVAR, TLAG);                            \
    };                                                                                                   \
    auto f2 = [&] (const std::vector<stan::math::var>& x) {                                              \
      return FUN(TIME, AMT, x, II, EVID, CMT, ADDL, SS, THETA, BIOVAR, TLAG);                            \
    };                                                                                                   \
    torsten::test::test_grad(f1, f2, RATE, H, EPS_VAL, EPS_RTOL, EPS_ATOL);                              \
  }

#define TORSTEN_CPT_ODE_GRAD_TEST(FUN_CPT, FUN_ODE, F,                                                   \
                                    NCMT, TIME, AMT, RATE, II, EVID, CMT, ADDL, SS, THETA, BIOVAR, TLAG, \
                                    EPS_VAL, EPS_GRAD)                                                   \
  {                                                                                                      \
    {                                                                                                    \
      auto theta_v = torsten::to_var(THETA);                                                             \
      auto x1 = FUN_CPT(TIME, AMT, RATE, II, EVID, CMT, ADDL, SS, theta_v, BIOVAR, TLAG);                \
      auto x2 = FUN_ODE(F, NCMT, TIME, AMT, RATE, II, EVID, CMT, ADDL, SS, theta_v, BIOVAR, TLAG, nullptr); \
      for (size_t i = 0; i < theta_v.size(); ++i) {                                                      \
        torsten::test::test_grad(theta_v[i], x1, x2, EPS_VAL, EPS_GRAD);                                 \
      }                                                                                                  \
    }                                                                                                    \
    {                                                                                                    \
      auto biovar_v = torsten::to_var(BIOVAR);                                                           \
      auto x1 = FUN_CPT(TIME, AMT, RATE, II, EVID, CMT, ADDL, SS, THETA, biovar_v, TLAG);                \
      auto x2 = FUN_ODE(F, NCMT, TIME, AMT, RATE, II, EVID, CMT, ADDL, SS, THETA, biovar_v, TLAG, nullptr);       \
      for (size_t i = 0; i < biovar_v.size(); ++i) {                                                     \
        torsten::test::test_grad(biovar_v[i], x1, x2, EPS_VAL, EPS_GRAD);                                \
      }                                                                                                  \
    }                                                                                                    \
    {                                                                                                    \
      auto tlag_v = torsten::to_var(TLAG);                                                               \
      auto x1 = FUN_CPT(TIME, AMT, RATE, II, EVID, CMT, ADDL, SS, THETA, BIOVAR, tlag_v);                \
      auto x2 = FUN_ODE(F, NCMT, TIME, AMT, RATE, II, EVID, CMT, ADDL, SS, THETA, BIOVAR, tlag_v, nullptr);       \
      for (size_t i = 0; i < tlag_v.size(); ++i) {                                                       \
        torsten::test::test_grad(tlag_v[i], x1, x2, EPS_VAL, EPS_GRAD);                                  \
      }                                                                                                  \
    }                                                                                                    \
  }

#define TORSTEN_ODE_GRAD_THETA_TEST(FUN, F, NCMT, TIME, AMT, RATE, II, EVID, CMT, ADDL, SS, THETA, BIOVAR, TLAG,  \
                                   RTOL, ATOL, NSTEP,                                                             \
                                   H, EPS_VAL, EPS_RTOL, EPS_ATOL)                                                \
  {                                                                                                               \
    auto f1 = [&] (const std::vector<std::vector<double> >& x) {                                                  \
      return FUN(F, NCMT, TIME, AMT, RATE, II, EVID, CMT, ADDL, SS, x, BIOVAR, TLAG, RTOL, ATOL, NSTEP, nullptr);       \
    };                                                                                                            \
    auto f2 = [&] (const std::vector<std::vector<stan::math::var> >& x) {                                         \
      return FUN(F, NCMT, TIME, AMT, RATE, II, EVID, CMT, ADDL, SS, x, BIOVAR, TLAG, RTOL, ATOL, NSTEP, nullptr);       \
    };                                                                                                            \
    torsten::test::test_grad(f1, f2, THETA, H, EPS_VAL, EPS_RTOL, EPS_ATOL);                                      \
  }

#define TORSTEN_ODE_GRAD_BIOVAR_TEST(FUN, F, NCMT, TIME, AMT, RATE, II, EVID, CMT, ADDL, SS, THETA, BIOVAR, TLAG, \
                                   RTOL, ATOL, NSTEP,                                                             \
                                   H, EPS_VAL, EPS_RTOL, EPS_ATOL)                                                \
  {                                                                                                               \
    auto f1 = [&] (const std::vector<std::vector<double> >& x) {                                                  \
      return FUN(F, NCMT, TIME, AMT, RATE, II, EVID, CMT, ADDL, SS, THETA, x, TLAG, RTOL, ATOL, NSTEP, nullptr);        \
    };                                                                                                            \
    auto f2 = [&] (const std::vector<std::vector<stan::math::var> >& x) {                                         \
      return FUN(F, NCMT, TIME, AMT, RATE, II, EVID, CMT, ADDL, SS, THETA, x, TLAG, RTOL, ATOL, NSTEP, nullptr);        \
    };                                                                                                            \
    torsten::test::test_grad(f1, f2, BIOVAR, H, EPS_VAL, EPS_RTOL, EPS_ATOL);                                     \
  }


#define TORSTEN_ODE_GRAD_TLAG_TEST(FUN, F, NCMT, TIME, AMT, RATE, II, EVID, CMT, ADDL, SS, THETA, BIOVAR, TLAG,   \
                                   RTOL, ATOL, NSTEP,                                                             \
                                   H, EPS_VAL, EPS_RTOL, EPS_ATOL)                                                \
  {                                                                                                               \
    auto f1 = [&] (const std::vector<std::vector<double> >& x) {                                                  \
      return FUN(F, NCMT, TIME, AMT, RATE, II, EVID, CMT, ADDL, SS, THETA, BIOVAR, x, RTOL, ATOL, NSTEP, nullptr);      \
    };                                                                                                            \
    auto f2 = [&] (const std::vector<std::vector<stan::math::var> >& x) {                                         \
      return FUN(F, NCMT, TIME, AMT, RATE, II, EVID, CMT, ADDL, SS, THETA, BIOVAR, x, RTOL, ATOL, NSTEP, nullptr);      \
    };                                                                                                            \
    torsten::test::test_grad(f1, f2, TLAG, H, EPS_VAL, EPS_RTOL, EPS_ATOL);                                       \
  }

#define TORSTEN_ODE_GRAD_RATE_TEST(FUN, F, NCMT, TIME, AMT, RATE, II, EVID, CMT, ADDL, SS, THETA, BIOVAR, TLAG,   \
                                   RTOL, ATOL, NSTEP,                                                             \
                                    H, EPS_VAL, EPS_RTOL, EPS_ATOL)                                               \
  {                                                                                                               \
    auto f1 = [&] (const std::vector<double>& x) {                                                                \
      return FUN(F, NCMT, TIME, AMT, x, II, EVID, CMT, ADDL, SS, THETA, BIOVAR, TLAG, RTOL, ATOL, NSTEP, nullptr);      \
    };                                                                                                            \
    auto f2 = [&] (const std::vector<stan::math::var>& x) {                                                       \
      return FUN(F, NCMT, TIME, AMT, x, II, EVID, CMT, ADDL, SS, THETA, BIOVAR, TLAG, RTOL, ATOL, NSTEP, nullptr);      \
    };                                                                                                            \
    torsten::test::test_grad(f1, f2, RATE, H, EPS_VAL, EPS_RTOL, EPS_ATOL);                                       \
  }

#define TORSTEN_ODE_GRAD_AMT_TEST(FUN, F, NCMT, TIME, AMT, RATE, II, EVID, CMT, ADDL, SS, THETA, BIOVAR, TLAG,    \
                                   RTOL, ATOL, NSTEP,                                                             \
                                    H, EPS_VAL, EPS_RTOL, EPS_ATOL)                                               \
  {                                                                                                               \
    auto f1 = [&] (const std::vector<double>& x) {                                                                \
      return FUN(F, NCMT, TIME, x, RATE, II, EVID, CMT, ADDL, SS, THETA, BIOVAR, TLAG, RTOL, ATOL, NSTEP, nullptr);     \
    };                                                                                                            \
    auto f2 = [&] (const std::vector<stan::math::var>& x) {                                                       \
      return FUN(F, NCMT, TIME, x, RATE, II, EVID, CMT, ADDL, SS, THETA, BIOVAR, TLAG, RTOL, ATOL, NSTEP, nullptr);     \
    };                                                                                                            \
    torsten::test::test_grad(f1, f2, AMT, H, EPS_VAL, EPS_RTOL, EPS_ATOL);                                        \
  }

#define TORSTEN_ODE_GRAD_TIME_TEST(FUN, F, NCMT, TIME, AMT, RATE, II, EVID, CMT, ADDL, SS, THETA, BIOVAR, TLAG,   \
                                   RTOL, ATOL, NSTEP,                                                             \
                                   H, EPS_VAL, EPS_RTOL, EPS_ATOL)                                                \
  {                                                                                                               \
    auto f1 = [&] (const std::vector<double>& x) {                                                                \
      return FUN(F, NCMT, x, AMT, RATE, II, EVID, CMT, ADDL, SS, THETA, BIOVAR, TLAG, RTOL, ATOL, NSTEP, nullptr);      \
    };                                                                                                            \
    auto f2 = [&] (const std::vector<stan::math::var>& x) {                                                       \
      return FUN(F, NCMT, x, AMT, RATE, II, EVID, CMT, ADDL, SS, THETA, BIOVAR, TLAG, RTOL, ATOL, NSTEP, nullptr);      \
    };                                                                                                            \
    torsten::test::test_grad(f1, f2, TIME, H, EPS_VAL, EPS_RTOL, EPS_ATOL);                                       \
  }

#define TORSTEN_ODE_GRAD_II_TEST(FUN, F, NCMT, TIME, AMT, RATE, II, EVID, CMT, ADDL, SS, THETA, BIOVAR, TLAG,     \
                                   RTOL, ATOL, NSTEP,                                                             \
                                    H, EPS_VAL, EPS_RTOL, EPS_ATOL)                                               \
   {                                                                                                              \
     auto f1 = [&] (const std::vector<double>& x) {                                                               \
       return FUN(F, NCMT, TIME, AMT, RATE, x, EVID, CMT, ADDL, SS, THETA, BIOVAR, TLAG, RTOL, ATOL, NSTEP, nullptr);   \
     };                                                                                                           \
     auto f2 = [&] (const std::vector<stan::math::var>& x) {                                                      \
       return FUN(F, NCMT, TIME, AMT, RATE, x, EVID, CMT, ADDL, SS, THETA, BIOVAR, TLAG, RTOL, ATOL, NSTEP, nullptr);   \
     };                                                                                                           \
     torsten::test::test_grad(f1, f2, II, H, EPS_VAL, EPS_RTOL, EPS_ATOL);                                        \
   }

#endif
