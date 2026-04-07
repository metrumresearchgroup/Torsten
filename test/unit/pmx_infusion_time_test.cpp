#include <stan/math/prim/fun/Eigen.hpp>
#include <stan/math/torsten/test/unit/pmx_ode_test_fixture.hpp>
#include <stan/math/torsten/test/unit/pmx_onecpt_test_fixture.hpp>
#include <stan/math/torsten/test/unit/pmx_twocpt_test_fixture.hpp>
#include <stan/math/torsten/test/unit/pmx_friberg_karlsson_test_fixture.hpp>
#include <stan/math/torsten/test/unit/expect_near_matrix_eq.hpp>
#include <stan/math/torsten/test/unit/expect_matrix_eq.hpp>
#include <stan/math/torsten/mpi/session_def.hpp>
#include <stan/math/torsten/pmx_solve_rk45.hpp>
#include <stan/math/torsten/pmx_solve_bdf.hpp>
#include <stan/math/torsten/pmx_solve_adams.hpp>
#include <stan/math/torsten/pmx_twocpt_model.hpp>
#include <stan/math/torsten/pmx_onecpt_model.hpp>
#include <stan/math/torsten/pmx_ode_model.hpp>
#include <gtest/gtest.h>

auto f_onecpt = torsten::PMXOneCptModel<double>::f_;
auto f_twocpt = torsten::PMXTwoCptModel<double>::f_;

using stan::math::var;
using std::vector;
using Eigen::Matrix;
using Eigen::MatrixXd;
using Eigen::Dynamic;

using torsten::pmx_solve_rk45;
using torsten::pmx_solve_bdf;
using torsten::pmx_solve_adams;
using torsten::NONMENEventsRecord;
using torsten::NonEventParameters;
using std::istream;
using std::string;
using std::stringstream;
using std::vector;
using std::pow;
using namespace stan::math;
using stan::math::pow; 

template <typename T0__, typename T1__, typename T2__, typename T3__>
Eigen::Matrix<stan::promote_args_t<T0__, T1__, T2__, T3__>, -1, 1>
memBODE(const T0__& t, const Eigen::Matrix<T1__, -1, 1>& x,
        std::ostream* pstream__,
        const std::vector<T2__>& parms, const std::vector<T3__>& rdummy,
        const std::vector<int>& idummy) {
  using local_scalar_t__ = stan::promote_args_t<T0__, T1__, T2__, T3__>;
  const static bool propto__ = true;
  (void) propto__;
  local_scalar_t__ DUMMY_VAR__(std::numeric_limits<double>::quiet_NaN());
  (void) DUMMY_VAR__;  // suppress unused var warning
  
  try {
    local_scalar_t__ memB0;
    memB0 = DUMMY_VAR__;
    
    memB0 = parms[(1 - 1)];
    local_scalar_t__ kM;
    kM = DUMMY_VAR__;
    
    
    kM = parms[(2 - 1)];
    local_scalar_t__ kP;
    kP = DUMMY_VAR__;
    
    
    kP = parms[(3 - 1)];
    local_scalar_t__ Emax;
    Emax = DUMMY_VAR__;
    
    
    Emax = parms[(4 - 1)];
    local_scalar_t__ EC50;
    EC50 = DUMMY_VAR__;
    
    
    EC50 = parms[(5 - 1)];
    local_scalar_t__ EmaxP;
    EmaxP = DUMMY_VAR__;
    
    
    EmaxP = parms[(6 - 1)];
    local_scalar_t__ EC50P;
    EC50P = DUMMY_VAR__;
    
    
    EC50P = parms[(7 - 1)];
    local_scalar_t__ Isc;
    Isc = DUMMY_VAR__;
    
    
    Isc = parms[(8 - 1)];
    local_scalar_t__ CLIV;
    CLIV = DUMMY_VAR__;
    
    
    CLIV = 215;
    local_scalar_t__ QIV;
    QIV = DUMMY_VAR__;
    
    
    QIV = 459;
    local_scalar_t__ VcenIV;
    VcenIV = DUMMY_VAR__;
    
    
    VcenIV = 2560;
    local_scalar_t__ VperIV;
    VperIV = DUMMY_VAR__;
    
    
    VperIV = 2730;
    local_scalar_t__ kaSC;
    kaSC = DUMMY_VAR__;
    
    
    kaSC = 0.235;
    local_scalar_t__ CLSC;
    CLSC = DUMMY_VAR__;
    
    
    CLSC = 204;
    local_scalar_t__ QSC;
    QSC = DUMMY_VAR__;
    
    
    QSC = 683;
    local_scalar_t__ VcenSC;
    VcenSC = DUMMY_VAR__;
    
    
    VcenSC = 2290;
    local_scalar_t__ VperSC;
    VperSC = DUMMY_VAR__;
    
    
    VperSC = 2650;
    local_scalar_t__ F1SC;
    F1SC = DUMMY_VAR__;
    
    
    F1SC = 0.752;
    local_scalar_t__ tLag1SC;
    tLag1SC = DUMMY_VAR__;
    
    
    tLag1SC = 0.179;
    local_scalar_t__ CL;
    CL = DUMMY_VAR__;
    
    
    CL = (((1 - Isc) * CLIV) + (Isc * CLSC));
    local_scalar_t__ Q;
    Q = DUMMY_VAR__;
    
    
    Q = (((1 - Isc) * QIV) + (Isc * QSC));
    local_scalar_t__ Vcen;
    Vcen = DUMMY_VAR__;
    
    
    Vcen = (((1 - Isc) * VcenIV) + (Isc * VcenSC));
    local_scalar_t__ Vper;
    Vper = DUMMY_VAR__;
    
    
    Vper = (((1 - Isc) * VperIV) + (Isc * VperSC));
    local_scalar_t__ ka;
    ka = DUMMY_VAR__;
    
    
    ka = ((1 - Isc) + (Isc * kaSC));
    local_scalar_t__ x4;
    x4 = DUMMY_VAR__;
    
    
    x4 = stan::math::fmax(stan::math::machine_precision(),
           (x[(4 - 1)] + memB0));
    local_scalar_t__ x5;
    x5 = DUMMY_VAR__;
    
    
    x5 = (x[(5 - 1)] + 1.0);
    local_scalar_t__ Cp;
    Cp = DUMMY_VAR__;
    
    
    Cp = ((1000 * x[(2 - 1)]) / Vcen);
    local_scalar_t__ ST;
    ST = DUMMY_VAR__;
    
    
    ST = (1 + ((Emax * Cp) / (EC50 + Cp)));
    local_scalar_t__ STP;
    STP = DUMMY_VAR__;
    
    
    STP = (1 + ((EmaxP * Cp) / (EC50P + Cp)));
    Eigen::Matrix<local_scalar_t__, -1, 1> dxdt(5);
    
    dxdt[0] = (-ka * x[(1 - 1)]);
    dxdt[1] = (((ka * x[(1 - 1)]) - (((CL / Vcen) + (Q / Vcen)) * x[(2 - 1)])) +
               ((Q / Vper) * x[(3 - 1)]));
    dxdt[2] = (((Q / Vcen) * x[(2 - 1)]) - ((Q / Vper) * x[(3 - 1)]));
    dxdt[3] = (kM * (((memB0 * ST) * x5) - x4));
    dxdt[4] = (kP * (1 - (STP * x5)));
    
    return dxdt;
  } catch (const std::exception& e) {
    throw;
  }
}

struct memBODE_functor__ {
template <typename T0__, typename T1__, typename T2__, typename T3__>
Eigen::Matrix<stan::promote_args_t<T0__, T1__, T2__, T3__>, -1, 1>
operator()(const T0__& t, const Eigen::Matrix<T1__, -1, 1>& x,
           std::ostream* pstream__,
           const std::vector<T2__>& parms, const std::vector<T3__>& rdummy,
           const std::vector<int>& idummy)  const 
{
  return memBODE(t, x, pstream__, parms, rdummy, idummy);
}
};

TEST(pmx_solve_ode, solver_iv_end_time_equals_to_next_ev_time) {
  // 2 subjects, 3 events per subject
  using stan::math::var;
  int nId = 1;
  int nt = 6;
  int nCmt = 5;
  std::vector<std::vector<var> > parms{{0.0208766,45439.7,68249.2,0.784342,1.22549,1.17677,0.889028,0}};
  std::vector<double> amt{0, 51, 51, 51, 0, 51};
  std::vector<double> rate{0, 3.64285714285714, 3.64285714285714, 1.82142857142857, 0, 1.82142857142857};
  std::vector<int> cmt{4, 2, 2, 2, 4, 2};
  std::vector<int> evid{0, 1, 1, 1, 0, 1};
  std::vector<double> time{0, 0, 14, 28, 56, 56.1};
  std::vector<int> addl(nt, 0);
  std::vector<int> ss(nt, 0);
  std::vector<double> ii(nt, 0.0);

  EXPECT_NO_THROW(pmx_solve_bdf(memBODE_functor__(), nCmt, time, amt, rate, ii, evid, cmt, addl, ss, parms, nullptr));
}
