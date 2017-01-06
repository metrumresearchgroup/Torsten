#ifndef TEST_UNIT_MATH_TORSTEN_REV_UTIL_TORSTEN_HPP
#define TEST_UNIT_MATH_TORSTEN_REV_UTIL_TORSTEN_HPP

#include <stan/math/rev/mat.hpp>
#include <gtest/gtest.h>
#include <test/unit/math/rev/arr/util.hpp>
#include <test/unit/util.hpp>

/*
 * Calculates finite difference for PKModelOneCpt with varying parameters. 
 */
Eigen::Matrix <double, Eigen::Dynamic, Eigen::Dynamic>
finite_diff_params(const std::vector<std::vector<double> >& pMatrix,
                   const std::vector<double>& time,
                   const std::vector<double>& amt,
                   const std::vector<double>& rate,
                   const std::vector<double>& ii,
                   const std::vector<int>& evid,
                   const std::vector<int>& cmt,
                   const std::vector<int>& addl,
                   const std::vector<int>& ss,
                   const size_t& param_row,
                   const size_t& param_col,
                   const double& diff) {
  using std::vector;
  using Eigen::Matrix;
  using Eigen::Dynamic;

  vector<double> parameters(pMatrix[0].size());
  vector<vector<double> > pMatrix_ub(pMatrix.size(), parameters);
  vector<vector<double> > pMatrix_lb(pMatrix.size(), parameters);
  for (size_t i = 0; i < pMatrix.size(); i++)
    for (size_t j = 0; j < pMatrix[0].size(); j++) {
      if (i == param_row && j == param_col) {        
        pMatrix_ub[i][j] = pMatrix[i][j] + diff;
        pMatrix_lb[i][j] = pMatrix[i][j] - diff;
      } else {
        pMatrix_ub[i][j] = pMatrix[i][j];
        pMatrix_lb[i][j] = pMatrix[i][j];
      }
    }

  Matrix<double, Dynamic, Dynamic> pk_res_ub;
  Matrix<double, Dynamic, Dynamic> pk_res_lb;
  pk_res_ub = PKModelTwoCpt(pMatrix_ub, time, amt, rate, ii, evid, cmt,
                            addl, ss);
  pk_res_lb = PKModelTwoCpt(pMatrix_lb, time, amt, rate, ii, evid, cmt,
                            addl, ss);

  return (pk_res_ub - pk_res_lb) / (2 * diff);
}

/*
 * Test PKModelTwoCpt with only pMatrix as vars and all other continuous
 * arguments as double.
 * Note: There is known issue when computing the derivative w.r.t the
 * lag time of a dosing compartment. The issue is reported on GitHub,
 * and the unit test overlooks it.
 */
void test_PKModelTwoCpt_finite_diff_v(
    const std::vector<std::vector<double> >& pMatrix,
    const std::vector<double>& time,
    const std::vector<double>& amt,
    const std::vector<double>& rate,
    const std::vector<double>& ii,
    const std::vector<int>& evid,
    const std::vector<int>& cmt,
    const std::vector<int>& addl,
    const std::vector<int>& ss,
    const double& diff,
    const double& diff2) {
  using std::vector;
  using Eigen::Matrix;
  using Eigen::Dynamic;
  using stan::math::var;

  size_t parmRows = pMatrix.size();
  size_t parmCols = pMatrix[0].size();
  size_t total_param = parmRows * parmCols;
  vector<vector<Matrix<double, Dynamic, Dynamic> > > finite_diff_res(parmRows);
  for(size_t i =  0; i < parmRows; i++)
    finite_diff_res[i].resize(parmCols);

  for(size_t i = 0; i < parmRows; i++) {
    for(size_t j = 0; j < parmCols; j++) {
      finite_diff_res[i][j] = finite_diff_params(pMatrix, time, amt, rate,
                                      ii, evid, cmt, addl, ss, i, j, diff);
    }
  }

  // Create pMatrix with vars
  vector<var> parameters(total_param);
  vector<vector<var> > pMatrix_v(parmRows);
  for (size_t i = 0; i < parmRows; i++) pMatrix_v[i].resize(parmCols);
  for (size_t i = 0; i < parmRows; i++) {
    for (size_t j = 0; j < parmCols; j++) {
      parameters[i * parmCols + j] = pMatrix[i][j];
      pMatrix_v[i][j] = parameters[i * parmCols + j];
    }
  }

  Matrix<var, Dynamic, Dynamic> ode_res;
  ode_res = PKModelTwoCpt(pMatrix_v, time, amt, rate, ii, evid, cmt, addl, ss);

  int nCmt = 3;
  size_t nEvent = time.size();

  // Identify dosing compartment
  vector<size_t> tlagIndexes(nCmt);
  for (int i = 0; i < nCmt; i++)
    tlagIndexes[i] = parmCols - nCmt + i;
  vector<bool> isDosingCmt(nCmt);
  for (size_t i = 0; i < nEvent; i++)
    if (evid[i] == 1 || evid[i] == 4) isDosingCmt[cmt[i] - 1] = true;

  vector<double> grads_eff(nEvent * nCmt);
  for (size_t i = 0; i < nEvent; i++)
    for (int j = 0; j < nCmt; j++) {
      grads_eff.clear();
      ode_res(i, j).grad(parameters, grads_eff);

      for (size_t k = 0; k < parmRows; k++)
        for (size_t l = 0; l < parmCols; l++) {

         bool discontinuous = false;
         for (int m = 0; m < nCmt; m++)
           if (l == tlagIndexes[m] && isDosingCmt[m]) discontinuous = true;

          if (discontinuous == false) {
            EXPECT_NEAR(grads_eff[k * parmCols + l],
              finite_diff_res[k][l](i, j), diff2)
              << "Gradient of PKModelTwoCpt failed with known"
              << " time, amt, rate, ii, evid, cmt, addl, ss "
              << " and unknown parameters at event " << i
              << ", in compartment " << j
              << ", and parameter index (" << k << ", " << l << ")";
          }
        }
      stan::math::set_zero_all_adjoints();
    }
}

void test_PKModelTwoCpt(const std::vector<std::vector<double> >& pMatrix,
                        const std::vector<double>& time,
                        const std::vector<double>& amt,
                        const std::vector<double>& rate,
                        const std::vector<double>& ii,
                        const std::vector<int>& evid,
                        const std::vector<int>& cmt,
                        const std::vector<int>& addl,
                        const std::vector<int>& ss,
                        const double& diff,
                        const double& diff2) {
  test_PKModelTwoCpt_finite_diff_v(pMatrix, time, amt, rate, ii, evid,
                                   cmt, addl, ss, diff, diff2);
}

void test_PKModelTwoCpt(const std::vector<double>& pMatrix_v,
                        const std::vector<double>& time,
                        const std::vector<double>& amt,
                        const std::vector<double>& rate,
                        const std::vector<double>& ii,
                        const std::vector<int>& evid,
                        const std::vector<int>& cmt,
                        const std::vector<int>& addl,
                        const std::vector<int>& ss,
                        const double& diff,
                        const double& diff2) {
  std::vector<std::vector<double> > pMatrix(1, pMatrix_v);
  test_PKModelTwoCpt(pMatrix, time, amt, rate, ii, evid, cmt, addl, ss,
                     diff, diff2);
}


// More tests
// test_ode_error_conditions
// test_ode_error_conditions_nan
// test_ode_error_conditions_inf
// test_ode_error_conditions_vd

#endif
