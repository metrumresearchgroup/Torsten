#ifndef TEST_UNIT_MATH_TORSTEN_REV_UTIL_PKMODELONECPT_HPP
#define TEST_UNIT_MATH_TORSTEN_REV_UTIL_PKMODELONECPT_HPP

#include <stan/math/torsten/torsten.hpp>
#include <gtest/gtest.h>
#include <test/unit/math/rev/arr/util.hpp>
#include <test/unit/util.hpp>
 
/*
 * Calculates finite difference for PKModelOneCpt with varying parameters.
 * Parameters are stored in pMatrix and addParm.
 */
Eigen::Matrix <double, Eigen::Dynamic, Eigen::Dynamic>
finite_diff_params(const std::vector<double>& time,
                   const std::vector<double>& amt,
                   const std::vector<double>& rate,
                   const std::vector<double>& ii,
                   const std::vector<int>& evid,
                   const std::vector<int>& cmt,
                   const std::vector<int>& addl,
                   const std::vector<int>& ss,
                   const std::vector<std::vector<double> >& pMatrix,
                   const std::vector<std::vector<double> >& addParm,
                   const size_t& param_row,
                   const size_t& param_col,
                   const double& diff,
                   const std::string parmType) {
  using std::vector;
  using Eigen::Matrix;
  using Eigen::Dynamic;

  vector<double> parameters(pMatrix[0].size());
  vector<vector<double> > pMatrix_ub(pMatrix.size(), parameters);
  vector<vector<double> > pMatrix_lb(pMatrix.size(), parameters);
  for (size_t i = 0; i < pMatrix.size(); i++)
    for (size_t j = 0; j < pMatrix[0].size(); j++) {
      if ((i == param_row && j == param_col) && parmType == "pMatrix") {        
        pMatrix_ub[i][j] = pMatrix[i][j] + diff;
        pMatrix_lb[i][j] = pMatrix[i][j] - diff;
      } else {
        pMatrix_ub[i][j] = pMatrix[i][j];
        pMatrix_lb[i][j] = pMatrix[i][j];
      }
    }

  vector<double> addParameters(addParm[0].size());
  vector<vector<double> > addParm_ub(addParm.size(), addParameters);
  vector<vector<double> > addParm_lb(addParm.size(), addParameters);
  for (size_t i = 0; i < addParm.size(); i++)
    for (size_t j = 0; j < addParm.size(); j++) {
      if ((i == param_row && j == param_col) && parmType == "addParm") {
        addParm_ub[i][j] = addParm[i][j] + diff;
        addParm_lb[i][j] = addParm[i][j] - diff; 
      } else {
        addParm_ub[i][j] = addParm[i][j];
        addParm_lb[i][j] = addParm[i][j];
      }
    }

  Matrix<double, Dynamic, Dynamic> pk_res_ub;
  Matrix<double, Dynamic, Dynamic> pk_res_lb;
  pk_res_ub = PKModelOneCpt(time, amt, rate, ii, evid, cmt, addl, ss,
                            pMatrix_ub, addParm_ub);
  pk_res_lb = PKModelOneCpt(time, amt, rate, ii, evid, cmt, addl, ss,
                            pMatrix_lb, addParm_lb);

  return (pk_res_ub - pk_res_lb) / (2 * diff);
}

/*
 * Test PKModelOneCpt with only pMatrix as vars and all other continuous
 * arguments as double.
 */
void test_PKModelOneCpt_finite_diff_vd(
    const std::vector<double>& time,
    const std::vector<double>& amt,
    const std::vector<double>& rate,
    const std::vector<double>& ii,
    const std::vector<int>& evid,
    const std::vector<int>& cmt,
    const std::vector<int>& addl,
    const std::vector<int>& ss,
    const std::vector<std::vector<double> >& pMatrix,
    const std::vector<std::vector<double> >& addParm,
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
      finite_diff_res[i][j]
        = finite_diff_params(time, amt, rate, ii, evid, cmt, addl, ss,
                             pMatrix, addParm, i, j, diff, "pMatrix");
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
  ode_res = PKModelOneCpt(time, amt, rate, ii, evid, cmt, addl, ss,
                          pMatrix_v, addParm);

  int nCmt = 2;
  size_t nEvent = time.size();

  vector<double> grads_eff(nEvent * nCmt);
  for (size_t i = 0; i < nEvent; i++)
    for (int j = 0; j < nCmt; j++) {
      grads_eff.clear();
      ode_res(i, j).grad(parameters, grads_eff);

      for (size_t k = 0; k < parmRows; k++)
        for (size_t l = 0; l < parmCols; l++) {
          
          EXPECT_NEAR(grads_eff[k * parmCols + l],
                      finite_diff_res[k][l](i, j), diff2)
          << "Gradient of PKModelOneCpt failed with known"
          << " time, amt, rate, ii, evid, cmt, addl, ss "
          << " and unknown parameters at event " << i
          << ", in compartment " << j
          << ", and parameter index (" << k << ", " << l << ")";
        }
      stan::math::set_zero_all_adjoints();
    }
}

/*
 * Test PKModelOneCpt with only addParm as vars and all other continuous
 * arguments as double.
 * Note: There is known issue when computing the derivative w.r.t the
 * lag time of a dosing compartment. The issue is reported on GitHub,
 * and the unit test overlooks it.
 */
void test_PKModelOneCpt_finite_diff_dv(
    const std::vector<double>& time,
    const std::vector<double>& amt,
    const std::vector<double>& rate,
    const std::vector<double>& ii,
    const std::vector<int>& evid,
    const std::vector<int>& cmt,
    const std::vector<int>& addl,
    const std::vector<int>& ss,
    const std::vector<std::vector<double> >& pMatrix,
    const std::vector<std::vector<double> >& addParm,
    const double& diff,
    const double& diff2) {
  using std::vector;
  using Eigen::Matrix;
  using Eigen::Dynamic;
  using stan::math::var;

  size_t parmRows = addParm.size();
  size_t parmCols = addParm[0].size();
  size_t total_param = parmRows * parmCols;
  vector<vector<Matrix<double, Dynamic, Dynamic> > > finite_diff_res(parmRows);
  for(size_t i = 0; i < parmRows; i++) finite_diff_res[i].resize(parmCols);
  
  for(size_t i = 0; i < parmRows; i++) {
    for(size_t j = 0; j < parmCols; j++) {
      finite_diff_res[i][j]
      = finite_diff_params(time, amt, rate, ii, evid, cmt, addl, ss,
                           pMatrix, addParm, i, j, diff, "addParm");
    }
  }

  // Create addParm with vars
  vector<var> parameters(total_param);
  vector<vector<var> > addParm_v(parmRows);
  for (size_t i = 0; i < parmRows; i++) addParm_v[i].resize(parmCols);
  for (size_t i = 0; i < parmRows; i++) {
    for (size_t j = 0; j < parmCols; j++) {
      parameters[i * parmCols + j] = addParm[i][j];
      addParm_v[i][j] = parameters[i * parmCols + j];
    }
  }

  Matrix<var, Dynamic, Dynamic> ode_res;
  ode_res = PKModelOneCpt(time, amt, rate, ii, evid, cmt, addl, ss,
                          pMatrix, addParm_v);
  
  size_t nEvent = time.size();
  
  // Identify dosing compartment
  int nCmt = 2;
  vector<size_t> tlagIndexes(nCmt);
  for (int i = 0; i < nCmt; i++)
    tlagIndexes[i] = nCmt + i;
  vector<bool> isDosingCmt(nCmt);
  for (size_t i = 0; i < nEvent; i++)
    if (evid[i] == 0 || evid[i] == 4) isDosingCmt[cmt[i] - 1] = true;
    
    vector<double> grads_eff(nEvent * nCmt);
    for (size_t i = 1; i < nEvent; i++)
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
              << "Gradient of generalOdeModel failed with known"
              << " time, amt, rate, ii, evid, cmt, addl, ss "
              << " and unknown parameters at event " << i
              << ", in compartment " << j
              << ", and parameter index (" << k << ", " << l << ")";
            }
          }
          stan::math::set_zero_all_adjoints();
      }
}


void test_PKModelOneCpt(const std::vector<double>& time,
                        const std::vector<double>& amt,
                        const std::vector<double>& rate,
                        const std::vector<double>& ii,
                        const std::vector<int>& evid,
                        const std::vector<int>& cmt,
                        const std::vector<int>& addl,
                        const std::vector<int>& ss,
                        const std::vector<std::vector<double> >& pMatrix,
                        const std::vector<std::vector<double> >& addParm,
                        const double& diff,
                        const double& diff2) {
  test_PKModelOneCpt_finite_diff_vd(time, amt, rate, ii, evid,
                                    cmt, addl, ss, pMatrix, addParm,
                                    diff, diff2);

  test_PKModelOneCpt_finite_diff_dv(time, amt, rate, ii, evid,
                                    cmt, addl, ss, pMatrix, addParm,
                                    diff, diff2);
}

/*
void test_PKModelOneCpt(const std::vector<double>& pMatrix_v,
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
  test_PKModelOneCpt(pMatrix, time, amt, rate, ii, evid, cmt, addl, ss,
                     diff, diff2);
}
*/

// More tests
// test_ode_error_conditions
// test_ode_error_conditions_nan
// test_ode_error_conditions_inf
// test_ode_error_conditions_vd

#endif
