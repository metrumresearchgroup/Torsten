#ifndef TEST_UNIT_MATH_TORSTEN_UTIL_LINODEMODEL2_HPP
#define TEST_UNIT_MATH_TORSTEN_UTIL_LINODEMODEL2_HPP

#include <stan/math/rev/mat.hpp>
#include <gtest/gtest.h>
#include <test/unit/math/rev/arr/util.hpp>
#include <test/unit/util.hpp>

/**
 * Calculates finite difference for generalOdeModel with varying parameters.
 * Parameters can be stored in biovar and/or system. Need to specify
 * the row and column of the parameter in these object. For system, also
 * need to specify the event number (recall system is a vector of matrix).  
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
                   const std::vector< Eigen::Matrix<double, Eigen::Dynamic,
                     Eigen::Dynamic> >& system,
                   const std::vector<std::vector<double> >& biovar,
                   const std::vector<std::vector<double> >& tlag,
                   const size_t& param_row,
                   const size_t& param_col,
		               const size_t& param_z,
                   const double& diff,
                   const std::string parmType) { 
  using std::vector;
  using Eigen::Matrix;
  using Eigen::Dynamic;

  Matrix<double, Dynamic, Dynamic> system_1(system[0].rows(), system[0].cols());
  vector<Matrix<double, Dynamic, Dynamic> > system_ub(system.size(), system_1);
  vector<Matrix<double, Dynamic, Dynamic> > system_lb(system.size(), system_1);
  for (size_t z = 0; z < system.size(); z++)
    for (int i = 0; i < system[0].rows(); i++)
      for (int j = 0; j < system[0].cols(); j++) {      
        if (((size_t) i == param_row
             && (size_t) j == param_col)
             && (z == param_z
             && parmType == "system")) {
	      system_ub[z](i, j) = system[z](i, j) + diff;
	      system_lb[z](i, j) = system[z](i, j) - diff;
        } else {
	      system_ub[z](i, j) = system[z](i, j);
	      system_lb[z](i, j) = system[z](i, j);
        }
      }
	
	vector<double> bio_para(biovar[0].size());
  vector<vector<double> > biovar_ub(biovar.size(), bio_para);
  vector<vector<double> > biovar_lb(biovar.size(), bio_para);
  for (size_t i = 0; i < biovar.size(); i++)
    for (size_t j = 0; j < biovar[i].size(); j++) {
      if ((i == param_row && j == param_col) && parmType == "biovar") {        
        biovar_ub[i][j] = biovar[i][j] + diff;
        biovar_lb[i][j] = biovar[i][j] - diff;
      } else {
        biovar_ub[i][j] = biovar[i][j];
        biovar_lb[i][j] = biovar[i][j];
      }
    }
    
  vector<double> tlag_par(tlag[0].size());
  vector<vector<double> > tlag_ub(tlag.size(), tlag_par);
  vector<vector<double> > tlag_lb(tlag.size(), tlag_par);
  for (size_t i = 0; i < tlag.size(); i++)
    for (size_t j = 0; j < tlag[i].size(); j++) {
      if ((i == param_row && j == param_col) && parmType == "tlag") {
        tlag_ub[i][j] = tlag[i][j] + diff;
        tlag_lb[i][j] = tlag[i][j] - diff;
      } else {
        tlag_ub[i][j] = tlag[i][j];
        tlag_lb[i][j] = tlag[i][j];
      }
    }

  Matrix<double, Dynamic, Dynamic> pk_res_ub;
  Matrix<double, Dynamic, Dynamic> pk_res_lb;
  pk_res_ub = torsten::pmx_solve_linode(time, amt, rate, ii, evid, cmt, addl, ss,
                          system_ub, biovar_ub, tlag_ub);
  pk_res_lb = torsten::pmx_solve_linode(time, amt, rate, ii, evid, cmt, addl, ss,
                          system_lb, biovar_lb, tlag_lb);

  return (pk_res_ub - pk_res_lb) / (2 * diff);
}

/**
 * Test linearOdeModel with only system as vars and all other
 * continuous arguments as double.
 */
void test_linearOdeModel_finite_diff_vdd(
     const std::vector<double>& time,
     const std::vector<double>& amt,
     const std::vector<double>& rate,
     const std::vector<double>& ii,
     const std::vector<int>& evid,
     const std::vector<int>& cmt,
     const std::vector<int>& addl,
     const std::vector<int>& ss,
     const std::vector<Eigen::Matrix<double, Eigen::Dynamic,
                                     Eigen::Dynamic> >& system,
     const std::vector<std::vector<double> >& biovar,
     const std::vector<std::vector<double> >& tlag,
     const double diff,
     const double diff2) {
  using std::vector;
  using Eigen::Matrix;
  using Eigen::Dynamic;
  using stan::math::var;

  size_t nCmt = system[0].rows();
  size_t nSystems = system.size();

  // Size object to store finite diff results
  vector<vector<vector<Matrix<double, Dynamic, Dynamic> > > >
    finite_diff_res(nSystems);
  for (size_t z = 0; z < nSystems; z++) {
    finite_diff_res[z].resize(nCmt);
    for (size_t i = 0; i < nCmt; i++)
      finite_diff_res[z][i].resize(nCmt);
  }

  // Compute finite diff results
  for (size_t z = 0; z < nSystems; z++)
    for (size_t k = 0; k < nCmt; k++)
      for (size_t l = 0; l < nCmt; l++)
        finite_diff_res[z][k][l]
          = finite_diff_params(time, amt, rate, ii, evid, cmt, addl, ss,
                               system, biovar, tlag, k, l, z, diff, "system");

  std::vector<double> grads_eff;

  // Create system with vars instead of doubles
  size_t totalParam = nCmt * nCmt * nSystems;
  vector<var> parameters(totalParam);
  vector<Matrix<var, Dynamic, Dynamic> > system_v(nSystems);
  for (size_t z = 0; z < nSystems; z++) {
    system_v[z].resize(nCmt, nCmt);
    for (size_t i = 0; i < nCmt; i++)
      for (size_t j = 0; j < nCmt; j++) {
        parameters[z * nCmt * nCmt + i * nCmt + j] = system[z](i, j);
        system_v[z](i, j) = parameters[z * nCmt * nCmt + i * nCmt + j];
      }
  }

  // Compute return of pmx_solve_linode
  Matrix<var, Dynamic, Dynamic> ode_res
    = torsten::pmx_solve_linode(time, amt, rate, ii, evid, cmt, addl, ss,
                  system_v, biovar, tlag);

  // Test auto-diff
  size_t nEvents = time.size();
  for (size_t i = 0; i < nEvents; i++)
    for (size_t j = 0; j < nCmt; j++) {
      grads_eff.clear();
      ode_res(i, j).grad(parameters, grads_eff);

  for (size_t z = 0; z < nSystems; z++)
    for (size_t k = 0; k < nCmt; k++)
      for (size_t l = 0; l < nCmt; l++)
        EXPECT_NEAR(grads_eff[z * nCmt * nCmt + k * nCmt + l],
          finite_diff_res[z][k][l](i, j), diff2)
        << "Gradient of generalOdeModel failed with known"
        << " event data, biovar, and tlag, "
        << " and unknown parameters at event " << i
        << ", in compartment " << j
        << ", and parameter index (" << z << ", " << k << ", "
        << l << ")";

        stan::math::set_zero_all_adjoints();
    }
}

/**
 * Test linearOdeModel with only biovar as vars and all other 
 * continuous arguments as double.
 */ 
void test_linearOdeModel_finite_diff_dvd(
    const std::vector<double>& time,
    const std::vector<double>& amt,
    const std::vector<double>& rate,
    const std::vector<double>& ii,
    const std::vector<int>& evid,
    const std::vector<int>& cmt,
    const std::vector<int>& addl,
    const std::vector<int>& ss,
    const std::vector<Eigen::Matrix<double, Eigen::Dynamic,
      Eigen::Dynamic> >& system,
    const std::vector<std::vector<double> >& biovar,
    const std::vector<std::vector<double> >& tlag,
    const double& diff,
    const double& diff2) {
  using std::vector;
  using Eigen::Matrix;
  using Eigen::Dynamic;
  using stan::math::var;

  size_t parmRows = biovar.size();
  size_t parmCols = biovar[0].size();
  size_t total_param = parmRows * parmCols;
  vector<vector<Matrix<double, Dynamic, Dynamic> > > finite_diff_res(parmRows);
  for(size_t i =  0; i < parmRows; i++)
    finite_diff_res[i].resize(parmCols);

  for(size_t i = 0; i < parmRows; i++) {
    for(size_t j = 0; j < parmCols; j++) {
      finite_diff_res[i][j] = finite_diff_params(time, amt, rate, ii, evid,
                                                 cmt, addl, ss, system,
                                                 biovar, tlag, i, j, 0,
                                                 diff, "biovar");
    }
  }

  // Create biovar with vars
  vector<var> parameters(total_param);
  vector<vector<var> > biovar_v(parmRows);
  for (size_t i = 0; i < parmRows; i++) {
    biovar_v[i].resize(parmCols);
    for (size_t j = 0; j < parmCols; j++) {
      parameters[i * parmCols + j] = biovar[i][j];
      biovar_v[i][j] = parameters[i * parmCols + j];
    }
  }

  Matrix<var, Dynamic, Dynamic> ode_res;
  ode_res = torsten::pmx_solve_linode(time, amt, rate, ii, evid, cmt, addl, ss,
                        system, biovar_v, tlag);

  int nCmt = system[0].cols();
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
            << "Gradient of generalOdeModel failed with known"
            << " event data, parameters, and tlag "
            << " and unknown biovar " << i
            << ", in compartment " << j
            << ", and parameter index (" << k << ", " << l << ")";
        }
      stan::math::set_zero_all_adjoints();
    }  
}

/**
 * Test linearOdeModel with only tlag as vars and all other 
 * continuous arguments as double.
 * Note: There is known issue when computing the derivative w.r.t the
 * lag time of a dosing compartment. The issue is reported on GitHub,
 * and the unit test overlooks it.
 */ 
void test_linearOdeModel_finite_diff_ddv(
    const std::vector<double>& time,
    const std::vector<double>& amt,
    const std::vector<double>& rate,
    const std::vector<double>& ii,
    const std::vector<int>& evid,
    const std::vector<int>& cmt,
    const std::vector<int>& addl,
    const std::vector<int>& ss,
    const std::vector<Eigen::Matrix<double, Eigen::Dynamic,
                                    Eigen::Dynamic> >& system,
    const std::vector<std::vector<double> >& biovar,
    const std::vector<std::vector<double> >& tlag,
    const double& diff,
    const double& diff2) {
  using std::vector;
  using Eigen::Matrix;
  using Eigen::Dynamic;
  using stan::math::var;

  size_t parmRows = tlag.size();
  size_t parmCols = tlag[0].size();
  size_t total_param = parmRows * parmCols;
  vector<vector<Matrix<double, Dynamic, Dynamic> > > finite_diff_res(parmRows);
  for(size_t i =  0; i < parmRows; i++)
    finite_diff_res[i].resize(parmCols);

  for(size_t i = 0; i < parmRows; i++) {
    for(size_t j = 0; j < parmCols; j++) {
      finite_diff_res[i][j] = finite_diff_params(time, amt, rate, ii, evid,
                                                 cmt, addl, ss, system,
                                                 biovar, tlag, i, j, 0,
                                                 diff, "tlag");
    }
  }
  
  // Create tlag with vars
  vector<var> parameters(total_param);
  vector<vector<var> > tlag_v(parmRows);
  for (size_t i = 0; i < parmRows; i++) {
    tlag_v[i].resize(parmCols);
    for (size_t j = 0; j < parmCols; j++) {
      parameters[i * parmCols + j] = tlag[i][j];
      tlag_v[i][j] = parameters[i * parmCols + j];
    }
  }
  
  Matrix<var, Dynamic, Dynamic> ode_res;
  ode_res = torsten::pmx_solve_linode(time, amt, rate, ii, evid, cmt, addl, ss,
                        system, biovar, tlag_v);
  
  int nCmt = system[0].cols();
  size_t nEvent = time.size();
  
  // Identify dosing compartment
  vector<bool> isDosingCmt(nCmt);
  for (size_t i = 0; i < nEvent; i++)
    if (evid[i] == 1 || evid[i] == 4) isDosingCmt[cmt[i] - 1] = true;

  vector<double> grads_eff(nEvent * nCmt);
  for (size_t i = 0; i < nEvent; i++)
    for (int j = 0; j < nCmt; j++) {
      grads_eff.clear();
      ode_res(i, j).grad(parameters, grads_eff);

      for (size_t k = 0; k < parmRows; k++) {
        for (size_t l = 0; l < parmCols; l++) {
          double tlag = parameters[k * parmCols + l].val();
          bool skip = false;

          // When tlag is zero, all the gradients w.r.t to tlag go to
          // 0, because of an if (tlag == 0) statement. This causes an
          // error if the lag time is in a dosing comopartment.
          if (tlag == 0)
            if (isDosingCmt[l]) skip = true;

          // When tlag is non-zero, the gradient will not properly
          // evaluate at an event which coincides with the time of
          // the dosing (lag time accounted for), in the dosing
          // compartment. The following IF cascade identifies such
          // entries.
          if (tlag != 0)
            for (size_t m = 0; m < nEvent; m++) {
              if ((evid[m] == 1 || evid[m] == 4)
                  && ((cmt[m] - 1) == (int) l)) {
                if ((time[m] + tlag) == time[i]) {
                  skip = true;
                }
              }
            }

          if (skip == false) {
            EXPECT_NEAR(grads_eff[k * parmCols + l],
                        finite_diff_res[k][l](i, j), diff2)
              << "Gradient of generalOdeModel failed with known"
              << " event data, parameters, and biovar, "
              << " and lag times at event " << i
              << ", in compartment " << j
              << ", and parameter index (" << k << ", " << l << ")";
          }
        }
      }
      stan::math::set_zero_all_adjoints();
    }  
}    

void test_pmx_solve_linode(const std::vector<double>& time,
                      const std::vector<double>& amt,
                      const std::vector<double>& rate,
                      const std::vector<double>& ii,
                      const std::vector<int>& evid,
                      const std::vector<int>& cmt,
                      const std::vector<int>& addl,
                      const std::vector<int>& ss,
                      const std::vector<Eigen::Matrix<double, Eigen::Dynamic,
                        Eigen::Dynamic> > system,
                      const std::vector<std::vector<double> >& biovar,
                      const std::vector<std::vector<double> >& tlag,
                      const double& diff,
                      const double& diff2) {
  test_linearOdeModel_finite_diff_vdd(time, amt, rate, ii, evid, cmt, addl, ss,
                                      system, biovar, tlag, diff, diff2);

  test_linearOdeModel_finite_diff_dvd(time, amt, rate, ii, evid, cmt, addl, ss,
                                      system, biovar, tlag, diff, diff2);

  test_linearOdeModel_finite_diff_ddv(time, amt, rate, ii, evid, cmt, addl, ss,
                                      system, biovar, tlag, diff, diff2);
}


// More tests
// test_ode_error_conditions
// test_ode_error_conditions_nan
// test_ode_error_conditions_inf
// test_ode_error_conditions_vd

#endif
