#include <stan/math/rev/mat.hpp>  // FIX ME - includes should be more specific
#include <test/unit/math/torsten/expect_near_matrix_eq.hpp>
#include <test/unit/math/torsten/expect_matrix_eq.hpp>
#include <test/unit/math/torsten/pmx_twocpt_mpi_test_fixture.hpp>
#include <test/unit/math/torsten/util_generalOdeModel.hpp>
#include <test/unit/math/torsten/test_util.hpp>
#include <stan/math/torsten/mpi/envionment.hpp>
#include <stan/math/torsten/pmx_solve_twocpt.hpp>
#include <stan/math/torsten/pmx_onecpt_model.hpp>
#include <stan/math/torsten/pmx_twocpt_model.hpp>
#include <stan/math/torsten/to_var.hpp>
#include <gtest/gtest.h>
#include <stan/math/rev/mat.hpp>  // FIX ME - include should be more specific
#include <test/unit/math/torsten/util_pmx_solve_twocpt.hpp>
#include <vector>

using std::vector;
using Eigen::Matrix;
using Eigen::MatrixXd;
using Eigen::Dynamic;
using stan::math::var;

TEST_F(TorstenPopulationPMXTwoCptTest, multiple_bolus_doses_data_only) {
  Matrix<double, Dynamic, Dynamic> x =
    torsten::pmx_solve_twocpt(time, amt, rate, ii, evid, cmt, addl, ss, pMatrix, biovar, tlag); // NOLINT

  Matrix<double, Dynamic, Dynamic> x_m =
    torsten::pmx_solve_group_twocpt(len, time_m, amt_m, rate_m, ii_m, evid_m, cmt_m, addl_m, ss_m, // NOLINT
                              pMatrix_m, biovar_m, tlag_m);

  int begin_i = 0;
  for (int i = 0; i < np; ++i) {
    MatrixXd x_i(x_m.rows(), len[i]);
    for (int j = 0; j < len[i]; ++j) {
      x_i.col(j) = x_m.col(begin_i + j);
    }
    torsten::test::test_val(x_i, x);
    begin_i += len[i];
  }
}

TEST_F(TorstenPopulationPMXTwoCptTest, multiple_IV_doses_data_only) {
  rate[0] = 300;
  for (int i = 0; i < np; ++i) {
    for (int j = 0; j < nt; ++j) {      
      rate_m[i * nt + j] = rate[j];
    }
  }

  Matrix<double, Dynamic, Dynamic> x =
    torsten::pmx_solve_twocpt(time, amt, rate, ii, evid, cmt, addl, ss, pMatrix, biovar, tlag); // NOLINT

  Matrix<double, Dynamic, Dynamic> x_m =
    torsten::pmx_solve_group_twocpt(len, time_m, amt_m, rate_m, ii_m, evid_m, cmt_m, addl_m, ss_m, // NOLINT
                                    pMatrix_m, biovar_m, tlag_m);

  int begin_i = 0;
  for (int i = 0; i < np; ++i) {
    MatrixXd x_i(x_m.rows(), len[i]);
    for (int j = 0; j < len[i]; ++j) {
      x_i.col(j) = x_m.col(begin_i + j);
    }
    torsten::test::test_val(x_i, x);
    begin_i += len[i];
  }
}

TEST_F(TorstenPopulationPMXTwoCptTest, multiple_bolus_doses_par_var) {
  vector<vector<var> > pMatrix_m_v(np);
  vector<vector<var> > pMatrix_v(torsten::to_var(pMatrix));
  for (int i = 0; i < np; ++i) {
    pMatrix_m_v[i] = stan::math::to_var(pMatrix[0]);
  }

  Matrix<var, Dynamic, Dynamic> x =
    torsten::pmx_solve_twocpt(time, amt, rate, ii, evid, cmt, addl, ss, pMatrix_v, biovar, tlag); // NOLINT

  Matrix<var, Dynamic, Dynamic> x_m =
    torsten::pmx_solve_group_twocpt(len, time_m, amt_m, rate_m, ii_m, evid_m, cmt_m, addl_m, ss_m, // NOLINT
                              pMatrix_m_v, biovar_m, tlag_m);

  int begin_i = 0;
  for (int i = 0; i < np; ++i) {
    Matrix<var, Dynamic, Dynamic> x_i(x_m.rows(), len[i]);
    for (int j = 0; j < len[i]; ++j) {
      x_i.col(j) = x_m.col(begin_i + j);
    }
    torsten::test::test_grad(pMatrix_m_v[i], pMatrix_v[0], x_i, x);
    begin_i += len[i];
  }
}

TEST_F(TorstenPopulationPMXTwoCptTest, multiple_IV_doses_par_var) {
  rate[0] = 300;
  for (int i = 0; i < np; ++i) {
    for (int j = 0; j < nt; ++j) {      
      rate_m[i * nt + j] = rate[j];
    }
  }

  vector<vector<var> > pMatrix_m_v(np);
  vector<vector<var> > pMatrix_v(torsten::to_var(pMatrix));
  for (int i = 0; i < np; ++i) {
    pMatrix_m_v[i] = stan::math::to_var(pMatrix[0]);
  }

  Matrix<var, Dynamic, Dynamic> x =
    torsten::pmx_solve_twocpt(time, amt, rate, ii, evid, cmt, addl, ss, pMatrix_v, biovar, tlag); // NOLINT

  Matrix<var, Dynamic, Dynamic> x_m =
    torsten::pmx_solve_group_twocpt(len, time_m, amt_m, rate_m, ii_m, evid_m, cmt_m, addl_m, ss_m, // NOLINT
                              pMatrix_m_v, biovar_m, tlag_m);

  int begin_i = 0;
  for (int i = 0; i < np; ++i) {
    Matrix<var, Dynamic, Dynamic> x_i(x_m.rows(), len[i]);
    for (int j = 0; j < len[i]; ++j) {
      x_i.col(j) = x_m.col(begin_i + j);
    }
    torsten::test::test_grad(pMatrix_m_v[i], pMatrix_v[0], x_i, x);
    begin_i += len[i];
  }
}
