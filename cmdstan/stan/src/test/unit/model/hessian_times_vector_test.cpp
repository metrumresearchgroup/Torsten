#include <stan/model/hessian_times_vector.hpp>
#include <test/test-models/good/model/valid.hpp>
#include <gtest/gtest.h>

TEST(ModelUtil, hessian_times_vector) {
  int dim = 5;

  Eigen::VectorXd x(dim);
  Eigen::VectorXd v(dim);
  double f;
  Eigen::VectorXd hess_f_dot_v(dim);

  std::fstream data_stream(std::string("").c_str(), std::fstream::in);
  stan::io::dump data_var_context(data_stream);
  data_stream.close();

  std::stringstream output;
  valid_model_namespace::valid_model valid_model(data_var_context, 0, &output);
  EXPECT_NO_THROW(
      stan::model::hessian_times_vector(valid_model, x, v, f, hess_f_dot_v));

  EXPECT_FLOAT_EQ(dim, x.size());
  EXPECT_FLOAT_EQ(dim, v.size());
  EXPECT_FLOAT_EQ(dim, hess_f_dot_v.size());

  EXPECT_EQ("", output.str());

  // Incorporate once operands and partials has been generalized
  // output.str("");
  // domain_fail_namespace::domain_fail domain_fail_model(data_var_context,
  // &output); EXPECT_THROW(stan::model::hessian_times_vector(domain_fail_model,
  // x, v, f, hess_f_dot_v),
  //             std::domain_error);
  // EXPECT_EQ("", output.str());
}
