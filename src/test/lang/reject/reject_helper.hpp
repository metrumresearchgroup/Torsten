#ifndef TEST_UNIT_LANG_REJECT_REJECT_HELPER_HPP
#define TEST_UNIT_LANG_REJECT_REJECT_HELPER_HPP

#include <gtest/gtest.h>
#include <stdexcept>
#include <sstream>
#include <fstream>
#include <boost/random/additive_combine.hpp>
#include <stan/io/dump.hpp>
#include <stan/callbacks/stream_writer.hpp>

void expect_substring(const std::string& msg,
                      const std::string& expected_substring) {
  if (msg.find(expected_substring) == std::string::npos)
    FAIL() << "expected to find substring=" << expected_substring
           << " in string=" << msg << std::endl;
}

template <class M, class E>
void reject_test(const std::string& expected_msg1 = "",
                 const std::string& expected_msg2 = "",
                 const std::string& expected_msg3 = "") {
  std::fstream empty_data_stream("");
  stan::io::dump empty_data_context(empty_data_stream);
  empty_data_stream.close();

  std::stringstream model_output;

  boost::ecuyer1988 base_rng;
  base_rng.seed(123456);

  std::stringstream out;
  try {
    M model(empty_data_context, 0, &model_output);
    std::vector<double> cont_vector(model.num_params_r(), 0.0);
    std::vector<int> disc_vector;
    double lp
        = model.template log_prob<false, false>(cont_vector, disc_vector, &out);
    (void)lp;
    stan::callbacks::stream_writer writer(out);
    std::vector<double> params;
    std::stringstream ss;
    model.write_array(base_rng, cont_vector, disc_vector, params, true, true,
                      &ss);
  } catch (const E& e) {
    expect_substring(e.what(), expected_msg1);
    expect_substring(e.what(), expected_msg2);
    expect_substring(e.what(), expected_msg3);
    EXPECT_EQ("", out.str());
    return;
  }
  EXPECT_EQ("", out.str());
  FAIL() << "model failed to reject" << std::endl;
}

template <class M, class E>
void print_reject_test(const std::string& expected_msg1 = "") {
  std::fstream empty_data_stream("");
  stan::io::dump empty_data_context(empty_data_stream);
  empty_data_stream.close();

  std::stringstream model_output;

  boost::ecuyer1988 base_rng;
  base_rng.seed(123456);

  std::stringstream ss;
  try {
    M model(empty_data_context, 0, &ss);
    std::vector<double> cont_vector(model.num_params_r(), 0.0);
    std::vector<int> disc_vector;
    double lp
        = model.template log_prob<false, false>(cont_vector, disc_vector, &ss);
    (void)lp;
    std::vector<double> params;
    model.write_array(base_rng, cont_vector, disc_vector, params, true, true,
                      &ss);
  } catch (const E& e) {
    expect_substring(ss.str(), expected_msg1);
    return;
  }
  FAIL() << "model failed to reject" << std::endl;
}

#endif
