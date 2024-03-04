#include <stan/services/util/run_adaptive_sampler.hpp>
#include <gtest/gtest.h>
#include <test/test-models/good/services/test_lp.hpp>
#include <stan/io/empty_var_context.hpp>
#include <stan/services/util/create_rng.hpp>
#include <test/unit/services/instrumented_callbacks.hpp>
#include <test/unit/mcmc/hmc/mock_hmc.hpp>
#include <stan/mcmc/hmc/nuts/adapt_unit_e_nuts.hpp>

class ServicesUtil : public testing::Test {
 public:
  ServicesUtil()
      : model(context, 0, &model_log),
        rng(stan::services::util::create_rng(0, 1)),
        sampler(model, rng),
        num_warmup(0),
        num_samples(0),
        num_thin(1),
        refresh(0),
        save_warmup(false) {
    cont_vector.push_back(0);
    cont_vector.push_back(0);
  }

  std::stringstream model_log;
  stan::io::empty_var_context context;
  stan_model model;
  std::vector<double> cont_vector;
  boost::ecuyer1988 rng;
  stan::test::unit::instrumented_interrupt interrupt;
  stan::test::unit::instrumented_writer sample_writer, diagnostic_writer;
  stan::test::unit::instrumented_logger logger;
  stan::mcmc::adapt_unit_e_nuts<stan_model, boost::ecuyer1988> sampler;
  int num_warmup, num_samples, num_thin, refresh;
  bool save_warmup;
};

TEST_F(ServicesUtil, all_zero) {
  stan::services::util::run_adaptive_sampler(
      sampler, model, cont_vector, num_warmup, num_samples, num_thin, refresh,
      save_warmup, rng, interrupt, logger, sample_writer, diagnostic_writer);
  EXPECT_EQ(0, interrupt.call_count());

  EXPECT_EQ(3 + 2, logger.call_count()) << "Writes the elapsed time";
  EXPECT_EQ(logger.call_count(), logger.call_count_info())
      << "No other calls to logger";

  EXPECT_EQ(8, sample_writer.call_count());
  EXPECT_EQ(1, sample_writer.call_count("vector_string")) << "header line";
  EXPECT_EQ(2 + 3, sample_writer.call_count("string"))
      << "adaptation info + elapsed time";
  EXPECT_EQ(2, sample_writer.call_count("empty")) << "blank lines";

  EXPECT_EQ(6, diagnostic_writer.call_count());
  EXPECT_EQ(1, diagnostic_writer.call_count("vector_string")) << "header line";
  EXPECT_EQ(3, diagnostic_writer.call_count("string")) << "elapsed time";
  EXPECT_EQ(2, diagnostic_writer.call_count("empty")) << "blank lines";
}

TEST_F(ServicesUtil, num_warmup_no_save) {
  num_warmup = 1000;
  stan::services::util::run_adaptive_sampler(
      sampler, model, cont_vector, num_warmup, num_samples, num_thin, refresh,
      save_warmup, rng, interrupt, logger, sample_writer, diagnostic_writer);
  EXPECT_EQ(num_warmup, interrupt.call_count());

  EXPECT_EQ(3 + 2, logger.call_count()) << "Writes the elapsed time";
  EXPECT_EQ(logger.call_count(), logger.call_count_info())
      << "No other calls to logger";

  EXPECT_EQ(8, sample_writer.call_count());
  EXPECT_EQ(1, sample_writer.call_count("vector_string")) << "header line";
  EXPECT_EQ(2 + 3, sample_writer.call_count("string"))
      << "adaptation info + elapsed time";
  EXPECT_EQ(2, sample_writer.call_count("empty")) << "blank lines";

  EXPECT_EQ(6, diagnostic_writer.call_count());
  EXPECT_EQ(1, diagnostic_writer.call_count("vector_string")) << "header line";
  EXPECT_EQ(3, diagnostic_writer.call_count("string")) << "elapsed time";
  EXPECT_EQ(2, diagnostic_writer.call_count("empty")) << "blank lines";
}

TEST_F(ServicesUtil, num_warmup_save) {
  num_warmup = 1000;
  save_warmup = true;
  stan::services::util::run_adaptive_sampler(
      sampler, model, cont_vector, num_warmup, num_samples, num_thin, refresh,
      save_warmup, rng, interrupt, logger, sample_writer, diagnostic_writer);
  EXPECT_EQ(num_warmup, interrupt.call_count());

  EXPECT_EQ(3 + 2, logger.call_count()) << "Writes the elapsed time";
  EXPECT_EQ(logger.call_count(), logger.call_count_info())
      << "No other calls to logger";

  EXPECT_EQ(num_warmup + 8, sample_writer.call_count());
  EXPECT_EQ(1, sample_writer.call_count("vector_string")) << "header line";
  EXPECT_EQ(2 + 3, sample_writer.call_count("string"))
      << "adaptation info + elapsed time";
  EXPECT_EQ(2, sample_writer.call_count("empty")) << "blank lines";
  EXPECT_EQ(num_warmup, sample_writer.call_count("vector_double"))
      << "warmup draws";

  EXPECT_EQ(num_warmup + 6, diagnostic_writer.call_count());
  EXPECT_EQ(1, diagnostic_writer.call_count("vector_string")) << "header line";
  EXPECT_EQ(3, diagnostic_writer.call_count("string")) << "elapsed time";
  EXPECT_EQ(2, diagnostic_writer.call_count("empty")) << "blank lines";
  EXPECT_EQ(num_warmup, diagnostic_writer.call_count("vector_double"))
      << "warmup draws";
}

TEST_F(ServicesUtil, num_samples) {
  num_samples = 1000;
  stan::services::util::run_adaptive_sampler(
      sampler, model, cont_vector, num_warmup, num_samples, num_thin, refresh,
      save_warmup, rng, interrupt, logger, sample_writer, diagnostic_writer);
  EXPECT_EQ(num_samples, interrupt.call_count());

  EXPECT_EQ(3 + 2, logger.call_count()) << "Writes the elapsed time";
  EXPECT_EQ(logger.call_count(), logger.call_count_info())
      << "No other calls to logger";

  EXPECT_EQ(num_samples + 8, sample_writer.call_count());
  EXPECT_EQ(1, sample_writer.call_count("vector_string")) << "header line";
  EXPECT_EQ(2 + 3, sample_writer.call_count("string"))
      << "adaptation info + elapsed time";
  EXPECT_EQ(2, sample_writer.call_count("empty")) << "blank lines";
  EXPECT_EQ(num_samples, sample_writer.call_count("vector_double"))
      << "num_samples draws";

  EXPECT_EQ(num_samples + 6, diagnostic_writer.call_count());
  EXPECT_EQ(1, diagnostic_writer.call_count("vector_string")) << "header line";
  EXPECT_EQ(3, diagnostic_writer.call_count("string")) << "elapsed time";
  EXPECT_EQ(2, diagnostic_writer.call_count("empty")) << "blank lines";
  EXPECT_EQ(num_samples, sample_writer.call_count("vector_double"))
      << "num_samples draws";
}

TEST_F(ServicesUtil, num_warmup_save_num_samples_num_thin) {
  num_warmup = 500;
  save_warmup = true;
  num_samples = 500;
  num_thin = 10;
  stan::services::util::run_adaptive_sampler(
      sampler, model, cont_vector, num_warmup, num_samples, num_thin, refresh,
      save_warmup, rng, interrupt, logger, sample_writer, diagnostic_writer);
  EXPECT_EQ(num_warmup + num_samples, interrupt.call_count());

  EXPECT_EQ(3 + 2, logger.call_count()) << "Writes the elapsed time";
  EXPECT_EQ(logger.call_count(), logger.call_count_info())
      << "No other calls to logger";

  EXPECT_EQ((num_warmup + num_samples) / num_thin + 8,
            sample_writer.call_count());
  EXPECT_EQ(1, sample_writer.call_count("vector_string")) << "header line";
  EXPECT_EQ(2 + 3, sample_writer.call_count("string")) << "elapsed time";
  EXPECT_EQ(2, sample_writer.call_count("empty")) << "blank lines";
  EXPECT_EQ((num_warmup + num_samples) / num_thin,
            sample_writer.call_count("vector_double"))
      << "thinned warmup and draws";

  EXPECT_EQ((num_warmup + num_samples) / num_thin + 6,
            diagnostic_writer.call_count());
  EXPECT_EQ(1, diagnostic_writer.call_count("vector_string")) << "header line";
  EXPECT_EQ(3, diagnostic_writer.call_count("string")) << "elapsed time";
  EXPECT_EQ(2, diagnostic_writer.call_count("empty")) << "blank lines";
  EXPECT_EQ((num_warmup + num_samples) / num_thin,
            diagnostic_writer.call_count("vector_double"))
      << "thinned warmup and draws";
}

TEST_F(ServicesUtil, num_warmup_num_samples_refresh) {
  num_warmup = 500;
  num_samples = 500;
  refresh = 10;
  stan::services::util::run_adaptive_sampler(
      sampler, model, cont_vector, num_warmup, num_samples, num_thin, refresh,
      save_warmup, rng, interrupt, logger, sample_writer, diagnostic_writer);
  EXPECT_EQ(num_warmup + num_samples, interrupt.call_count());

  EXPECT_EQ((num_warmup + num_samples) / refresh + 2 + 3 + 2,
            logger.call_count())
      << "Writes 1 to start warmup, 1 to start post-warmup, and "
      << "(num_warmup + num_samples) / refresh, then the elapsed time";
  EXPECT_EQ(logger.call_count(), logger.call_count_info())
      << "No other calls to logger";

  EXPECT_EQ(num_samples + 8, sample_writer.call_count());
  EXPECT_EQ(1, sample_writer.call_count("vector_string")) << "header line";
  EXPECT_EQ(2 + 3, sample_writer.call_count("string")) << "elapsed time";
  EXPECT_EQ(2, sample_writer.call_count("empty")) << "blank lines";
  EXPECT_EQ(num_samples, sample_writer.call_count("vector_double")) << "draws";

  EXPECT_EQ(num_samples + 6, diagnostic_writer.call_count());
  EXPECT_EQ(1, diagnostic_writer.call_count("vector_string")) << "header line";
  EXPECT_EQ(3, diagnostic_writer.call_count("string")) << "elapsed time";
  EXPECT_EQ(2, diagnostic_writer.call_count("empty")) << "blank lines";
  EXPECT_EQ(num_samples, diagnostic_writer.call_count("vector_double"))
      << "draws";
}
