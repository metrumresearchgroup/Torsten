#ifndef STAN_MCMC_HMC_MPI_CROSS_CHAIN_ADAPTER_HPP
#define STAN_MCMC_HMC_MPI_CROSS_CHAIN_ADAPTER_HPP

#include <stan/callbacks/writer.hpp>
#include <stan/callbacks/interrupt.hpp>
#include <stan/mcmc/hmc/base_hmc.hpp>
#include <stan/mcmc/stepsize_covar_adapter.hpp>
#include <stan/mcmc/cross_chain/mpi_var_adaptation.hpp>
#include <stan/mcmc/cross_chain/mpi_covar_adaptation.hpp>
#include <stan/services/util/mcmc_writer.hpp>
#include <stan/analyze/mcmc/autocovariance.hpp>
#include <stan/analyze/mcmc/compute_effective_sample_size.hpp>
#include <boost/accumulators/accumulators.hpp>
#include <boost/accumulators/statistics/stats.hpp>
#include <boost/accumulators/statistics/mean.hpp>
#include <boost/accumulators/statistics/variance.hpp>
#include <boost/accumulators/statistics/count.hpp>
#include <string>
#include <chrono>

#ifdef MPI_ADAPTED_WARMUP
#include <stan/math/torsten/mpi/session.hpp>
#endif

namespace stan {
namespace mcmc {

#ifdef MPI_ADAPTED_WARMUP
  class mpi_cross_chain_adapter {
  protected:
    bool is_adapted_;
    int init_buffer_;
    int term_buffer_;
    int term_buffer_counter_;
    int window_size_; 
    int num_chains_;
    int max_num_windows_;
    double target_rhat_;
    double target_ess_;
    std::vector<double> lp_draws_;
    Eigen::MatrixXd all_lp_draws_;
    std::vector<boost::accumulators::accumulator_set<double,
                                                     boost::accumulators::stats<boost::accumulators::tag::mean, // NOLINT
                                                                                boost::accumulators::tag::variance>>> lp_acc_; // NOLINT
    boost::accumulators::accumulator_set<int,
                                         boost::accumulators::features<boost::accumulators::tag::count> > draw_count_acc_;
    Eigen::ArrayXd rhat_;
    Eigen::ArrayXd ess_;
    mpi_metric_adaptation* var_adapt;
  public:
    const static int nd_win = 2;

    mpi_cross_chain_adapter() = default;

    mpi_cross_chain_adapter(int init_buffer, int term_buffer,
                            int num_iterations, int window_size,
                            int num_chains,
                            double target_rhat, double target_ess) :
      is_adapted_(false),
      init_buffer_(init_buffer),
      term_buffer_(term_buffer),
      term_buffer_counter_(0),
      window_size_(window_size),
      num_chains_(num_chains),
      max_num_windows_(num_iterations / window_size),
      target_rhat_(target_rhat),
      target_ess_(target_ess),
      lp_draws_(window_size),
      all_lp_draws_(),
      lp_acc_(max_num_windows_),
      draw_count_acc_(),
      rhat_(Eigen::ArrayXd::Zero(max_num_windows_)),
      ess_(Eigen::ArrayXd::Zero(max_num_windows_))
    {}

    inline void set_cross_chain_metric_adaptation(mpi_metric_adaptation* ptr) {var_adapt = ptr;}

    inline void set_cross_chain_adaptation_params(int init_buffer, int term_buffer,
                                                  int num_iterations,
                                                  int window_size,
                                                  int num_chains,
                                                  double target_rhat, double target_ess) {
      is_adapted_ = false;
      init_buffer_ = init_buffer;
      term_buffer_ = term_buffer;
      term_buffer_counter_ = 0,
      window_size_ = window_size;
      num_chains_ = num_chains;
      max_num_windows_ = num_iterations / window_size;
      target_rhat_ = target_rhat;
      target_ess_ = target_ess;
      lp_draws_.resize(window_size);
      lp_acc_.clear();
      lp_acc_.resize(max_num_windows_);
      draw_count_acc_ = {};
      rhat_ = Eigen::ArrayXd::Zero(max_num_windows_);
      ess_ = Eigen::ArrayXd::Zero(max_num_windows_);
    }

    inline int max_num_windows() {return max_num_windows_;}

    /*
     * current number of draws.
     */
    inline int num_cross_chain_draws() {
      return boost::accumulators::count(draw_count_acc_);
    }

    inline void
    write_num_cross_chain_warmup(callbacks::writer& sample_writer,
                                 int num_thin, int num_warmup) {
      if (use_cross_chain_adapt()) {
        size_t n = num_cross_chain_draws();
        sample_writer("num_warmup = " + std::to_string(n / num_thin));
      } else {
        sample_writer("num_warmup = " + std::to_string(num_warmup));
      }
    }

    /*
     * The number of active windows, also equals to the window
     * to which the most recent
     * sample belongs. That is, if window_size = 10, then the 10th
     * draw belongs to window 1, but 11th draw belongs to window 2.
     */
    inline int num_active_cross_chain_windows() {
      size_t n = num_cross_chain_draws() - 1;
      return n / window_size_ + 1;
    }

    /** 
     * Add samples for current chain. Only add samples to inter-chain ranks
     * 
     * @param s lp___
     * @param q current sample vector
     */
    inline void add_cross_chain_sample(double lp, const Eigen::VectorXd& q) {
      using stan::math::mpi::Session;

      int i = num_cross_chain_draws() % window_size_;
      draw_count_acc_(0);     // move counter
      int n_win = num_active_cross_chain_windows();

      // For term buffer we only need to add dummy sample by move the counter
      if (is_adapted_) {
        term_buffer_counter_++;
      } else if (Session::is_in_inter_chain_comm(num_chains_)) {
        lp_draws_[i] = lp;
        for (int win = 0; win < n_win; ++win) {
          lp_acc_[win](lp);
        }

        // add current samples to var/covar adapter after init buffer
        if (num_cross_chain_draws() > init_buffer_) {
          var_adapt -> add_sample(q, n_win);          
        }
      }
    }

    /*
     * Computes the effective sample size (ESS) for the specified
     * parameter across all kept samples.  The value returned is the
     * minimum of ESS and the number_total_draws *
     * log10(number_total_draws).
     */
    inline double compute_effective_sample_size(int win, int win_count,
                                                const Eigen::MatrixXd& all_lp_draws) {
      std::vector<const double*> draws(num_chains_);
      size_t num_draws = (win_count - win) * window_size_;
      for (int chain = 0; chain < num_chains_; ++chain) {
        draws[chain] = &all_lp_draws(win * window_size_, chain);
      }
      return stan::analyze::compute_effective_sample_size(draws, num_draws);
    }

    inline const Eigen::ArrayXd& cross_chain_adapt_rhat() {
      return rhat_;
    }

    inline const Eigen::ArrayXd& cross_chain_adapt_ess() {
      return ess_;
    }

    inline bool is_cross_chain_adapt_window_end() {
      size_t n = num_cross_chain_draws();
      return n > 0 && (n % window_size_ == 0);
    }

    inline bool is_cross_chain_adapted() { return is_adapted_; }

    inline void set_cross_chain_adapted(bool a) { is_adapted_ = a; }

    /*
     * if transition needs to be ended in @c generate_transitions()
     * Once warmup converges and term buffer is finished we
     * further increase counter by 1 so this function never
     * returns true.
     */
    inline bool end_transitions() {
      if (is_adapted_ && (term_buffer_counter_ == term_buffer_)) {
        term_buffer_counter_++;
        return true;
      } else {
        return false;
      }
    }

    inline void msg_adaptation(int win, callbacks::logger& logger) {
      std::stringstream message;
      message << "iteration: ";
      message << std::setw(3);
      message << num_cross_chain_draws();
      message << " window: " << win + 1 << " / " << num_active_cross_chain_windows();
      message << std::setw(7) << std::setprecision(4);
      message << " Rhat: " << std::fixed << cross_chain_adapt_rhat()[win];
      const Eigen::ArrayXd& ess(cross_chain_adapt_ess());
      message << " ESS: " << std::fixed << ess_[win];

      logger.info(message);
    }

    /**
     * Gather data from chains. Get from each chain:
     * - mean for window 1, 1+2, 1+2+3, ...
     * - variance for window 1, 1+2, 1+2+3, ...
     * - lp__ draws in the latest window
     * 
     * @param all_chain_gather vector to be filled with gathered data in rank 0.
     * @param all_lp_draws lp__ draws from all the chains, filled by rank 0
     * 
     * @return number of gathered data from each chain. For
     *         non-inter-chain rank this is 0.
     */
    inline int cross_chain_gather(std::vector<double>& all_chain_gather,
                                  Eigen::MatrixXd& all_lp_draws) {
      using boost::accumulators::tag::mean;
      using boost::accumulators::tag::variance;

      using stan::math::mpi::Session;
      using stan::math::mpi::Communicator;

      int n_gather = 0;

      /// only update adaptation at the end of each window
      if ((!is_adapted_) && is_cross_chain_adapt_window_end()) {
        if (Session::is_in_inter_chain_comm(num_chains_)) {
          const Communicator& comm = Session::inter_chain_comm(num_chains_);
          const int win_count = num_active_cross_chain_windows();

          n_gather = nd_win * win_count + window_size_;

          /**
           * prepare data to be gathered in each chain
           */
          std::vector<double> chain_gather(n_gather, 0.0);
          for (int win = 0; win < win_count; ++win) {
            int num_draws = (win_count - win) * window_size_;
            double unbiased_var_scale = num_draws / (num_draws - 1.0);
            chain_gather[nd_win * win] = boost::accumulators::mean(lp_acc_[win]);
            chain_gather[nd_win * win + 1] = boost::accumulators::variance(lp_acc_[win]) *
              unbiased_var_scale;
          }
          std::copy(lp_draws_.begin(), lp_draws_.end(),
                    chain_gather.begin() + nd_win * win_count);

          if (comm.rank() == 0) {
            /**
             * rank 0 gather data from all the chains
             */
            all_chain_gather.resize(n_gather * num_chains_);
            MPI_Gather(chain_gather.data(), n_gather, MPI_DOUBLE,
                       all_chain_gather.data(), n_gather, MPI_DOUBLE, 0, comm.comm());
            all_lp_draws.resize(window_size_ * max_num_windows_, num_chains_);
            int begin_row = (win_count - 1) * window_size_;
            for (int chain = 0; chain < num_chains_; ++chain) {
              int j = n_gather * chain + nd_win * win_count;
              for (int i = 0; i < window_size_; ++i) {
                all_lp_draws(begin_row + i, chain) = all_chain_gather[j + i];
              }
            }
          } else {
            MPI_Gather(chain_gather.data(), n_gather, MPI_DOUBLE,
                       NULL, 0, MPI_DOUBLE, 0, comm.comm());
          }
        }
      }
      return n_gather;
    }

    /** 
     * new cross-chain stepsize is the harmonic mean of chain stepsize
     * 
     * @param chain_stepsize stepsize of current chain
     * 
     * @return new stepsize is the harmonic mean of chain stepsizes
     */
    inline double cross_chain_stepsize(double chain_stepsize) {
      using stan::math::mpi::Session;
      using stan::math::mpi::Communicator;

      double new_stepsize = chain_stepsize;
      if(is_cross_chain_adapt_window_end()) {
        const Communicator& comm = Session::inter_chain_comm(num_chains_);
        if (Session::is_in_inter_chain_comm(num_chains_)) {
          chain_stepsize = 1.0/chain_stepsize;
          MPI_Allreduce(&chain_stepsize, &new_stepsize, 1, MPI_DOUBLE, MPI_SUM, comm.comm());
          new_stepsize = num_chains_ / new_stepsize;
        }
        const Communicator& intra_comm = Session::intra_chain_comm(num_chains_);
        MPI_Bcast(&new_stepsize, 1, MPI_DOUBLE, 0, intra_comm.comm());
      }
      return new_stepsize;
    }

    /** 
     * find adapted window that has maximum ESS and also meets rhat &
     * ess target
     * 
     * @param all_chain_gather all chain information used to calculate
     * rhat & ess. This calculation only occurs in rank 0.
     * 
     * @return adapted window id. A nonnegative id indicates actual
     * widow that meets convergence condition and has the max ESS.
     * A negative id indicates non-convergence, and abs(id) - 1 is the
     * actual id that has the maximum ESS.
     */
    inline int cross_chain_adapted_window(int n_gather,
                                          const std::vector<double>& all_chain_gather,
                                          const Eigen::MatrixXd& all_lp_draws,
                                          callbacks::logger& logger) {
      using boost::accumulators::accumulator_set;
      using boost::accumulators::stats;
      using boost::accumulators::tag::mean;
      using boost::accumulators::tag::variance;
      using math::mpi::Communicator;
      using math::mpi::Session;

      const Communicator& comm = Session::inter_chain_comm(num_chains_);
      const int win_count = num_active_cross_chain_windows();
      int adapted_win = -999;
      if (comm.rank() == 0) {
        rhat_.setZero();
        ess_.setZero();
        double max_ess = 0.0;

        /**
         * loop through windows to check convergence
         */
        for (int win = 0; win < win_count; ++win) {
          bool win_adapted;
          accumulator_set<double, stats<variance>> acc_chain_mean;
          accumulator_set<double, stats<mean>> acc_chain_var;
          accumulator_set<double, stats<mean>> acc_step;
          Eigen::VectorXd chain_mean(num_chains_);
          Eigen::VectorXd chain_var(num_chains_);
          for (int chain = 0; chain < num_chains_; ++chain) {
            chain_mean(chain) = all_chain_gather[chain * n_gather + nd_win * win];
            acc_chain_mean(chain_mean(chain));
            chain_var(chain) = all_chain_gather[chain * n_gather + nd_win * win + 1];
            acc_chain_var(chain_var(chain));
          }
          size_t num_draws = (win_count - win) * window_size_;
          double var_between = num_draws * boost::accumulators::variance(acc_chain_mean)
            * num_chains_ / (num_chains_ - 1);
          double var_within = boost::accumulators::mean(acc_chain_var);
          rhat_(win) = sqrt((var_between / var_within + num_draws - 1) / num_draws);
          ess_[win] = compute_effective_sample_size(win, win_count, all_lp_draws);
          win_adapted = rhat_(win) < target_rhat_ && ess_[win] > target_ess_;

          msg_adaptation(win, logger);

          /// get the win with the largest ESS
          if(ess_[win] > max_ess) {
            max_ess = ess_[win];
            adapted_win = -(win + 1);
            if (win_adapted) {
              adapted_win = std::abs(adapted_win) - 1;
            }
          }
        }
      }
      MPI_Bcast(&adapted_win, 1, MPI_INT, 0, comm.comm());      
      return adapted_win;
    }

    /*
     * @tparam T_metric metric type, <code>Eigen::VectorXd</code> or <code>Eigen::MatrixXd</code>
     * @param[in] m_win number of windows
     * @param[in] window_size window size
     * @param[in] num_chains number of chains
     * @param[in,out] chain_gather gathered information from each chain,
     *                must have enough capacity to store up to
     *                maximum windows for all chains.
     # @return vector {stepsize, rhat(only in rank 0)}
    */
    template<typename T_metric>
    inline bool cross_chain_adaptation(T_metric& inv_e_metric,
                                       callbacks::logger& logger) {
      using boost::accumulators::accumulator_set;
      using boost::accumulators::stats;
      using boost::accumulators::tag::mean;
      using boost::accumulators::tag::variance;

      using stan::math::mpi::Session;
      using stan::math::mpi::Communicator;

      auto t0 = std::chrono::high_resolution_clock::now();
      bool update = false;

      std::vector<double> all_chain_gather;
      int n_gather = cross_chain_gather(all_chain_gather, all_lp_draws_);

      if ((!is_adapted_) && is_cross_chain_adapt_window_end()) {
        const Communicator& comm = Session::inter_chain_comm(num_chains_);
        if (Session::is_in_inter_chain_comm(num_chains_)) {
          const int win_count = num_active_cross_chain_windows();

          // test convergence and return adapt window if any.
          int adapted_win = cross_chain_adapted_window(n_gather,
                                                       all_chain_gather,
                                                       all_lp_draws_,
                                                       logger);
          set_cross_chain_adapted((adapted_win >= 0));

          /// learn metric based on the window with max ESS
          int max_ess_win = is_adapted_ ? adapted_win : (-adapted_win - 1);
          var_adapt -> learn_metric(inv_e_metric, max_ess_win, win_count, comm);
        }

        /// send info to intra-chain nodes
        const Communicator& intra_comm = Session::intra_chain_comm(num_chains_);
        MPI_Bcast(&is_adapted_, 1, MPI_C_BOOL, 0, intra_comm.comm());
        MPI_Bcast(inv_e_metric.data(), inv_e_metric.size(), MPI_DOUBLE, 0, intra_comm.comm());

        if (is_adapted_) {
          // set_cross_chain_stepsize();
        }
        update = true;

        // print aggregating time cost
        if (comm.rank() == 0) {
          auto t1 = std::chrono::high_resolution_clock::now();
          std::stringstream message;
          message << "cross-chain adaptation time: ";
          message << std::chrono::duration_cast<std::chrono::seconds>(t1 - t0).count();
          message << " seconds";
          logger.info(message);
        }
      }
      return update;
    }

    // inline void cross_chain_term_buffer() {
    //   if (is_adapted_) {
    //     term_buffer_counter_++;
    //   }
    // }

    inline bool use_cross_chain_adapt() {
      return num_chains_ > 1;
    }
  };

#else  // sequential version
  class mpi_cross_chain_adapter {
  public:
    inline void set_cross_chain_metric_adaptation(mpi_metric_adaptation* ptr) {}

    inline void set_cross_chain_adaptation_params(int init_buffer, int term_buffer,
                                                  int num_iterations,
                                                  int window_size,
                                                  int num_chains,
                                                  double target_rhat, double target_ess) {}
    inline void add_cross_chain_sample(double lp, const Eigen::VectorXd& q) {}

    template<typename T_metric>
    inline bool cross_chain_adaptation(T_metric& inv_e_metric,
                                       callbacks::logger& logger) { return false; }

    inline bool is_cross_chain_adapted() { return false; }

    inline double cross_chain_stepsize(double chain_stepsize) { return 0.0; }

    inline bool use_cross_chain_adapt() { return false; }
  };
  
#endif

}
}
#endif


