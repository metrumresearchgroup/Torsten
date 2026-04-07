#ifndef STAN_CALLBACKS_CONCURRENT_WRITER_HPP
#define STAN_CALLBACKS_CONCURRENT_WRITER_HPP

#include <stan/math/prim/fun/Eigen.hpp>
#include <stan/math/prim/meta.hpp>
#include <tbb/concurrent_queue.h>
#include <condition_variable>
#include <functional>
#include <string>
#include <thread>
#include <vector>

namespace stan::callbacks {
#ifdef STAN_THREADS
/**
 * Enables thread-safe writing of numeric values to a writer.
 * On construction, a thread is spawned to write to the writer.
 * This class uses an `std::thread` instead of a tbb task graph because
 * of deadlocking issues. A deadlock can happen in two major cases.
 * 1. If TBB gives all threads a task, and all threads hit an instance of max
 * capacity. TBB can choose to wait for a thread to finish instead of spinning
 * up the write thread. So to circumvent that issue, we use an std::thread.
 * 2. If the bounded queues are full but the consumer thread is not scheduled
 *  because there are more busy threads than the number of threads available.
 * The producer threads are blocked because the queues are full. The consumer
 * thread is blocked because the producer thread is spinning. Then we have a
 * deadlock because the consumer thread is blocked because the producer threads
 * are blocked.
 * i.e. queue(full)->producer(blocked)->consumer(blocked)->producer(blocked)
 * To circumvent this issue, we check in the producer threads if
 * the queues are almost* full and if they are, we make a lock and wait for
 * the consumer thread to signal it's queues are not longer at capacity. This
 * frees a thread for the consumer thread to write to the writer. Once the
 * consumer thread is finished writing, it will notify all the producer threads
 * to continue sending data. The check for the queues being almost full is
 * done by checking if the size of the queue is greater than the max capacity
 * minus the number of threads times 2. This is a heuristic to make sure that,
 * if some threads slip past the check and write to the queue, the bounded queue
 * will still not be full.
 *
 * @tparam Writer A type that inherits from `writer`
 */
template <typename Writer>
struct concurrent_writer {
  // A reference to the writer to write to
  std::reference_wrapper<Writer> writer;
  // Queue for Eigen vector messages
  tbb::concurrent_bounded_queue<Eigen::RowVectorXd> eigen_messages_{};
  // Block threads from writing to queues if the queues are full
  std::mutex block_{};
  // The writing thread
  std::thread thread_;
  // Condition variable to signal the writing thread to continue
  std::condition_variable cv;
  // Maximum number of threads that can be in use
  std::size_t max_threads{tbb::global_control::max_allowed_parallelism};
  // Max capacity of queue
  std::size_t max_capacity{1000 + max_threads};
  // Threshold where the writing threads will wait for the queues to empty
  std::size_t wait_threshold{max_capacity - max_threads - 1};
  // Flag to stop the writing thread once all queues are empty
  bool continue_writing_{true};

  /**
   * Constructs a concurrent writer from a writer.
   * @note This will start a thread to write to the writer.
   * @param writer A writer to write to
   */
  explicit concurrent_writer(Writer& writer) : writer(writer) {
    eigen_messages_.set_capacity(max_capacity);
    thread_ = std::thread([&]() {
      Eigen::RowVectorXd eigen;
      while (continue_writing_ || !eigen_messages_.empty()) {
        while (eigen_messages_.try_pop(eigen)) {
          writer(eigen);
        }
        if (this->empty()) {
          cv.notify_all();
          std::this_thread::yield();
        }
      }
    });
  }

  /**
   * Checks if all queues are empty
   */
  inline bool empty() { return eigen_messages_.empty(); }

  /**
   * Check if any of the queues are at capacity
   */
  inline bool hit_capacity() {
    return eigen_messages_.size() >= wait_threshold;
  }

  /**
   * Place a value in a queue for writing.
   * @note If any of the queues are at capacity, the thread yields itself until
   * the the queues empty. In the case of spurious startups the wait just checks
   *  that the queues are not full.
   * @tparam T An Eigen vector
   * @param t A value to put on a queue
   */
  template <typename T>
  void operator()(T&& t) {
    bool pushed = false;
    if (this->hit_capacity()) {
      std::unique_lock lk(block_);
      cv.wait(lk, [this_ = this] { return !(this_->hit_capacity()); });
    }
    while (!pushed) {
      if constexpr (stan::is_std_vector<T>::value) {
        pushed = eigen_messages_.try_push(
            Eigen::RowVectorXd::Map(t.data(), t.size()));
      } else if constexpr (stan::is_eigen_vector<T>::value) {
        pushed = eigen_messages_.try_push(std::forward<T>(t));
      } else {
        constexpr bool is_numeric_std_vector
            = stan::is_std_vector<T>::value
              && std::is_arithmetic_v<stan::value_type_t<T>>;
        static_assert(
            (!is_numeric_std_vector && !stan::is_eigen_vector<T>::value),
            "Unsupported type passed to concurrent_writer. This is an "
            "internal error. Please file an issue on the stan github "
            "repository with the error log from the compiler.\n"
            "https://github.com/stan-dev/stan/issues/new?template=Blank+issue");
      }
      if (!pushed) {
        std::this_thread::yield();
      }
    }
  }

  /**
   * Waits till all writes are finished on the thread
   */
  void wait() {
    continue_writing_ = false;
    if (thread_.joinable()) {
      // If any threads are waiting for the queues to empty, notify them
      cv.notify_all();
      thread_.join();
    }
  }
  /**
   * Destructor makes sure the thread is joined before destruction
   */
  ~concurrent_writer() { wait(); }
};
#else
/**
 * When STAN_THREADS is not defined, the concurrent writer is just a wrapper
 */
template <typename Writer>
struct concurrent_writer {
  std::reference_wrapper<Writer> writer;
  explicit concurrent_writer(Writer& writer) : writer(writer) {}
  template <typename T>
  void operator()(T&& t) {
    writer(std::forward<T>(t));
  }
  inline static constexpr void wait() {}
};
#endif
}  // namespace stan::callbacks
#endif
