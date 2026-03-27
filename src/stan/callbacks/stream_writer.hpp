#ifndef STAN_CALLBACKS_STREAM_WRITER_HPP
#define STAN_CALLBACKS_STREAM_WRITER_HPP

#include <stan/callbacks/writer.hpp>
#include <ostream>
#include <vector>
#include <string>

namespace stan {
namespace callbacks {

/**
 * <code>stream_writer</code> is an implementation
 * of <code>writer</code> that writes to a stream.
 */
class stream_writer : public writer {
 public:
  /**
   * Constructs a stream writer with an output stream
   * and an optional prefix for comments.
   *
   * @param[in, out] output stream to write
   * @param[in] comment_prefix string to stream before
   *   each comment line. Default is "".
   */
  explicit stream_writer(std::ostream& output,
                         const std::string& comment_prefix = "")
      : output_(output), comment_prefix_(comment_prefix) {}

  /**
   * Virtual destructor
   */
  virtual ~stream_writer() {}

  /**
   * Writes a set of names on a single line in csv format followed
   * by a newline.
   *
   * Note: the names are not escaped.
   *
   * @param[in] names Names in a std::vector
   */
  void operator()(const std::vector<std::string>& names) {
    write_vector(names);
  }

  /**
   * Writes a set of values in csv format followed by a newline.
   *
   * Note: the precision of the output is determined by the settings
   *  of the stream on construction.
   *
   * @param[in] state Values in a std::vector
   */
  void operator()(const std::vector<double>& state) { write_vector(state); }

  /**
   * Writes the comment_prefix to the stream followed by a newline.
   */
  void operator()() { output_ << comment_prefix_ << std::endl; }

  /**
   * Writes the comment_prefix then the message followed by a newline.
   *
   * @param[in] message A string
   */
  void operator()(const std::string& message) {
    output_ << comment_prefix_ << message << std::endl;
  }

  /**
   * Writes multiple rows and columns of values in csv format.
   *
   * Note: the precision of the output is determined by the settings
   *  of the stream on construction.
   *
   * @param[in] values A matrix of values. The input is expected to have
   * parameters in the rows and samples in the columns. The matrix is then
   * transposed for the output.
   */
  void operator()(const Eigen::Matrix<double, -1, -1>& values) {
    output_ << values.transpose().format(CommaInitFmt);
  }
  /**
   * Write a row of values in csv format.
   *
   * Note: the precision of the output is determined by the settings
   *  of the stream on construction.
   *
   * @param[in] values A column vector of values.
   */
  void operator()(const Eigen::Matrix<double, -1, 1>& values) {
    output_ << values.transpose().format(CommaInitFmt);
  }

  /**
   * Write a row of values in csv format
   *
   * Note: the precision of the output is determined by the settings
   *  of the stream on construction.
   *
   * @param[in] values A row vector of values.
   */
  void operator()(const Eigen::Matrix<double, 1, -1>& values) {
    output_ << values.format(CommaInitFmt);
  }

  /**
   * Checks if stream is valid.
   */
  virtual bool is_valid() const noexcept { return output_.good(); }

 private:
  /**
   * Comma formatter for writing Eigen matrices
   */
  Eigen::IOFormat CommaInitFmt{
      Eigen::StreamPrecision, Eigen::DontAlignCols, ", ", "", "", "\n", "", ""};

  /**
   * Output stream
   */
  std::ostream& output_;

  /**
   * Comment prefix to use when printing comments: strings and blank lines
   */
  std::string comment_prefix_;

  /**
   * Writes a set of values in csv format followed by a newline.
   *
   * Note: the precision of the output is determined by the settings
   *  of the stream on construction.
   *
   * @param[in] v Values in a std::vector
   */
  template <class T>
  void write_vector(const std::vector<T>& v) {
    if (v.empty())
      return;

    typename std::vector<T>::const_iterator last = v.end();
    --last;

    for (typename std::vector<T>::const_iterator it = v.begin(); it != last;
         ++it)
      output_ << *it << ",";
    output_ << v.back() << std::endl;
  }
};

}  // namespace callbacks
}  // namespace stan
#endif
