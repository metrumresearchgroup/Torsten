#ifndef STAN_IO_CHAINED_VAR_CONTEXT_HPP
#define STAN_IO_CHAINED_VAR_CONTEXT_HPP

#include <stan/io/var_context.hpp>
#include <stan/io/validate_dims.hpp>
#include <string>
#include <vector>

namespace stan {
namespace io {

/**
 * A chained_var_context object represents two objects of var_context
 * as one.
 */
class chained_var_context : public var_context {
 private:
  const var_context& vc1_;
  const var_context& vc2_;

 public:
  chained_var_context(const var_context& v1, const var_context& v2)
      : vc1_(v1), vc2_(v2) {}

  bool contains_i(const std::string& name) const {
    return vc1_.contains_i(name) || vc2_.contains_i(name);
  }

  bool contains_r(const std::string& name) const {
    return vc1_.contains_r(name) || vc2_.contains_r(name);
  }

  std::vector<double> vals_r(const std::string& name) const {
    return vc1_.contains_r(name) ? vc1_.vals_r(name) : vc2_.vals_r(name);
  }

  std::vector<int> vals_i(const std::string& name) const {
    return vc1_.contains_i(name) ? vc1_.vals_i(name) : vc2_.vals_i(name);
  }

  std::vector<size_t> dims_r(const std::string& name) const {
    return vc1_.contains_r(name) ? vc1_.dims_r(name) : vc2_.dims_r(name);
  }

  std::vector<size_t> dims_i(const std::string& name) const {
    return vc1_.contains_r(name) ? vc1_.dims_i(name) : vc2_.dims_i(name);
  }

  void names_r(std::vector<std::string>& names) const {
    vc1_.names_r(names);
    std::vector<std::string> names2;
    vc2_.names_r(names2);
    names.insert(names.end(), names2.begin(), names2.end());
  }

  void names_i(std::vector<std::string>& names) const {
    vc1_.names_i(names);
    std::vector<std::string> names2;
    vc2_.names_i(names2);
    names.insert(names.end(), names2.begin(), names2.end());
  }

  /**
   * Check variable dimensions against variable declaration.
   * Only used for data read in from file.
   *
   * @param stage stan program processing stage
   * @param name variable name
   * @param base_type declared stan variable type
   * @param dims variable dimensions
   * @throw std::runtime_error if mismatch between declared
   *        dimensions and dimensions found in context.
   */
  void validate_dims(const std::string& stage, const std::string& name,
                     const std::string& base_type,
                     const std::vector<size_t>& dims_declared) const {
    stan::io::validate_dims(*this, stage, name, base_type, dims_declared);
  }
};
}  // namespace io
}  // namespace stan

#endif
