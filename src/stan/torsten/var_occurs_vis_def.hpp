#ifndef STAN_LANG_TORSTEN_AST_FUN_VAR_OCCURS_VIS_DEF_HPP
#define STAN_LANG_TORSTEN_AST_FUN_VAR_OCCURS_VIS_DEF_HPP

#include <stan/lang/ast.hpp>

namespace stan {
  namespace lang {

    bool var_occurs_vis::operator()(const univariate_integral_control& e) const {
      return false;  // no refs persist out of univariate_integral_control() call
    }
    bool var_occurs_vis::operator()(const generalOdeModel_control& e) const {
      return false;  // no refs persist out of generalOdeModel_control() call
    }
    bool var_occurs_vis::operator()(const generalOdeModel& e) const {
      return false;  // no refs persist out of generalOdeModel() call
    }

  }
}

#endif
