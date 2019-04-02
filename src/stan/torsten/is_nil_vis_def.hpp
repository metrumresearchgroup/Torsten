#ifndef STAN_LANG_TORSTEN_AST_FUN_IS_NIL_VIS_DEF_HPP
#define STAN_LANG_TORSTEN_AST_FUN_IS_NIL_VIS_DEF_HPP

namespace stan {
  namespace lang {

    bool is_nil_vis::operator()(const univariate_integral_control& /* x */) const {
      return false;
    }
    bool is_nil_vis::operator()(const generalOdeModel_control& /* x */) const {
      return false;
    }
    bool is_nil_vis::operator()(const generalOdeModel& /* x */) const {
      return false;
    }
    bool is_nil_vis::operator()(const pmx_solve_group& /* x */) const {
      return false;
    }
  }
}

#endif
