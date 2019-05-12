#ifndef STAN_LANG_TORSTEN_AST_FUN_HAS_NON_PARAM_VAR_VIS_DEF_HPP
#define STAN_LANG_TORSTEN_AST_FUN_HAS_NON_PARAM_VAR_VIS_DEF_HPP

#include <stan/lang/ast.hpp>
#include <boost/variant/apply_visitor.hpp>

namespace stan {
  namespace lang {

    bool has_non_param_var_vis::operator()(const univariate_integral_control& e)
      const {
      return boost::apply_visitor(*this, e.t0_.expr_)
        || boost::apply_visitor(*this, e.t1_.expr_)
        || boost::apply_visitor(*this, e.theta_.expr_);
    }

    bool has_non_param_var_vis::operator()(const generalOdeModel_control& e)
      const {
      // CHECK - anything to do with nonlinearity ?
      // Putting in the variables that may contain var types
      return ((((((boost::apply_visitor(*this, e.time_.expr_)
        || boost::apply_visitor(*this, e.amt_.expr_))
        || boost::apply_visitor(*this, e.rate_.expr_))
        || boost::apply_visitor(*this, e.ii_.expr_))
        || boost::apply_visitor(*this, e.pMatrix_.expr_))
        || boost::apply_visitor(*this, e.biovar_.expr_))
        || boost::apply_visitor(*this, e.tlag_.expr_));
    }

    bool has_non_param_var_vis::operator()(const generalOdeModel& e)
      const {
      // CHECK - anything to do with nonlinearity ?
      // Putting in the variables that may contain var types
      return ((((((boost::apply_visitor(*this, e.time_.expr_)
        || boost::apply_visitor(*this, e.amt_.expr_))
        || boost::apply_visitor(*this, e.rate_.expr_))
        || boost::apply_visitor(*this, e.ii_.expr_))
        || boost::apply_visitor(*this, e.pMatrix_.expr_))
        || boost::apply_visitor(*this, e.biovar_.expr_))
        || boost::apply_visitor(*this, e.tlag_.expr_));
    }

    bool has_non_param_var_vis::operator()(const pmx_solve_group_control& e)
      const {
      // CHECK - anything to do with nonlinearity ?
      // Putting in the variables that may contain var types
      return ((((((boost::apply_visitor(*this, e.time_.expr_)
        || boost::apply_visitor(*this, e.amt_.expr_))
        || boost::apply_visitor(*this, e.rate_.expr_))
        || boost::apply_visitor(*this, e.ii_.expr_))
        || boost::apply_visitor(*this, e.pMatrix_.expr_))
        || boost::apply_visitor(*this, e.biovar_.expr_))
        || boost::apply_visitor(*this, e.tlag_.expr_));
    }

    bool has_non_param_var_vis::operator()(const pmx_solve_group& e)
      const {
      // CHECK - anything to do with nonlinearity ?
      // Putting in the variables that may contain var types
      return ((((((boost::apply_visitor(*this, e.time_.expr_)
        || boost::apply_visitor(*this, e.amt_.expr_))
        || boost::apply_visitor(*this, e.rate_.expr_))
        || boost::apply_visitor(*this, e.ii_.expr_))
        || boost::apply_visitor(*this, e.pMatrix_.expr_))
        || boost::apply_visitor(*this, e.biovar_.expr_))
        || boost::apply_visitor(*this, e.tlag_.expr_));
    }

    bool has_non_param_var_vis::operator()(const pmx_integrate_ode_group_control& e) const {
      // if any vars, return true because integration will be nonlinear
      return boost::apply_visitor(*this, e.y0_.expr_)
        || boost::apply_visitor(*this, e.theta_.expr_);
    }

    bool has_non_param_var_vis::operator()(const pmx_integrate_ode_group& e) const {
      // if any vars, return true because integration will be nonlinear
      return boost::apply_visitor(*this, e.y0_.expr_)
        || boost::apply_visitor(*this, e.theta_.expr_);
    }

    bool has_non_param_var_vis::operator()(const pmx_integrate_ode_control& e) const {
      // if any vars, return true because integration will be nonlinear
      return boost::apply_visitor(*this, e.y0_.expr_)
        || boost::apply_visitor(*this, e.ts_.expr_)
        || boost::apply_visitor(*this, e.theta_.expr_);
    }

    bool has_non_param_var_vis::operator()(const pmx_integrate_ode& e) const {
      // if any vars, return true because integration will be nonlinear
      return boost::apply_visitor(*this, e.y0_.expr_)
        || boost::apply_visitor(*this, e.ts_.expr_)
        || boost::apply_visitor(*this, e.theta_.expr_);
    }
  }
}


#endif
