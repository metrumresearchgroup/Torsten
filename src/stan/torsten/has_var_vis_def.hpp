#ifndef STAN_LANG_TORSTEN_AST_FUN_HAS_VAR_VIS_DEF_HPP
#define STAN_LANG_TORSTEN_AST_FUN_HAS_VAR_VIS_DEF_HPP

#include <stan/lang/ast.hpp>
#include <iostream>

namespace stan {
  namespace lang {

    bool has_var_vis::operator()(const univariate_integral_control& e) const {
      return boost::apply_visitor(*this, e.t0_.expr_)
        || boost::apply_visitor(*this, e.t1_.expr_)
        || boost::apply_visitor(*this, e.theta_.expr_);
    }

    bool has_var_vis::operator()(const generalOdeModel_control& e) const {
      return ((((((boost::apply_visitor(*this, e.time_.expr_)
                   || boost::apply_visitor(*this, e.amt_.expr_))
                  || boost::apply_visitor(*this, e.rate_.expr_))
                 || boost::apply_visitor(*this, e.ii_.expr_))
                || boost::apply_visitor(*this, e.pMatrix_.expr_))
               || boost::apply_visitor(*this, e.biovar_.expr_))
              || boost::apply_visitor(*this, e.tlag_.expr_));
    }

    bool has_var_vis::operator()(const generalOdeModel& e) const {
      return ((((((boost::apply_visitor(*this, e.time_.expr_)
                   || boost::apply_visitor(*this, e.amt_.expr_))
                  || boost::apply_visitor(*this, e.rate_.expr_))
                 || boost::apply_visitor(*this, e.ii_.expr_))
                || boost::apply_visitor(*this, e.pMatrix_.expr_))
               || boost::apply_visitor(*this, e.biovar_.expr_))
              || boost::apply_visitor(*this, e.tlag_.expr_));
    }

    bool has_var_vis::operator()(const pmx_solve_group_control& e) const {
      return ((((((boost::apply_visitor(*this, e.time_.expr_)
                   || boost::apply_visitor(*this, e.amt_.expr_))
                  || boost::apply_visitor(*this, e.rate_.expr_))
                 || boost::apply_visitor(*this, e.ii_.expr_))
                || boost::apply_visitor(*this, e.pMatrix_.expr_))
               || boost::apply_visitor(*this, e.biovar_.expr_))
              || boost::apply_visitor(*this, e.tlag_.expr_));
    }

    bool has_var_vis::operator()(const pmx_solve_group& e) const {
      return ((((((boost::apply_visitor(*this, e.time_.expr_)
                   || boost::apply_visitor(*this, e.amt_.expr_))
                  || boost::apply_visitor(*this, e.rate_.expr_))
                 || boost::apply_visitor(*this, e.ii_.expr_))
                || boost::apply_visitor(*this, e.pMatrix_.expr_))
               || boost::apply_visitor(*this, e.biovar_.expr_))
              || boost::apply_visitor(*this, e.tlag_.expr_));
    }

    bool has_var_vis::operator()(const pmx_integrate_ode_group_control& e) const {
      // only init state and params may contain vars
      return boost::apply_visitor(*this, e.y0_.expr_)
        || boost::apply_visitor(*this, e.theta_.expr_);
    }

    bool has_var_vis::operator()(const pmx_integrate_ode_group& e) const {
      // only init state and params may contain vars
      return boost::apply_visitor(*this, e.y0_.expr_)
        || boost::apply_visitor(*this, e.theta_.expr_);
    }

    bool has_var_vis::operator()(const pmx_integrate_ode_control& e) const {
      // possible vars: init state, time steps, and params
      return boost::apply_visitor(*this, e.y0_.expr_)
        || boost::apply_visitor(*this, e.ts_.expr_)
        || boost::apply_visitor(*this, e.theta_.expr_);
    }

    bool has_var_vis::operator()(const pmx_integrate_ode& e) const {
      // possible vars: init state, time steps, and params
      return boost::apply_visitor(*this, e.y0_.expr_)
        || boost::apply_visitor(*this, e.ts_.expr_)
        || boost::apply_visitor(*this, e.theta_.expr_);
    }

  }
}
#endif
