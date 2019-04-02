#ifndef STAN_LANG_TORSTEN_GRAMMARS_SEMANTIC_ACTIONS_DATA_ONLY_EXPRESSION_HPP
#define STAN_LANG_TORSTEN_GRAMMARS_SEMANTIC_ACTIONS_DATA_ONLY_EXPRESSION_HPP

      bool operator()(const univariate_integral_control& /*x*/) const;
      bool operator()(const generalOdeModel_control& /*x*/) const;
      bool operator()(const generalOdeModel& /*x*/) const;
      bool operator()(const pmx_solve_group& /*x*/) const;
#endif
