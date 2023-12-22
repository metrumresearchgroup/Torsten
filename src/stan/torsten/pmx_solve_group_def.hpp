#ifndef STAN_LANG_AST_NODE_TORSTEN_PMX_SOLVE_GROUP_DEF_HPP
#define STAN_LANG_AST_NODE_TORSTEN_PMX_SOLVE_GROUP_DEF_HPP

#include <stan/lang/ast.hpp>
#include <string>

namespace stan {
  namespace lang {

    pmx_solve_group::pmx_solve_group() { }

    pmx_solve_group::pmx_solve_group (const std::string& integration_function_name,
                                                                              const std::string& system_function_name,
                                                                              const expression& nCmt,
                                                                              const expression& len,
                                                                              const expression& time,
                                                                              const expression& amt,
                                                                              const expression& rate,
                                                                              const expression& ii,
                                                                              const expression& evid,
                                                                              const expression& cmt,
                                                                              const expression& addl,
                                                                              const expression& ss,
                                                                              const expression& pMatrix,
                                                                              const expression& biovar,
                                                                              const expression& tlag)
    : integration_function_name_(integration_function_name),
      system_function_name_(system_function_name),
      nCmt_(nCmt),
      len_(len),
      time_(time), amt_(amt), rate_(rate), ii_(ii),
      evid_(evid), cmt_(cmt), addl_(addl), ss_(ss),
      pMatrix_(pMatrix),
      biovar_(biovar),
      tlag_(tlag)
    {}
  }
}
#endif
