#ifndef STAN_LANG_TORSTEN_GRAMMARS_TERM_GRAMMAR_DEF_HPP
#define STAN_LANG_TORSTEN_GRAMMARS_TERM_GRAMMAR_DEF_HPP

univariate_integral_control_r.name("expression");
univariate_integral_control_r
%= ( (string("univariate_integral_rk45") >>
      no_skip[!char_("a-zA-Z0-9_")])
     | (string("univariate_integral_bdf") >>
        no_skip[!char_("a-zA-Z0-9_")]) )
  > lit('(')
  > identifier_r          // 1) system function name (function only)
  > lit(',')
  > expression_g(_r1)     // 2) t0 (data only)
  > lit(',')
  > expression_g(_r1)     // 2) t1 (data only)
  > lit(',')
  > expression_g(_r1)     // 3) theta
  > lit(',')
  > expression_g(_r1)     // 4) x_r (data only)
  > lit(',')
  > expression_g(_r1)     // 5) x_i (data only)
  > lit(')')
  [validate_univariate_integral_control_f(_val,
                                          boost::phoenix::ref(var_map_),
                                          _pass,
                                          boost::phoenix::ref(error_msgs_))];

generalOdeModel_control_r.name("expression");
generalOdeModel_control_r
%= ( (string("generalOdeModel_bdf") >> no_skip[!char_("a-zA-Z0-9_")])
     | (string("generalOdeModel_rk45") >> no_skip[!char_("a-zA-Z0-9_")])
     | (string("mixOde1CptModel_rk45") >> no_skip[!char_("a-zA-Z0-9_")])
     | (string("mixOde1CptModel_bdf") >> no_skip[!char_("a-zA-Z0-9_")])
     | (string("mixOde2CptModel_rk45") >> no_skip[!char_("a-zA-Z0-9_")])
     | (string("mixOde2CptModel_bdf") >> no_skip[!char_("a-zA-Z0-9_")])
     | (string("pmx_solve_adams") >> no_skip[!char_("a-zA-Z0-9_")])
     | (string("pmx_solve_bdf") >> no_skip[!char_("a-zA-Z0-9_")])
     | (string("pmx_solve_rk45") >> no_skip[!char_("a-zA-Z0-9_")])
     | (string("pmx_solve_onecpt_adams") >> no_skip[!char_("a-zA-Z0-9_")])
     | (string("pmx_solve_onecpt_bdf") >> no_skip[!char_("a-zA-Z0-9_")])
     | (string("pmx_solve_onecpt_rk45") >> no_skip[!char_("a-zA-Z0-9_")])
     | (string("pmx_solve_twocpt_adams") >> no_skip[!char_("a-zA-Z0-9_")])
     | (string("pmx_solve_twocpt_bdf") >> no_skip[!char_("a-zA-Z0-9_")])
     | (string("pmx_solve_twocpt_rk45") >> no_skip[!char_("a-zA-Z0-9_")])
     )
  >> lit('(')            // >> allows backtracking to non-control
  >> identifier_r        // 1) system function name (function only)
  >> lit(',')
  >> expression_g(_r1)   // 2) nCmt
  >> lit(',')
  >> expression_g(_r1)   // 3) time
  >> lit(',')
  >> expression_g(_r1)   // 4) amt
  >> lit(',')
  >> expression_g(_r1)   // 5) rate
  >> lit(',')
  >> expression_g(_r1)   // 6) ii
  >> lit(',')
  >> expression_g(_r1)   // 7) evid (data only)
  >> lit(',')
  >> expression_g(_r1)   // 8) cmt (data only)
  >> lit(',')
  >> expression_g(_r1)   // 9) addl (data only)
  >> lit(',')
  >> expression_g(_r1)   // 10) ss (data only)
  >> lit(',')
  >> expression_g(_r1)   // 11) pMatrix
  >> lit(',')
  >> expression_g(_r1)   // 12) biovar
  >> lit(',')
  >> expression_g(_r1)   // 13) tlag
  >> lit(',')
  >> expression_g(_r1)   // 14) relative tolerance (data only)
  >> lit(',')
  >> expression_g(_r1)   // 15) absolute tolerance (data only)
  >> lit(',')
  >> expression_g(_r1)   // 16) maximum number of steps
  > lit(')')
  [validate_generalOdeModel_control_f(_val,
                                      boost::phoenix::ref(var_map_), _pass,
                                      boost::phoenix::ref(error_msgs_))];

generalOdeModel_r.name("expression");
generalOdeModel_r
%= ( (string("generalOdeModel_bdf") >> no_skip[!char_("a-zA-Z0-9_")])
     | (string("generalOdeModel_rk45") >> no_skip[!char_("a-zA-Z0-9_")])
     | (string("mixOde1CptModel_rk45") >> no_skip[!char_("a-zA-Z0-9_")])
     | (string("mixOde1CptModel_bdf") >> no_skip[!char_("a-zA-Z0-9_")])
     | (string("mixOde2CptModel_rk45") >> no_skip[!char_("a-zA-Z0-9_")])
     | (string("mixOde2CptModel_bdf") >> no_skip[!char_("a-zA-Z0-9_")])
     | (string("pmx_solve_adams") >> no_skip[!char_("a-zA-Z0-9_")])
     | (string("pmx_solve_bdf") >> no_skip[!char_("a-zA-Z0-9_")])
     | (string("pmx_solve_rk45") >> no_skip[!char_("a-zA-Z0-9_")])
     | (string("pmx_solve_onecpt_adams") >> no_skip[!char_("a-zA-Z0-9_")])
     | (string("pmx_solve_onecpt_bdf") >> no_skip[!char_("a-zA-Z0-9_")])
     | (string("pmx_solve_onecpt_rk45") >> no_skip[!char_("a-zA-Z0-9_")])
     | (string("pmx_solve_twocpt_adams") >> no_skip[!char_("a-zA-Z0-9_")])
     | (string("pmx_solve_twocpt_bdf") >> no_skip[!char_("a-zA-Z0-9_")])
     | (string("pmx_solve_twocpt_rk45") >> no_skip[!char_("a-zA-Z0-9_")])
     )
  > lit('(')
  > identifier_r        // 1) system function name (function only)
  > lit(',')
  > expression_g(_r1)   // 2) nCmt
  > lit(',')
  > expression_g(_r1)   // 3) time
  > lit(',')
  > expression_g(_r1)   // 4) amt
  > lit(',')
  > expression_g(_r1)   // 5) rate
  > lit(',')
  > expression_g(_r1)   // 6) ii
  > lit(',')
  > expression_g(_r1)   // 7) evid (data only)
  > lit(',')
  > expression_g(_r1)   // 8) cmt (data only)
  > lit(',')
  > expression_g(_r1)   // 9) addl (data only)
  > lit(',')
  > expression_g(_r1)   // 10) ss (data only)
  > lit(',')
  > expression_g(_r1)   // 11) pMatrix
  > lit(',')
  > expression_g(_r1)   // 12) biovar
  > lit(',')
  > expression_g(_r1)   // 13) tlag
  > lit(')')
  [validate_generalOdeModel_f(_val,
                              boost::phoenix::ref(var_map_), _pass,
                              boost::phoenix::ref(error_msgs_))];

pmx_solve_group_control_r.name("expression");
pmx_solve_group_control_r
%= (   (string("pmx_solve_group_adams") >> no_skip[!char_("a-zA-Z0-9_")])
     | (string("pmx_solve_group_bdf")   >> no_skip[!char_("a-zA-Z0-9_")])
     | (string("pmx_solve_group_rk45")  >> no_skip[!char_("a-zA-Z0-9_")])
     )
  >> lit('(')            // >> allows backtracking to non-control
  >> identifier_r        // 1) system function name (function only)
  >> lit(',')
  >> expression_g(_r1)   // 2) nCmt
  >> lit(',')
  >> expression_g(_r1)   // 3) len
  >> lit(',')
  >> expression_g(_r1)   // 4) time
  >> lit(',')
  >> expression_g(_r1)   // 5) amt
  >> lit(',')
  >> expression_g(_r1)   // 6) rate
  >> lit(',')
  >> expression_g(_r1)   // 7) ii
  >> lit(',')
  >> expression_g(_r1)   // 8) evid (data only)
  >> lit(',')
  >> expression_g(_r1)   // 9) cmt (data only)
  >> lit(',')
  >> expression_g(_r1)   // 10) addl (data only)
  >> lit(',')
  >> expression_g(_r1)   // 11) ss (data only)
  >> lit(',')
  >> expression_g(_r1)   // 12) pMatrix
  >> lit(',')
  >> expression_g(_r1)   // 13) biovar
  >> lit(',')
  >> expression_g(_r1)   // 14) tlag
  >> lit(',')
  >> expression_g(_r1)   // 15) relative tolerance (data only)
  >> lit(',')
  >> expression_g(_r1)   // 16) absolute tolerance (data only)
  >> lit(',')
  >> expression_g(_r1)   // 17) maximum number of steps
  > lit(')')
  [validate_pmx_solve_group_control_f(_val,
                                      boost::phoenix::ref(var_map_), _pass,
                                      boost::phoenix::ref(error_msgs_))];

pmx_solve_group_r.name("expression");
pmx_solve_group_r
%= ( (string("pmx_solve_group_rk45") >> no_skip[!char_("a-zA-Z0-9_")])
     | (string("pmx_solve_group_adams") >> no_skip[!char_("a-zA-Z0-9_")])
     | (string("pmx_solve_group_bdf") >> no_skip[!char_("a-zA-Z0-9_")])
     | (string("pmx_solve_group_onecpt_bdf") >> no_skip[!char_("a-zA-Z0-9_")])
     | (string("pmx_solve_group_onecpt_rk45") >> no_skip[!char_("a-zA-Z0-9_")])
     | (string("pmx_solve_group_twocpt_bdf") >> no_skip[!char_("a-zA-Z0-9_")])
     | (string("pmx_solve_group_twocpt_rk45") >> no_skip[!char_("a-zA-Z0-9_")])
     )
  > lit('(')
  > identifier_r        // 1) system function name (function only)
  > lit(',')
  > expression_g(_r1)   // 2) nCmt
  > lit(',')
  > expression_g(_r1)   // 3) len
  > lit(',')
  > expression_g(_r1)   // 4) time
  > lit(',')
  > expression_g(_r1)   // 5) amt
  > lit(',')
  > expression_g(_r1)   // 6) rate
  > lit(',')
  > expression_g(_r1)   // 7) ii
  > lit(',')
  > expression_g(_r1)   // 8) evid (data only)
  > lit(',')
  > expression_g(_r1)   // 9) cmt (data only)
  > lit(',')
  > expression_g(_r1)   // 10) addl (data only)
  > lit(',')
  > expression_g(_r1)   // 11) ss (data only)
  > lit(',')
  > expression_g(_r1)   // 13) pMatrix
  > lit(',')
  > expression_g(_r1)   // 15) biovar
  > lit(',')
  > expression_g(_r1)   // 17) tlag
  > lit(')')
  [validate_pmx_solve_group_f(_val,
                              boost::phoenix::ref(var_map_), _pass,
                              boost::phoenix::ref(error_msgs_))];

pmx_integrate_ode_group_control_r.name("expression");
pmx_integrate_ode_group_control_r
%= (   (string("pmx_integrate_ode_group_rk45")  >> no_skip[!char_("a-zA-Z0-9_")])
     | (string("pmx_integrate_ode_group_bdf")   >> no_skip[!char_("a-zA-Z0-9_")])
     | (string("pmx_integrate_ode_group_adams") >> no_skip[!char_("a-zA-Z0-9_")])
     )
  >> lit('(')              // >> allows backtracking to non-control
  >> identifier_r          // 1) system function name (function only)
  >> lit(',')
  >> expression_g(_r1)     // 2) y0
  >> lit(',')
  >> expression_g(_r1)     // 3) t0 (data only)
  >> lit(',')
  >> expression_g(_r1)     // 4) len(data only)
  >> lit(',')
  >> expression_g(_r1)     // 5) ts (data only)
  >> lit(',')
  >> expression_g(_r1)     // 6) theta
  >> lit(',')
  >> expression_g(_r1)     // 7) x (data only)
  >> lit(',')
  >> expression_g(_r1)     // 8) x_int (data only)
  >> lit(',')
  >> expression_g(_r1)     // 9) relative tolerance (data only)
  >> lit(',')
  >> expression_g(_r1)     // 10) absolute tolerance (data only)
  >> lit(',')
  >> expression_g(_r1)     // 11) maximum number of steps (data only)
  > lit(')')
  [validate_pmx_integrate_ode_group_control_f(_val, boost::phoenix::ref(var_map_),
                                        _pass, boost::phoenix::ref(error_msgs_))];

pmx_integrate_ode_group_r.name("expression");
pmx_integrate_ode_group_r
%= (   (string("pmx_integrate_ode_group_rk45")  >> no_skip[!char_("a-zA-Z0-9_")])
     | (string("pmx_integrate_ode_group_bdf")   >> no_skip[!char_("a-zA-Z0-9_")])
     | (string("pmx_integrate_ode_group_adams") >> no_skip[!char_("a-zA-Z0-9_")])
     )
  > lit('(')
  > identifier_r          // 1) system function name (function only)
  > lit(',')
  > expression_g(_r1)     // 2) y0
  > lit(',')
  > expression_g(_r1)     // 3) t0 (data only)
  > lit(',')
  > expression_g(_r1)     // 4) len(data only)
  > lit(',')
  > expression_g(_r1)     // 5) ts (data only)
  > lit(',')
  > expression_g(_r1)     // 6) theta
  > lit(',')
  > expression_g(_r1)     // 7) x (data only)
  > lit(',')
  > expression_g(_r1)     // 8) x_int (data only)
  > lit(')')
  [validate_pmx_integrate_ode_group_f(_val, boost::phoenix::ref(var_map_),
                                      _pass, boost::phoenix::ref(error_msgs_))];

pmx_integrate_ode_control_r.name("expression");
pmx_integrate_ode_control_r
%= (   (string("pmx_integrate_ode_rk45")  >> no_skip[!char_("a-zA-Z0-9_")])
     | (string("pmx_integrate_ode_bdf")   >> no_skip[!char_("a-zA-Z0-9_")])
     | (string("pmx_integrate_ode_adams") >> no_skip[!char_("a-zA-Z0-9_")])
     )
  >> lit('(')              // >> allows backtracking to non-control
  >> identifier_r          // 1) system function name (function only)
  >> lit(',')
  >> expression_g(_r1)     // 2) y0
  >> lit(',')
  >> expression_g(_r1)     // 3) t0 (data only)
  >> lit(',')
  >> expression_g(_r1)     // 4) ts
  >> lit(',')
  >> expression_g(_r1)     // 5) theta
  >> lit(',')
  >> expression_g(_r1)     // 6) x (data only)
  >> lit(',')
  >> expression_g(_r1)     // 7) x_int (data only)
  >> lit(',')
  >> expression_g(_r1)     // 8) relative tolerance (data only)
  >> lit(',')
  >> expression_g(_r1)     // 9) absolute tolerance (data only)
  >> lit(',')
  >> expression_g(_r1)     // 10) maximum number of steps (data only)
  > lit(')')
  [validate_pmx_integrate_ode_control_f(_val, boost::phoenix::ref(var_map_),
                                        _pass, boost::phoenix::ref(error_msgs_))];

pmx_integrate_ode_r.name("expression");
pmx_integrate_ode_r
%= (   (string("pmx_integrate_ode_rk45")  >> no_skip[!char_("a-zA-Z0-9_")])
     | (string("pmx_integrate_ode_bdf")   >> no_skip[!char_("a-zA-Z0-9_")])
     | (string("pmx_integrate_ode_adams") >> no_skip[!char_("a-zA-Z0-9_")])
     )
  > lit('(')
  > identifier_r          // 1) system function name (function only)
  > lit(',')
  > expression_g(_r1)     // 2) y0
  > lit(',')
  > expression_g(_r1)     // 3) t0 (data only)
  > lit(',')
  > expression_g(_r1)     // 4) ts
  > lit(',')
  > expression_g(_r1)     // 5) theta
  > lit(',')
  > expression_g(_r1)     // 6) x (data only)
  > lit(',')
  > expression_g(_r1)     // 7) x_int (data only)
  > lit(')')
  [validate_pmx_integrate_ode_f(_val, boost::phoenix::ref(var_map_),
                                _pass, boost::phoenix::ref(error_msgs_))];


#endif
