#include <stan/math.hpp>
#include <stan/math/rev/core.hpp>
#include <test/unit/math/rev/mat/fun/util.hpp>
#include <test/unit/math/torsten/pk_cpt_model_test_fixture.hpp>
#include <test/unit/util.hpp>
#include <gtest/gtest.h>

TEST_F(TorstenCptOdeModelTest, 2_cpt_rate_dbl) {
  using stan::math::var;
  using stan::math::to_var;
  using refactor::PKTwoCptModel;
  using torsten::dsolve::pk_integrate_ode_bdf;
  using stan::math::integrate_ode_bdf;
  using refactor::PKTwoCptODE;
  using refactor::PKOdeFunctorRateAdaptor;

  rate[0] = 1200;
  rate[1] = 200;
  rate[2] = 300;
  using model_t = PKTwoCptModel<double, double, double, double>;
  model_t model(t0, y0, rate, CL, Q, V2, V3, ka);
  std::vector<double> yvec(y0.data(), y0.data() + y0.size());
  PKOdeFunctorRateAdaptor<PKTwoCptODE, double> f1(model.f());

  std::vector<double> y = f1(t0, yvec, model.par(), rate, x_i, msgs);
  EXPECT_FLOAT_EQ(y[0], rate[0]);
  EXPECT_FLOAT_EQ(y[1], rate[1]);
  EXPECT_FLOAT_EQ(y[2], rate[2]);
}

TEST_F(TorstenCptOdeModelTest, 2_cpt_rate_var) {
  using stan::math::var;
  using stan::math::to_var;
  using refactor::PKTwoCptModel;
  using torsten::dsolve::pk_integrate_ode_bdf;
  using stan::math::integrate_ode_bdf;
  using refactor::PKTwoCptODE;
  using refactor::PKOdeFunctorRateAdaptor;

  rate[0] = 1200;
  rate[1] = 200;
  rate[2] = 100;
  var CLv = to_var(CL);
  var Qv  = to_var(Q);
  var V2v = to_var(V2);
  var V3v = to_var(V3);
  var kav = to_var(ka);
  std::vector<stan::math::var> rate_var{to_var(rate)};
  using model_t = PKTwoCptModel<double, double, var, var>;
  model_t model(t0, y0, rate_var, CLv, Qv, V2v, V3v, kav);
  std::vector<stan::math::var> theta(model.par());
  std::vector<double> yvec(y0.data(), y0.data() + y0.size());
  PKOdeFunctorRateAdaptor<PKTwoCptODE, var> f1(model.f(), theta.size());
  theta.insert(theta.end(), rate_var.begin(), rate_var.end());

  std::vector<var> y = f1(t0, yvec, theta, x_r, x_i, msgs);
  EXPECT_FLOAT_EQ(y[0].val(), rate[0]);
  EXPECT_FLOAT_EQ(y[1].val(), rate[1]);
  EXPECT_FLOAT_EQ(y[2].val(), rate[2]);
}

TEST_F(TorstenCptOdeModelTest, 2_cpt_solver) {
  using stan::math::var;
  using stan::math::to_var;
  using refactor::PKTwoCptModel;
  using torsten::dsolve::pk_integrate_ode_bdf;
  using stan::math::integrate_ode_bdf;
  using refactor::PKTwoCptODE;
  using refactor::PKOdeFunctorRateAdaptor;

  rate[0] = 1200;
  rate[1] = 200;
  rate[2] = 100;
  y0[0] = 150;
  y0[1] = 50;
  y0[2] = 60;
  ts[0] = 10.0;
  ts.resize(1);
  var CLv = to_var(CL);
  var Qv  = to_var(Q);
  var V2v = to_var(V2);
  var V3v = to_var(V3);
  var kav = to_var(ka);
  std::vector<var> theta{CLv, Qv, V2v, V3v, kav};
  std::vector<stan::math::var> rate_var{to_var(rate)};
  using model_t = PKTwoCptModel<double, double, var, var>;
  model_t model(t0, y0, rate_var, CLv, Qv, V2v, V3v, kav);
  std::vector<double> yvec(y0.data(), y0.data() + y0.size());
  PKOdeFunctorRateAdaptor<PKTwoCptODE, var> f1(model.f(), theta.size());
  theta.insert(theta.end(), rate_var.begin(), rate_var.end());

  auto y1 = pk_integrate_ode_bdf(f1, yvec, t0, ts, theta, x_r, x_i, msgs);
  auto y2 = model.solve(ts[0]);
  EXPECT_FLOAT_EQ(y1[0][0].val(), y2(0).val());
  EXPECT_FLOAT_EQ(y1[0][1].val(), y2(1).val());
  EXPECT_FLOAT_EQ(y1[0][2].val(), y2(2).val());

  std::vector<double> g1, g2;
  for (int i = 0; i < y0.size(); ++i) {
    stan::math::set_zero_all_adjoints();    
    y1[0][i].grad(theta, g1);
    stan::math::set_zero_all_adjoints();    
    y2(i).grad(theta, g2);
    for (size_t j = 0; j < theta.size(); ++j) {
      EXPECT_FLOAT_EQ(g1[j], g2[j]);
    }
  }

  for (int i = 0; i < y0.size(); ++i) {
    stan::math::set_zero_all_adjoints();    
    y1[0][i].grad(rate_var, g1);
    stan::math::set_zero_all_adjoints();    
    y2(i).grad(rate_var, g2);
    for (size_t j = 0; j < rate.size(); ++j) {
      EXPECT_FLOAT_EQ(g1[j], g2[j]);
    }
  }
}

TEST_F(TorstenCptOdeModelTest, 2_cpt_ss_solver_bolus) {
  using stan::math::var;
  using stan::math::to_var;
  using refactor::PKTwoCptModel;
  using torsten::dsolve::pk_integrate_ode_bdf;
  using stan::math::integrate_ode_bdf;
  using refactor::PKTwoCptODE;
  using refactor::PKOdeFunctorRateAdaptor;

  rate[0] = 0;
  rate[1] = 0;
  rate[2] = 0;
  y0[0] = 150;
  y0[1] = 250;
  y0[2] = 350;
  ts[0] = 10.0;
  ts.resize(1);
  var CLv = to_var(CL);
  var Qv  = to_var(Q);
  var V2v = to_var(V2);
  var V3v = to_var(V3);
  var kav = to_var(ka);
  std::vector<stan::math::var> rate_var{to_var(rate)};
  using model_t = PKTwoCptModel<double, double, var, var>;
  model_t model(t0, y0, rate_var, CLv, Qv, V2v, V3v, kav);
  std::vector<var> theta{model.par()};

  //  std::cout.precision(12);

  double amt = 1800;
  int cmt = 1;
  double ii = 12.5;
  
  auto y1 = model.solve(amt, rate_var[cmt - 1], ii, cmt);
  EXPECT_FLOAT_EQ(y1(0).val(), 0.0); // TODO: check solver implementation, why @c y1(0) is always zero
  EXPECT_FLOAT_EQ(y1(1).val(), 738.870108248);
  EXPECT_FLOAT_EQ(y1(2).val(), 763.661542764);

  std::vector<double> g1, g2;
  stan::math::set_zero_all_adjoints();
  y1(0).grad(theta, g1);
  EXPECT_FLOAT_EQ(g1[0], 0.0);
  EXPECT_FLOAT_EQ(g1[1], 0.0);
  EXPECT_FLOAT_EQ(g1[2], 0.0);
  EXPECT_FLOAT_EQ(g1[3], 0.0);
  EXPECT_FLOAT_EQ(g1[4], 0.0);
  stan::math::set_zero_all_adjoints();
  y1(1).grad(theta, g1);
  EXPECT_FLOAT_EQ(g1[0], -106.766204214);
  EXPECT_FLOAT_EQ(g1[1], 1.71716124388 );
  EXPECT_FLOAT_EQ(g1[2], 11.6448699701 );
  EXPECT_FLOAT_EQ(g1[3], 1.25702756718 );
  EXPECT_FLOAT_EQ(g1[4], -33.4223234266);
  stan::math::set_zero_all_adjoints();
  y1(2).grad(theta, g1);
  EXPECT_FLOAT_EQ(g1[0],  -97.6726925784);
  EXPECT_FLOAT_EQ(g1[1],  -2.69921403074);
  EXPECT_FLOAT_EQ(g1[2], 1.70090555664  );
  EXPECT_FLOAT_EQ(g1[3], 13.0890353445  );
  EXPECT_FLOAT_EQ(g1[4],  -34.1903160204);

  cmt = 2;
  y1 = model.solve(amt, rate_var[cmt - 1], ii, cmt);
  
  // std::cout << "taki test: " << y1(0).val() << " " << y1(1).val() << " " << y1(2).val() << "\n";
  // std::cout << "taki test: " << y1(0).val() << " " << y2(1).val() << " " << y2(2).val() << "\n";
  EXPECT_FLOAT_EQ(y1(0).val(), 0.0);
  EXPECT_FLOAT_EQ(y1(1).val(), 700.955088612);
  EXPECT_FLOAT_EQ(y1(2).val(), 763.661017168);

  stan::math::set_zero_all_adjoints();
  y1(0).grad(theta, g1);
  EXPECT_FLOAT_EQ(g1[0], 0.0);
  EXPECT_FLOAT_EQ(g1[1], 0.0);
  EXPECT_FLOAT_EQ(g1[2], 0.0);
  EXPECT_FLOAT_EQ(g1[3], 0.0);
  EXPECT_FLOAT_EQ(g1[4], 0.0);
  stan::math::set_zero_all_adjoints();
  y1(1).grad(theta, g1);
  EXPECT_FLOAT_EQ(g1[0], -104.758709533);
  EXPECT_FLOAT_EQ(g1[1], 1.53235878944 );
  EXPECT_FLOAT_EQ(g1[2], 11.2582785156 );
  EXPECT_FLOAT_EQ(g1[3], 1.48598239969 );
  EXPECT_FLOAT_EQ(g1[4], 0);
  stan::math::set_zero_all_adjoints();
  y1(2).grad(theta, g1);
  EXPECT_FLOAT_EQ(g1[0], -97.6727045238);
  EXPECT_FLOAT_EQ(g1[1], -2.69926906473);
  EXPECT_FLOAT_EQ(g1[2], 1.70091989107 );
  EXPECT_FLOAT_EQ(g1[3], 13.0890426823 );
  EXPECT_FLOAT_EQ(g1[4], -34.1823623977);

  cmt = 3;
  y1 = model.solve(amt, rate_var[cmt - 1], ii, cmt);
  
  // std::cout << "taki test: " << y1(0).val() << " " << y1(1).val() << " " << y1(2).val() << "\n";
  // std::cout << "taki test: " << y1(0).val() << " " << y2(1).val() << " " << y2(2).val() << "\n";
  EXPECT_FLOAT_EQ(y1(0).val(), 0.0);
  EXPECT_FLOAT_EQ(y1(1).val(), 828.127401998);
  EXPECT_FLOAT_EQ(y1(2).val(), 856.228976486);

  stan::math::set_zero_all_adjoints();
  y1(0).grad(theta, g1);
  EXPECT_FLOAT_EQ(g1[0], 0.0);
  EXPECT_FLOAT_EQ(g1[1], 0.0);
  EXPECT_FLOAT_EQ(g1[2], 0.0);
  EXPECT_FLOAT_EQ(g1[3], 0.0);
  EXPECT_FLOAT_EQ(g1[4], 0.0);
  stan::math::set_zero_all_adjoints();
  y1(1).grad(theta, g1);
  EXPECT_FLOAT_EQ(g1[0], -110.006024321);
  EXPECT_FLOAT_EQ(g1[1], -3.07803718766);
  EXPECT_FLOAT_EQ(g1[2], 12.4505184404 );
  EXPECT_FLOAT_EQ(g1[3], 2.71719727476 );
  EXPECT_FLOAT_EQ(g1[4], 0.0);
  stan::math::set_zero_all_adjoints();
  y1(2).grad(theta, g1);
  EXPECT_FLOAT_EQ(g1[0], -99.5058577806);
  EXPECT_FLOAT_EQ(g1[1], -8.28726650911);
  EXPECT_FLOAT_EQ(g1[2], 1.30023459973 );
  EXPECT_FLOAT_EQ(g1[3], 16.044046744  );
  EXPECT_FLOAT_EQ(g1[4], 0.0);
}

TEST_F(TorstenCptOdeModelTest, 2_cpt_ss_solver_multi_trunc_infusion) {
  using stan::math::var;
  using stan::math::to_var;
  using refactor::PKTwoCptModel;
  using torsten::dsolve::pk_integrate_ode_bdf;
  using stan::math::integrate_ode_bdf;
  using refactor::PKTwoCptODE;
  using refactor::PKOdeFunctorRateAdaptor;

  rate[0] = 1200;
  rate[1] = 1100;
  rate[2] = 800;
  y0[0] = 150;
  y0[1] = 250;
  y0[2] = 350;
  ts[0] = 10.0;
  ts.resize(1);
  var CLv = to_var(CL);
  var Qv  = to_var(Q);
  var V2v = to_var(V2);
  var V3v = to_var(V3);
  var kav = to_var(ka);
  std::vector<stan::math::var> rate_var{to_var(rate)};
  using model_t = PKTwoCptModel<double, double, var, var>;
  model_t model(t0, y0, rate_var, CLv, Qv, V2v, V3v, kav);
  std::vector<var> theta(model.par());

  std::cout.precision(12);

  double amt = 1800;
  int cmt = 1;
  double ii = 12.5;
  std::vector<double> g1, g2;

  auto y1 = model.solve(amt, rate_var[cmt - 1], ii, cmt);
  EXPECT_FLOAT_EQ(y1(0).val(), 0.00154469934961);
  EXPECT_FLOAT_EQ(y1(1).val(), 774.104413167);
  EXPECT_FLOAT_EQ(y1(2).val(), 799.879747792);
  // std::cout << "taki test: " << y1(0).val() << " " << y1(1).val() << " " << y1(2).val() << "\n";
  // std::cout << "taki test: " << y1(0).val() << " " << y2(1).val() << " " << y2(2).val() << "\n";

  stan::math::set_zero_all_adjoints();
  y1(0).grad(theta, g1);
  EXPECT_FLOAT_EQ(g1[0], 0.0);
  EXPECT_FLOAT_EQ(g1[1], 0.0);
  EXPECT_FLOAT_EQ(g1[2], 0.0);
  EXPECT_FLOAT_EQ(g1[3], 0.0);
  EXPECT_FLOAT_EQ(g1[4], -0.0178200945892);
  stan::math::set_zero_all_adjoints();
  y1(1).grad(theta, g1);
  EXPECT_FLOAT_EQ(g1[0], -108.530709442);
  EXPECT_FLOAT_EQ(g1[1], 1.88466941803 );
  EXPECT_FLOAT_EQ(g1[2], 11.9988317957 );
  EXPECT_FLOAT_EQ(g1[3], 1.0375686723  );
  EXPECT_FLOAT_EQ(g1[4], -35.1732711317);
  stan::math::set_zero_all_adjoints();
  y1(2).grad(theta, g1);
  EXPECT_FLOAT_EQ(g1[0], -98.8831379736);
  EXPECT_FLOAT_EQ(g1[1], -2.69526994326);
  EXPECT_FLOAT_EQ(g1[2], 1.56750688457 );
  EXPECT_FLOAT_EQ(g1[3], 13.4128341054 );
  EXPECT_FLOAT_EQ(g1[4], -35.6693548262);

  cmt = 2;
  y1 = model.solve(amt, rate_var[cmt - 1], ii, cmt);
  
  EXPECT_FLOAT_EQ(y1(0).val(), 0.0);
  EXPECT_FLOAT_EQ(y1(1).val(), 737.454228957);
  EXPECT_FLOAT_EQ(y1(2).val(), 762.268194176);
  // std::cout << "taki test: " << y1(0).val() << " " << y1(1).val() << " " << y1(2).val() << "\n";
  // std::cout << "taki test: " << y1(0).val() << " " << y2(1).val() << " " << y2(2).val() << "\n";

  stan::math::set_zero_all_adjoints();
  y1(0).grad(theta, g1);
  EXPECT_FLOAT_EQ(g1[0], 0.0);
  EXPECT_FLOAT_EQ(g1[1], 0.0);
  EXPECT_FLOAT_EQ(g1[2], 0.0);
  EXPECT_FLOAT_EQ(g1[3], 0.0);
  EXPECT_FLOAT_EQ(g1[4], 0.0);
  stan::math::set_zero_all_adjoints();
  y1(1).grad(theta, g1);
  EXPECT_FLOAT_EQ(g1[0], -106.758381233);
  EXPECT_FLOAT_EQ(g1[1], 1.71458102804 );
  EXPECT_FLOAT_EQ(g1[2], 11.6337986397 );
  EXPECT_FLOAT_EQ(g1[3], 1.26959503378 );
  EXPECT_FLOAT_EQ(g1[4], 0);
  stan::math::set_zero_all_adjoints();
  y1(2).grad(theta, g1);
  EXPECT_FLOAT_EQ(g1[0], -97.6910811146);
  EXPECT_FLOAT_EQ(g1[1], -2.70645214819);
  EXPECT_FLOAT_EQ(g1[2], 1.71099901438 );
  EXPECT_FLOAT_EQ(g1[3], 13.0830221449 );
  EXPECT_FLOAT_EQ(g1[4], 0.0);

  cmt = 3;
  y1 = model.solve(amt, rate_var[cmt - 1], ii, cmt);
  
  // std::cout << "taki test: " << y1(0).val() << " " << y1(1).val() << " " << y1(2).val() << "\n";
  // std::cout << "taki test: " << y1(0).val() << " " << y2(1).val() << " " << y2(2).val() << "\n";
  EXPECT_FLOAT_EQ(y1(0).val(), 0.0);
  EXPECT_FLOAT_EQ(y1(1).val(), 888.053512871);
  EXPECT_FLOAT_EQ(y1(2).val(), 918.312019204);

  stan::math::set_zero_all_adjoints();
  y1(0).grad(theta, g1);
  EXPECT_FLOAT_EQ(g1[0], 0.0);
  EXPECT_FLOAT_EQ(g1[1], 0.0);
  EXPECT_FLOAT_EQ(g1[2], 0.0);
  EXPECT_FLOAT_EQ(g1[3], 0.0);
  EXPECT_FLOAT_EQ(g1[4], 0.0);
  stan::math::set_zero_all_adjoints();
  y1(1).grad(theta, g1);
  EXPECT_FLOAT_EQ(g1[0], -112.218811747);
  EXPECT_FLOAT_EQ(g1[1], -3.09372149524);
  EXPECT_FLOAT_EQ(g1[2], 12.9947786493 );
  EXPECT_FLOAT_EQ(g1[3], 2.41757181997 );
  EXPECT_FLOAT_EQ(g1[4], 0.0);
  stan::math::set_zero_all_adjoints();
  y1(2).grad(theta, g1);
  EXPECT_FLOAT_EQ(g1[0], -100.772892327);
  EXPECT_FLOAT_EQ(g1[1], -8.7063180167 );
  EXPECT_FLOAT_EQ(g1[2], 1.03257281914 );
  EXPECT_FLOAT_EQ(g1[3], 16.69857146   );
  EXPECT_FLOAT_EQ(g1[4], 0.0);
}

// TEST_F(TorstenCptOdeModelTest, 2_cpt_ss_solver_const_infusion) {
//   using stan::math::var;
//   using stan::math::to_var;
//   using refactor::PKTwoCptModel;
//   using torsten::dsolve::pk_integrate_ode_bdf;
//   using stan::math::integrate_ode_bdf;
//   using refactor::PKTwoCptODE;
//   using refactor::PKOdeFunctorRateAdaptor;

//   rate[0] = 1200;
//   rate[1] = 1100;
//   rate[2] = 800;
//   y0[0] = 150;
//   y0[1] = 250;
//   y0[2] = 350;
//   ts[0] = 10.0;
//   ts.resize(1);
//   var CLv = to_var(CL);
//   var Qv  = to_var(Q);
//   var V2v = to_var(V2);
//   var V3v = to_var(V3);
//   var kav = to_var(ka);
//   std::vector<stan::math::var> rate_var{to_var(rate)};
//   using model_t = PKTwoCptModel<double, double, var, var>;
//   model_t model(t0, y0, rate_var, CLv, Qv, V2v, V3v, kav);
//   std::vector<var> theta(model.par());

//   std::cout.precision(12);

//   double amt = 1800;
//   int cmt = 1;
//   double ii = 0.0;
//   std::vector<double> g1, g2;

//   auto y1 = model.solve(amt, rate_var[cmt - 1], ii, cmt);
//   EXPECT_FLOAT_EQ(y1(0).val(), 1000);
//   EXPECT_FLOAT_EQ(y1(1).val(), 9600);
//   EXPECT_FLOAT_EQ(y1(2).val(), 8400);
//   // std::cout << "taki test: " << y1(0).val() << " " << y1(1).val() << " " << y1(2).val() << "\n";
//   // std::cout << "taki test: " << y1(0).val() << " " << y2(1).val() << " " << y2(2).val() << "\n";

//   stan::math::set_zero_all_adjoints();
//   y1(0).grad(theta, g1);
//   EXPECT_FLOAT_EQ(g1[0], 0.0);
//   EXPECT_FLOAT_EQ(g1[1], 0.0);
//   EXPECT_FLOAT_EQ(g1[2], 0.0);
//   EXPECT_FLOAT_EQ(g1[3], 0.0);
//   EXPECT_FLOAT_EQ(g1[4], -833.333333333);
//   stan::math::set_zero_all_adjoints();
//   y1(1).grad(theta, g1);
//   EXPECT_FLOAT_EQ(g1[0], -960);
//   EXPECT_FLOAT_EQ(g1[1], 2.27373675443e-14 );
//   EXPECT_FLOAT_EQ(g1[2], 120 );
//   EXPECT_FLOAT_EQ(g1[3], -1.81898940355e-14 );
//   EXPECT_FLOAT_EQ(g1[4], 0.0);
//   stan::math::set_zero_all_adjoints();
//   y1(2).grad(theta, g1);
//   EXPECT_FLOAT_EQ(g1[0], -840);
//   EXPECT_FLOAT_EQ(g1[1], 1.70530256582e-13);
//   EXPECT_FLOAT_EQ(g1[2], -7.1054273576e-14);
//   EXPECT_FLOAT_EQ(g1[3], 120 );
//   EXPECT_FLOAT_EQ(g1[4], 1.02318153949e-12);

//   // cmt = 2;
//   // y1 = model.solve(amt, rate_var[cmt - 1], ii, cmt);

//   EXPECT_FLOAT_EQ(y1(0).val(), 1000.0);
//   EXPECT_FLOAT_EQ(y1(1).val(), 9600);
//   EXPECT_FLOAT_EQ(y1(2).val(), 8400);
//   // std::cout << "taki test: " << y1(0).val() << " " << y1(1).val() << " " << y1(2).val() << "\n";
//   // std::cout << "taki test: " << y1(0).val() << " " << y2(1).val() << " " << y2(2).val() << "\n";

//   stan::math::set_zero_all_adjoints();
//   y1(0).grad(theta, g1);
//   EXPECT_FLOAT_EQ(g1[0], 0.0);
//   EXPECT_FLOAT_EQ(g1[1], 0.0);
//   EXPECT_FLOAT_EQ(g1[2], 0.0);
//   EXPECT_FLOAT_EQ(g1[3], 0.0);
//   EXPECT_FLOAT_EQ(g1[4], -833.333333333);
//   stan::math::set_zero_all_adjoints();
//   y1(1).grad(theta, g1);
//   EXPECT_FLOAT_EQ(g1[0], -960);
//   EXPECT_FLOAT_EQ(g1[1], 2.27373675443e-14);
//   EXPECT_FLOAT_EQ(g1[2], 120);
//   EXPECT_FLOAT_EQ(g1[3], -1.81898940355e-14);
//   EXPECT_FLOAT_EQ(g1[4], 0);
//   stan::math::set_zero_all_adjoints();
//   y1(2).grad(theta, g1);
//   EXPECT_FLOAT_EQ(g1[0], -840);
//   EXPECT_FLOAT_EQ(g1[1], 1.70530256582e-13);
//   EXPECT_FLOAT_EQ(g1[2], -7.1054273576e-14);
//   EXPECT_FLOAT_EQ(g1[3], 120);
//   EXPECT_FLOAT_EQ(g1[4], 1.02318153949e-12);

//   cmt = 3;
//   y1 = model.solve(amt, rate_var[cmt - 1], ii, cmt);
  
//   // std::cout << "taki test: " << y1(0).val() << " " << y1(1).val() << " " << y1(2).val() << "\n";
//   // std::cout << "taki test: " << y1(0).val() << " " << y2(1).val() << " " << y2(2).val() << "\n";
//   EXPECT_FLOAT_EQ(y1(0).val(), 0.0);
//   EXPECT_FLOAT_EQ(y1(1).val(), 6400);
//   EXPECT_FLOAT_EQ(y1(2).val(), 7600);

//   stan::math::set_zero_all_adjoints();
//   y1(0).grad(theta, g1);
//   EXPECT_FLOAT_EQ(g1[0], 0.0);
//   EXPECT_FLOAT_EQ(g1[1], 0.0);
//   EXPECT_FLOAT_EQ(g1[2], 0.0);
//   EXPECT_FLOAT_EQ(g1[3], 0.0);
//   EXPECT_FLOAT_EQ(g1[4], 0.0);
//   stan::math::set_zero_all_adjoints();
//   y1(1).grad(theta, g1);
//   EXPECT_FLOAT_EQ(g1[0], -640);
//   EXPECT_FLOAT_EQ(g1[1], -1.18559130767e-13);
//   EXPECT_FLOAT_EQ(g1[2], 80);
//   EXPECT_FLOAT_EQ(g1[3], 1.55913377447e-14);
//   EXPECT_FLOAT_EQ(g1[4], 0);
//   stan::math::set_zero_all_adjoints();
//   y1(2).grad(theta, g1);
//   EXPECT_FLOAT_EQ(g1[0], -560);
//   EXPECT_FLOAT_EQ(g1[1], -71.4285714286);
//   EXPECT_FLOAT_EQ(g1[2], 5.68434188608e-14);
//   EXPECT_FLOAT_EQ(g1[3], 108.571428571);
//   EXPECT_FLOAT_EQ(g1[4], 0.0);
// }
