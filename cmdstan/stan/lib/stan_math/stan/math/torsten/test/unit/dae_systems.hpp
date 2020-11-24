#ifndef TORSTEN_TEST_DAE_SYSTEMS_HPP
#define TORSTEN_TEST_DAE_SYSTEMS_HPP

struct chemical_kinetics {
  template <typename T0, typename TYY, typename TYP, typename TPAR>
  inline std::vector<typename stan::return_type<TYY, TYP, TPAR>::type>
  operator()(const T0& t_in, const std::vector<TYY>& yy,
             const std::vector<TYP>& yp, const std::vector<TPAR>& theta,
             const std::vector<double>& x_r, const std::vector<int>& x_i,
             std::ostream* msgs) const {
    if (yy.size() != 3 || yp.size() != 3)
      throw std::domain_error(
          "this function was called with inconsistent state");

    std::vector<typename stan::return_type<TYY, TYP, TPAR>::type> res(3);

    auto yy1 = yy.at(0);
    auto yy2 = yy.at(1);
    auto yy3 = yy.at(2);

    auto yp1 = yp.at(0);
    auto yp2 = yp.at(1);

    auto p1 = theta.at(0);
    auto p2 = theta.at(1);
    auto p3 = theta.at(2);

    res[0] = yp1 + p1 * yy1 - p2 * yy2 * yy3;
    res[1] = yp2 - p1 * yy1 + p2 * yy2 * yy3 + p3 * yy2 * yy2;
    res[2] = yy1 + yy2 + yy3 - 1.0;

    return res;
  }
};

struct prey_predator_harvest {
  template <typename T0, typename TYY, typename TYP, typename TPAR>
  inline std::vector<typename stan::return_type<TYY, TYP, TPAR>::type>
  operator()(const T0& t_in, const std::vector<TYY>& yy,
             const std::vector<TYP>& yp, const std::vector<TPAR>& theta,
             const std::vector<double>& x_r, const std::vector<int>& x_i,
             std::ostream* msgs) const {
    if (yy.size() != 3 || yp.size() != 3)
      throw std::domain_error(
          "this function was called with inconsistent state");

    std::vector<typename stan::return_type<TYY, TYP, TPAR>::type> res(3);

    auto yy1 = yy.at(0);
    auto yy2 = yy.at(1);
    auto yy3 = yy.at(2);

    auto yp1 = yp.at(0);
    auto yp2 = yp.at(1);

    constexpr auto r1 = 1.0;
    constexpr auto r2 = 3.0;
    constexpr auto p = 2.0;
    auto mu = theta.at(0);

    res[0] = yp1 - yy1 * (r1 - yy2);
    res[1] = yp2 - yy2 * (r2 - yy2 / yy1 - yy3);
    res[2] = yy3 * (p * yy2 - 1.0) - mu;

    return res;
  }
};

#endif
