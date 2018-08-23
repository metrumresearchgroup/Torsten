#ifndef TORSTEN_IS_DETECTED_HPP
#define TORSTEN_IS_DETECTED_HPP

/* torsten's own @c void_t type and detect idiom */
namespace torsten {

  template <class... Ts>
  using void_t = void;

  struct nonesuch {
    nonesuch() = delete;
    ~nonesuch() = delete;
    nonesuch(nonesuch const&) = delete;
    void operator=(nonesuch const&) = delete;
  };

  namespace detail {
    template <class Default, class AlwaysVoid,
              template<class...> class Op, class... Args>
    struct detector {
      using value_t = std::false_type;
      using type = Default;
    };
 
    template <class Default, template<class...> class Op, class... Args>
    struct detector<Default, void_t<Op<Args...>>, Op, Args...> {
      // Note that std::void_t is a C++17 feature
      using value_t = std::true_type;
      using type = Op<Args...>;
    };
 
  } // namespace detail
 
  template <template<class...> class Op, class... Args>
  using is_detected = typename detail::detector<nonesuch, void, Op, Args...>::value_t;
 
  template <template<class...> class Op, class... Args>
  using detected_t = typename detail::detector<nonesuch, void, Op, Args...>::type;
 
  template <class Default, template<class...> class Op, class... Args>
  using detected_or = detail::detector<Default, void, Op, Args...>;  

}

#endif
