namespace waveblocks {
  namespace utilities {
    template<int D, class T>
    struct Squeeze {
      static const T& apply(const T& in) {return in;}
      static T& apply(T& in) {return in;}
    };

    template<class T>
    struct Squeeze<1,T> {
      static const typename std::remove_reference<decltype(std::declval<T>()[0])>::type& apply(const T& in) {return in[0];}
      static decltype(std::declval<T>()[0]) apply(T& in) {return in[0];}
    };

    template<int D>
    struct Squeeze<D, CMatrix<D,Eigen::Dynamic>> {
      static CMatrix<D,1> apply(const CMatrix<D,Eigen::Dynamic>& nodes, int l) {
        return nodes.template block<D,1>(0,l);
      }
    };

    template<>
    struct Squeeze<1, CMatrix<1,Eigen::Dynamic>> {
      static const complex_t& apply(const CMatrix<1,Eigen::Dynamic>& nodes, int l) {
        return nodes(0,l);
      }
    };

    template<int D, class T>
    struct Unsqueeze {
      static const T& apply(const T& t) {
        return t;
      }
    };

    template<class T>
    struct Unsqueeze<1,T> {
      static T apply(const typename std::remove_reference<decltype(std::declval<T>()[0])>::type& t) {
        T result;
        result[0] = t;
        return result;
      }
    };
  }
}
