#include "../types.hpp"
#include "../hawp_commons.hpp"
#include "../homogeneous_inner_product.hpp"
#include "../hawp_gradient_operator.hpp"
#include "../hawp_gradient_evaluator.hpp"

namespace waveblocks {

  template<int D, class T>
  struct helper {
    static const T& apply(const T& in) {return in;}
  };

  template<class T>
  struct helper<1,T> {
    static const typename std::remove_reference<decltype(std::declval<T>()[0])>::type& apply(const T& in) {return in[0];}
  };

  template<class Potential, int D, class MultiIndex, class TQR>
  real_t potential_energy(const ScalarHaWp<D, MultiIndex>& packet, const Potential& V) {
    InhomogeneousInnerProduct<D, MultiIndex, TQR> ip;
    auto INHOM =  ip.quadrature(packet, packet, [&V] (const CMatrix<D,Eigen::Dynamic>& nodes, const RMatrix<D,1>&) -> CMatrix<1,Eigen::Dynamic> {
      const dim_t n_nodes = nodes.cols();
      CMatrix<1,Eigen::Dynamic> result(1, n_nodes);
      for(int i = 0; i < n_nodes; ++i)  {
        result(0,i) = V.evaluate_at(HelperArg<D>::first(nodes,i));
      }
      return result;
    }).real();
    HomogeneousInnerProduct<D, MultiIndex, TQR> ip2;
    auto HOM = ip2.quadrature(packet, [&V] (const CMatrix<D,Eigen::Dynamic>& nodes, const CMatrix<D,1>&) -> CMatrix<1,Eigen::Dynamic> {
      const dim_t n_nodes = nodes.cols();
      CMatrix<1,Eigen::Dynamic> result(1, n_nodes);
      for(int i = 0; i < n_nodes; ++i)  {
        result(0,i) = V.evaluate_at(HelperArg<D>::first(nodes,i));
      }
      return result;
    }).real();
    std::cout<<"HOM - INHOM: " << HOM - INHOM << std::endl;
    return HOM;
  }

  template<int D, class MultiIndex>
  real_t kinetic_energy(const ScalarHaWp<D, MultiIndex>& packet) {
    HaWpGradientOperator<D,MultiIndex> nabla;
    HaWpGradient<D,MultiIndex> gradwp = nabla(packet);
    complex_t result(0,0);
    for (size_t i = 0 ; i < gradwp.n_components(); ++i) {
        result += gradwp.component(i).coefficients().dot(gradwp.component(i).coefficients());
        std::cout<< gradwp.component(i).coefficients().rows() << ", " << gradwp.component(i).coefficients().cols() << std::endl;
    }
    return 0.5 * result.real();
  }
}
