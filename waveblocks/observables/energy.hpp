#include "../types.hpp"
#include "../hawp_commons.hpp"
#include "../homogeneous_inner_product.hpp"
#include "../hawp_gradient_operator.hpp"
#include "../hawp_gradient_evaluator.hpp"

namespace waveblocks {
    using utilities::Squeeze;
    
    /**
    * \brief Computes potential energy of a Hagedorn Wavepacket.
    *
    * \param
    * 
    * \tparam Potential
    * Needs to implement evaluation::Abstract interface
    * \tparam D
    * Dimension of argument space
    * \tparam MultiIndex
    * \tparam TQR    *
    */
    template<class Potential, int D, class MultiIndex, class TQR>
    real_t potential_energy(const ScalarHaWp<D, MultiIndex>& packet,
                            const Potential& V) {
        HomogeneousInnerProduct<D, MultiIndex, TQR> ip;
        return ip.quadrature(packet,
                                  [&V] (const CMatrix<D,Eigen::Dynamic>& nodes, const CMatrix<D,1>&) -> CMatrix<1,Eigen::Dynamic> {
                                      const dim_t n_nodes = nodes.cols();
                                      CMatrix<1,Eigen::Dynamic> result(1, n_nodes);
                                      for(int i = 0; i < n_nodes; ++i)  {
                                          result(0,i) = V.evaluate_at(Squeeze<D,CMatrix<D,Eigen::Dynamic>>::apply(nodes,i));
                                      }
                                      return result;
                                  }).real();
    }

    /**
    * \brief Computes kinetic energy of a Hagedorn Wavepacket.
    * 
    * \tparam D
    * Dimension of argument space
    * \tparam MultiIndex
    */
    template<int D, class MultiIndex>
    real_t kinetic_energy(const ScalarHaWp<D, MultiIndex>& packet) {
        HaWpGradientOperator<D,MultiIndex> nabla;
        HaWpGradient<D,MultiIndex> gradwp = nabla(packet);
        complex_t result(0,0);
        for (size_t i = 0 ; i < gradwp.n_components(); ++i) {
            result += gradwp.component(i).coefficients().dot(gradwp.component(i).coefficients());
        }
        return 0.5 * result.real();
    }
  
}
