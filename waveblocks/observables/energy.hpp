#pragma once

#include "../types.hpp"
#include "../wavepackets/hawp_commons.hpp"
#include "../wavepackets/hawp_gradient_operator.hpp"
#include "../wavepackets/hawp_gradient_evaluator.hpp"
#include "../innerproducts/homogeneous_inner_product.hpp"


namespace waveblocks {
    namespace observables {
        using wavepackets::ScalarHaWp;
        using wavepackets::HaWpGradient;
        using wavepackets::HaWpGradientOperator;
        using innerproducts::HomogeneousInnerProduct;
        using utilities::Squeeze;

        /**
         * \brief Computes potential energy of a Hagedorn Wavepacket.
         *
         * \param packet
         * \param V potential
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
                                 [&V] (const CMatrix<D,Eigen::Dynamic>& nodes,
                                       const RMatrix<D,1>& pos)
                                 -> CMatrix<1,Eigen::Dynamic> {
                                     (void)pos; // Unused
                                     const dim_t n_nodes = nodes.cols();
                                     CMatrix<1,Eigen::Dynamic> result(1, n_nodes);
                                     for(int i = 0; i < n_nodes; ++i) {
                                         result(0,i) = V.evaluate_at(Squeeze<D,CMatrix<D,Eigen::Dynamic>>::apply(nodes,i));
                                     }
                                     return result;
                                 }
                                 ).real();
        }

        /**
         * \brief Computes kinetic energy of a Hagedorn Wavepacket.
         *
         * \param packet
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

        template<int D,class MultiIndex>
        real_t norm(const ScalarHaWp<D,MultiIndex>& packet)
        {
            complex_t result(0,0);
            result += packet.coefficients().dot(packet.coefficients());
            return 0.5*result.real();
        }
    }
}
