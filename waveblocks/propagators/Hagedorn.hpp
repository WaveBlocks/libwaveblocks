#pragma once

#include "../types.hpp"
#include "../matrixPotentials/potentials.hpp"
#include "../matrixPotentials/bases.hpp"
#include "../homogeneous_inner_product.hpp"
#include "../vector_inner_product.hpp"
#include "basic_steps.hpp"


namespace waveblocks
{
    namespace propagators
    {
        using innerproducts::HomogeneousInnerProduct;
        using innerproducts::VectorInnerProduct;
        using steps::StepT;
        using steps::StepU;
        using steps::StepW;

        /**
         * \brief Implements the Hagedorn propagator for
         * vector valued wavepackets. Offers a method for
         * time propagation.
         *
         * \tparam N
         * Number of levels
         * \tparam D
         * Dimension of space
         * \tparam MultiIndex
         * Type of multi index used in the basis shape
         * \tparam MDQR
         * Multi-dimensional quadrature rule
         */
        template <int N, int D, class MultiIndex, class MDQR>
        struct Hagedorn {
            // Inhomogeneous wavepackets
            template<class Potential>
            static void propagate(InhomogeneousHaWp<D, MultiIndex> &packet,
                                  const real_t &delta_t,
                                  const Potential &V) {
                int i = 0;
                for (auto& component : packet.components()) {
                    auto& params = component.parameters();
                    StepT<N,D>::apply(params, delta_t);
                    StepU<N,D>::inhomogeneous(i, V, params, delta_t);
                    i++;
                }
                StepW<InhomogeneousHaWp<D, MultiIndex>, Potential, N, D, VectorInnerProduct<D, MultiIndex, MDQR> >::apply(packet, V, delta_t);
                for (auto& component : packet.components()) {
                    auto& params = component.parameters();
                    StepT<N,D>::apply(params, delta_t);
                }
            }

            // Homogeneous wavepackets
            template<class Potential>
            static void propagate(HomogeneousHaWp<D,MultiIndex> &packet,
                                  const real_t &delta_t,
                                  const Potential &V) {
                auto& params = packet.parameters();
                StepT<N,D>::apply(params, delta_t);
                StepU<N,D>::homogeneous(V, params, delta_t);
                StepW<HomogeneousHaWp<D, MultiIndex>, Potential, N, D, VectorInnerProduct<D, MultiIndex, MDQR> >::apply(packet, V, delta_t);
                StepT<N,D>::apply(params, delta_t);
            }
        };

        /**
         * \brief Implements the Hagedorn propagator for
         * scalar wavepackets. Offers a method for
         * time propagation.
         *
         * \tparam D
         * Dimension of space
         * \tparam MultiIndex
         * Type of multi index used in the basis shape
         * \tparam MDQR
         * Multi-dimensional quadrature rule
         */
        template<int D, class MultiIndex, class MDQR>
        struct Hagedorn<1, D, MultiIndex, MDQR> {
            // Scalar wavepackets
            template<class Potential>
            static void propagate(ScalarHaWp<D, MultiIndex> &packet,
                                  const real_t &delta_t,
                                  const Potential &V) {
                auto& params = packet.parameters();
                StepT<1,D>::apply(params, delta_t);
                StepU<1,D>::homogeneous(V, params, delta_t);
                StepW<ScalarHaWp<D, MultiIndex>, Potential, 1, D, HomogeneousInnerProduct<D, MultiIndex, MDQR> >::apply(packet, V, delta_t);
                StepT<1,D>::apply(params, delta_t);
            }
        };
    }
}
