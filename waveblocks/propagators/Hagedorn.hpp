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
        using steps::Step1;
        using steps::Step2;
        using steps::Step3;


        /**
         * \brief Models the Hagedorn propagator. Offers propagation method.
         */
        template <int N, int D, class MultiIndex, class TQR>
        struct Hagedorn {
            template<class Potential>
            static void propagate(InhomogeneousHaWp<D,MultiIndex> &packet,
                                  const real_t &delta_t,
                                  const Potential &V) {
                int i = 0;
                for (auto& component : packet.components()) {

                    auto& params = component.parameters();
                    Step1<N,D>::apply(params,delta_t);
                    Step2<N,D>::inhomogeneous(i,V,params,delta_t);
                    i++;
                }
                Step3<InhomogeneousHaWp<D,MultiIndex>,Potential,N,D,VectorInnerProduct<D,MultiIndex,TQR>>::apply(packet, V, delta_t);

                i = 0;
                for (auto& component : packet.components()) {
                    auto& params = component.parameters();
                    Step1<N,D>::apply(params,delta_t);
                    i++;
                }
            }

            template<class Potential>
            static void propagate(HomogeneousHaWp<D,MultiIndex> &packet,
                                  const real_t &delta_t,
                                  const Potential &V
                                  ) {
                auto& params = packet.parameters();
                Step1<N,D>::apply(params, delta_t);
                Step2<N,D>::homogeneous(V,params,delta_t);
                Step3<HomogeneousHaWp<D,MultiIndex>,Potential,N,D,VectorInnerProduct<D,MultiIndex,TQR>>::apply(packet, V, delta_t);
                Step1<N,D>::apply(params,delta_t);

            }
        };
        template<int D, class MultiIndex, class TQR>
        struct Hagedorn<1,D,MultiIndex,TQR> {
            template<class Potential>
            static void propagate(ScalarHaWp<D, MultiIndex> &packet,
                                  const real_t &delta_t,
                                  const Potential &V) {
                auto& params = packet.parameters();
                Step1<1,D>::apply(params, delta_t);
                Step2<1,D>::homogeneous(V,params,delta_t);
                Step3<ScalarHaWp<D,MultiIndex>,Potential,1,D,HomogeneousInnerProduct<D,MultiIndex,TQR>>::apply(packet, V, delta_t);
                Step1<1,D>::apply(params,delta_t);
            }
        };
    }
}
