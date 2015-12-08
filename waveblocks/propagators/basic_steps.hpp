#pragma once

#include "../types.hpp"
#include "../utilities/squeeze.hpp"
#include "../utilities/adaptors.hpp"


namespace waveblocks
{
    namespace propagators
    {
        namespace steps
        {
            using utilities::Squeeze;
            using utilities::PacketToCoefficients;
            using utilities::Unsqueeze;


            template<int N, int D>
            struct Step1 {
                static void apply(HaWpParamSet<D>& params,
                                  const real_t &delta_t) {
                    params.updateq( 0.5 * delta_t * params.p() );
                    params.updateQ( 0.5 * delta_t * params.P() );
                    params.updateS( 0.25 * delta_t * params.p().dot(params.p()) );
                }
            };

            /**
             * \brief Helper class for Step2. If Mode then level is supscripted in i.
             *
             * \tparam Mode flag true for inhomogeneous computation
             */
            template<bool Mode>
            struct HelperL {
                template<class L>
                static typename decltype(std::declval<L>()[0])::type apply(const L& level, int i) {
                    return level[i];
                }
            };

            template<>
            struct HelperL<false> {
                template<class L>
                static const L& apply(const L& level, int) {
                    return level;
                }
            };

            /**
             * \brief Performs commong code of Step2 for all specializations.
             *
             * \tparam Mode flag true for inhomogeneous computation
             */
            template<int N, int D, bool MODE>
            struct HelperA {
                template<class Potential>
                static void apply(int i, const Potential &V, HaWpParamSet<D> &params, const real_t &delta_t) {
                    // inefficient since all levels are evaluated for each q
                    const auto& leading_level_taylor = V.get_leading_level().taylor_at(complex_t(1,0) * Squeeze<D,RVector<D>>::apply(params.q()));
                    params.updatep( -delta_t * Unsqueeze<D,RVector<D>>::apply(HelperL<MODE>::apply(std::get<1>(leading_level_taylor),i).real()));
                    params.updateP( -delta_t * HelperL<MODE>::apply(std::get<2>(leading_level_taylor),i) * params.Q() );
                    params.updateS( -delta_t * HelperL<MODE>::apply(std::get<0>(leading_level_taylor),i));
                }
            };

            template<int N, int D>
            struct Step2 {
                template<class Potential>
                static void inhomogeneous(int i, const Potential &V, HaWpParamSet<D> &params, const real_t &delta_t) {
                    HelperA<N,D,true>::apply(i,V,params,delta_t);
                }
                template<class Potential>
                static void homogeneous(const Potential &V, HaWpParamSet<D> &params, const real_t &delta_t) {
                    HelperA<N,D,false>::apply(-1,V,params,delta_t);
                }
            };


            /**
             * \brief Builds the inner product matrix
             */
            template<class Packet, class Potential, class IP, int N, int D>
            struct HelperF {
                static void build(CMatrix<Eigen::Dynamic, Eigen::Dynamic>& F,
                                  const Packet& packet,
                                  const Potential& V) {
                    int size = 0;
                    for (auto& component : packet.components()) {
                        size += component.coefficients().size();
                    }
                    F.resize(size,size);

                    // Operator
                    auto op =
                        [&V] (const CMatrix<D,Eigen::Dynamic> &nodes,
                              const RMatrix<D,1> &pos,
                              dim_t i,
                              dim_t j)
                        {
                            const dim_t n_nodes = nodes.cols();
                            CMatrix<1,Eigen::Dynamic> result(1, n_nodes);

                            #pragma omp parallel for schedule(guided)
                            for(int l = 0; l < n_nodes; ++l) {
                                // SUUUUPER INEFFICIENT. COMPUTE N x N Matrix when we only want one entry ...
                                result(0, l) = V.evaluate_local_remainder_at(Squeeze<D, CMatrix<D,Eigen::Dynamic>>::apply(nodes,l),
                                                                             Squeeze<D, CVector<D>>::apply(complex_t(1,0)*pos)) (i,j);
                            }
                            return result;
                        };

                    // Build matrix
                    F = IP::build_matrix(packet, op);
                }
            };

            /**
             * \brief Specialization for N = 1
             */
            template<class Packet, class Potential, class IP,int D>
            struct HelperF<Packet,Potential,IP,1,D> {
                static void build(CMatrix<Eigen::Dynamic, Eigen::Dynamic>& F,
                                  const Packet& packet,
                                  const Potential &V) {

                    // Operator
                    auto op =
                        [&V] (const CMatrix<D,Eigen::Dynamic> &nodes,
                              const RMatrix<D,1> &pos)
                        {
                            const dim_t n_nodes = nodes.cols();
                            CMatrix<1,Eigen::Dynamic> result(1, n_nodes);

                            #pragma omp parallel for schedule(guided)
                            for(int l = 0; l < n_nodes; ++l) {
                                result(0, l) = V.evaluate_local_remainder_at(Squeeze<D, CMatrix<D,Eigen::Dynamic>>::apply(nodes,l),
                                                                             Squeeze<D, CVector<D>>::apply(complex_t(1,0)*pos));
                            }
                            return result;
                        };

                    // Build matrix
                    F = IP::build_matrix(packet, op);
                }
            };

            template<class Packet, class Potential, int N, int D, class IP>
            struct Step3 {

                static void apply(Packet &packet, const Potential& V, const real_t& delta_t) {

                    CMatrix<Eigen::Dynamic,Eigen::Dynamic> F;
                    HelperF<Packet,Potential,IP,N,D>::build(F, packet, V);

                    complex_t factor = -delta_t * complex_t(0,1) / (packet.eps()*packet.eps());

                    // Put all coefficients into a vector
                    CVector<Eigen::Dynamic> coefficients = PacketToCoefficients<Packet>::to(packet);

                    // Compute product of exponential with coefficients
                    coefficients = (factor * F).exp() * coefficients;

                    // Put all coefficients back
                    PacketToCoefficients<Packet>::from(coefficients, packet);
                }
            };

        }
    }
}
