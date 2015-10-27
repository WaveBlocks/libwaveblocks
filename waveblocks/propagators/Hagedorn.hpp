#include <vector>
#include <iostream>
#include "../types.hpp"
#include "../matrixPotentials/bases.hpp"
#include "../matrixPotentials/potentials.hpp"
#include "../hawp_commons.hpp"
#include "../inhomogeneous_inner_product.hpp"
#include "../homogeneous_inner_product.hpp"
#include "../tensor_product_qr.hpp"
#include "../gauss_hermite_qr.hpp"
#include "../utilities/adaptors.hpp"
#include "../utilities/squeeze.hpp"


namespace waveblocks
{
  namespace propagators {
    namespace helper {
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
        static void build(CMatrix<Eigen::Dynamic, Eigen::Dynamic>& F, const Packet& packet, const Potential& V) {
          IP ip;
            int size = 0;
            for (auto& component : packet.components()) {
              size += component.coefficients().size();
            }
            F.resize(size,size);
            
            int i_offset = 0;
            #pragma omp parallel for schedule(guided)
            for (int i = 0; i < N; ++i){
              int i_size = packet.component(i).coefficients().size();
              int j_offset = 0;
              
              #pragma omp parallel for schedule(guided)
              for (int j = 0; j < N; ++j) {
                int j_size = packet.component(j).coefficients().size();
                // Set up operator
                auto op =
                    [&V, i, j] (const CMatrix<D,Eigen::Dynamic> &nodes, const RMatrix<D,1> &pos)
                {
                    const dim_t n_nodes = nodes.cols();
                    CMatrix<1,Eigen::Dynamic> result(1, n_nodes);
                    
                    #pragma omp parallel for schedule(guided)
                    for(int l = 0; l < n_nodes; ++l) {
                      result(0, l) = V.evaluate_local_remainder_at(Squeeze<D, CMatrix<D,Eigen::Dynamic>>::apply(nodes,l), Squeeze<D, CVector<D>>::apply(complex_t(1, 0) * pos))(i,j); // SUUUUPER INEFFICIENT. COMPUTE NxN Matrix when we only want one entry...
                    }
                    return result;
                };

                // Build matrix
                F.block(i_offset,j_offset,i_size,j_size) = ip.build_matrix( packet.component(i), packet.component(j), op );
                j_offset += j_size;
              }
              i_offset += i_size;
            }
        }
      };

    /**
      * \brief Specialization for N = 1
      */ 
      template<class Packet, class Potential, class IP,int D>
      struct HelperF<Packet,Potential,IP,1,D> {
        static void build(CMatrix<Eigen::Dynamic, Eigen::Dynamic>& F, const Packet& packet, const Potential &V) {
          IP ip;
          auto op =
              [&V] (const CMatrix<D,Eigen::Dynamic> &nodes, const CMatrix<D,1> &pos) {
              const dim_t n_nodes = nodes.cols();
              CMatrix<1,Eigen::Dynamic> result(n_nodes);

              #pragma omp parallel for schedule(guided)
              for(int l = 0; l < n_nodes; ++l) {
                result(0, l) = V.evaluate_local_remainder_at(Squeeze<D, CMatrix<D,Eigen::Dynamic>>::apply(nodes,l), Squeeze<D, CVector<D>>::apply(pos));
              }
              return result;
          };

          // Build matrix
          F = ip.build_matrix(packet, op);
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

      
  
    using helper::Step1;
    using helper::Step2;
    using helper::Step3;

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
        Step3<InhomogeneousHaWp<D,MultiIndex>,Potential,N,D,InhomogeneousInnerProduct<D,MultiIndex,TQR>>::apply(packet, V, delta_t);
        
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
        Step3<HomogeneousHaWp<D,MultiIndex>,Potential,N,D,InhomogeneousInnerProduct<D,MultiIndex,TQR>>::apply(packet, V, delta_t);
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



