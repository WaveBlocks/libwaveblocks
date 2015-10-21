#include <vector>
#include <iostream>
#include "../types.hpp"
#include "../matrixPotentials/bases.hpp"
#include "../matrixPotentials/potentials.hpp"
#include "../hawp_commons.hpp"
#include "../inhomogeneous_inner_product.hpp"
#include "../homogeneous_inner_product.hpp"
#include "../vector_inner_product.hpp"
#include "../tensor_product_qr.hpp"
#include "../gauss_hermite_qr.hpp"
#include "../utilities/adaptors.hpp"


namespace waveblocks
{
    template<int N, int D>
    struct Step1 {
        static void apply(HaWpParamSet<D>& params,
                          const real_t &delta_t) {
            params.updateq( 0.5 * delta_t * params.p() );
            params.updateQ( 0.5 * delta_t * params.P() );
            params.updateS( 0.25 * delta_t * params.p().dot(params.p()) );
        }
    };
    template<int N>
    struct Step1<N,1> {
        static void apply(HaWpParamSet<1>& params,
                          const real_t &delta_t) {
            params.updateq( 0.5 * delta_t * params.p() );
            params.updateQ( 0.5 * delta_t * params.P() );
            params.updateS( 0.25 * delta_t * params.p().dot(params.p()) );
        }
    };

    template<int N, int D>
    struct Step2 {
        template<class Potential>
        static void inhomogenous(int i, const Potential &V, HaWpParamSet<D> &params, const real_t &delta_t) {
            // inefficient since all levels are evaluated for each q
            const auto& leading_level_taylor = V.get_leading_level().taylor_at(complex_t(1,0) * params.q());
            params.updatep( -delta_t * std::get<1>(leading_level_taylor)[i].real() );
            params.updateP( -delta_t * std::get<2>(leading_level_taylor)[i] * params.Q() );
            params.updateS( -delta_t * std::get<0>(leading_level_taylor)[i] );
        }
        template<class Potential>
        static void homogenous(const Potential &V, HaWpParamSet<D> &params, const real_t &delta_t) {
            const auto& leading_level_taylor = V.get_leading_level().taylor_at(complex_t(1,0) * params.q());
            params.updatep( -delta_t * std::get<1>(leading_level_taylor).real() );
            params.updateP( -delta_t * std::get<2>(leading_level_taylor) * params.Q() );
            params.updateS( -delta_t * std::get<0>(leading_level_taylor) );
        }
    };

    template<int N>
    struct Step2<N,1> {
        template<class Potential>
        static void inhomogenous(int i, const Potential &V, HaWpParamSet<1> &params, const real_t &delta_t) {
            // inefficient since all levels are evaluated for each q
            const auto& leading_level_taylor = V.get_leading_level().taylor_at(complex_t(1,0) * params.q()[0]);
            params.updatep( -delta_t * std::get<1>(leading_level_taylor)[i].real() * RMatrix<1,1>::Identity() );
            params.updateP( -delta_t * std::get<2>(leading_level_taylor)[i] * params.Q() );
            params.updateS( -delta_t * std::get<0>(leading_level_taylor)[i] );
        }
        template<class Potential>
        static void homogenous(const Potential &V, HaWpParamSet<1> &params, const real_t &delta_t) {
            const auto& leading_level_taylor = V.get_leading_level().taylor_at(complex_t(1,0) * params.q()[0]);
            params.updatep( -delta_t * std::get<1>(leading_level_taylor).real() * RMatrix<1,1>::Identity() );
            params.updateP( -delta_t * std::get<2>(leading_level_taylor) * params.Q() );
            params.updateS( -delta_t * std::get<0>(leading_level_taylor) );
        }
    };

  template<int D>
  struct HelperArg {
    static CMatrix<D,1> first(const CMatrix<D,Eigen::Dynamic>& nodes, int l) {
      return nodes.template block<D,1>(0,l);
    }
    static const CMatrix<D,1>& second(const CMatrix<D,1> &pos) {
      return pos;
    }
  };

  template<>
  struct HelperArg<1> {
    static const complex_t& first(const CMatrix<1,Eigen::Dynamic>& nodes, int l) {
      return nodes(0,l);
    }
    static const complex_t& second(const CMatrix<1,1> &pos) {
      return pos[0];
    }
  };

  template<class Packet, class Potential, class IP, int N, int D>
  struct HelperF {
    static void build(CMatrix<Eigen::Dynamic, Eigen::Dynamic>& F, const Packet& packet, const Potential& V) {
        int size = 0;
        for (auto& component : packet.components()) {
          size += component.coefficients().size();
        }
        F.resize(size,size);
        
        //int i_offset = 0;
        //for (int i = 0; i < N; ++i){
        //  int i_size = packet.component(i).coefficients().size();
        //  int j_offset = 0;
        //  
        //  for (int j = 0; j < N; ++j) {
        //    int j_size = packet.component(j).coefficients().size();
        //    // Set up operator
        //    auto op =
        //        [&V, i, j] (const CMatrix<D,Eigen::Dynamic> &nodes, const RMatrix<D,1> &pos)
        //    {
        //        const dim_t n_nodes = nodes.cols();
        //        CMatrix<1,Eigen::Dynamic> result(1, n_nodes);
        //        
        //        #pragma omp parallel for schedule(guided)
        //        for(int l = 0; l < n_nodes; ++l) {
        //          result(0, l) = V.evaluate_local_remainder_at(HelperArg<D>::first(nodes,l), HelperArg<D>::second(complex_t(1, 0) * pos))(i,j); // SUUUUPER INEFFICIENT. COMPUTE NxN Matrix when we only want one entry...
        //        }
        //        return result;
        //    };

        //    // Build matrix
        //    F.block(i_offset,j_offset,i_size,j_size) = IP::build_matrix( packet.component(i), packet.component(j), op );
        //    j_offset += j_size;
        //  }
        //  i_offset += i_size;
        //}

        auto op =
            [&V] (const CMatrix<D,Eigen::Dynamic> &nodes, const RMatrix<D,1> &pos, dim_t i, dim_t j)
        {
            const dim_t n_nodes = nodes.cols();
            CMatrix<1,Eigen::Dynamic> result(1, n_nodes);
            
            #pragma omp parallel for schedule(guided)
            for(int l = 0; l < n_nodes; ++l) {
              result(0, l) = V.evaluate_local_remainder_at(HelperArg<D>::first(nodes,l), HelperArg<D>::second(complex_t(1, 0) * pos))(i,j); // SUUUUPER INEFFICIENT. COMPUTE NxN Matrix when we only want one entry...
            }
            return result;
        };
        F = IP::build_matrix(packet, op);
    }
  };

  template<class Packet, class Potential, class IP,int D>
  struct HelperF<Packet,Potential,IP,1,D> {
    static void build(CMatrix<Eigen::Dynamic, Eigen::Dynamic>& F, const Packet& packet, const Potential &V) {
      auto op =
          [&V] (const CMatrix<D,Eigen::Dynamic> &nodes, const CMatrix<D,1> &pos) {
          const dim_t n_nodes = nodes.cols();
          CMatrix<1,Eigen::Dynamic> result(n_nodes);
          for(int l = 0; l < n_nodes; ++l) {
            result(0, l) = V.evaluate_local_remainder_at(HelperArg<D>::first(nodes,l), HelperArg<D>::second(pos));
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

      auto M = -delta_t *  complex_t(0,1) / (packet.eps()*packet.eps()) * F;

      // Exponential
      CMatrix<Eigen::Dynamic,Eigen::Dynamic> expM;
      (Eigen::MatrixExponential<CMatrix<Eigen::Dynamic, Eigen::Dynamic> >( M )).compute( expM );
      
      // Put all coefficients into a vector
      CVector<Eigen::Dynamic> coefficients = utilities::PacketToCoefficients<Packet>::to(packet);

      // Compute product
      coefficients = expM * coefficients;
      
      utilities::PacketToCoefficients<Packet>::from(coefficients, packet);
    }
  };

    
  
  namespace propagators
  {
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
          Step2<N,D>::inhomogenous(i,V,params,delta_t);
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
        Step2<N,D>::homogenous(V,params,delta_t);
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
          Step2<1,D>::homogenous(V,params,delta_t);
          Step3<ScalarHaWp<D,MultiIndex>,Potential,1,D,HomogeneousInnerProduct<D,MultiIndex,TQR>>::apply(packet, V, delta_t);
          Step1<1,D>::apply(params,delta_t);
        }
    };
  }
}



