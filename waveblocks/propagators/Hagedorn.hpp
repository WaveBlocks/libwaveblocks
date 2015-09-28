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



namespace waveblocks
{
  template<int N, int D>
  struct Step1 {
    static void apply(HaWpParamSet<D>& params,
                        const real_t &delta_t,
                        complex_t &S) {
      params.q += 0.5 * delta_t * params.p;
      params.Q +=  0.5 * delta_t * params.P;
      S += 0.25 * delta_t * params.p.dot(params.p);
    }
  };
  template<int N>
  struct Step1<N,1> {
    static void apply(HaWpParamSet<1>& params,
                        const real_t &delta_t,
                        complex_t &S) {
      params.q[0] += 0.5 * delta_t * params.p[0];
      params.Q[0] +=  0.5 * delta_t * params.P[0];
      S += 0.25 * delta_t * params.p.dot(params.p);
    }
  };

  template<int N, int D>
  struct Step2 {
    static void inhomogenous(int i, const InhomogenousMatrixPotential<N, D> &V, HaWpParamSet<D> &params, const real_t &delta_t, complex_t &S) {
      const auto& leading_level_taylor = V.get_leading_level().taylor_at(complex_t(1,0) * params.q); // inefficient since all levels are evaluated for each q
          params.p -= delta_t*std::get<1>(leading_level_taylor)[i].real();
          params.P -= delta_t*std::get<2>(leading_level_taylor)[i];
          S -= delta_t*std::get<0>(leading_level_taylor)[i];
    }
    static void homogenous(const HomogenousMatrixPotential<N, D> &V, HaWpParamSet<D> &params, const real_t &delta_t, complex_t &S) {
      const auto& leading_level_taylor = V.get_leading_level().taylor_at(complex_t(1,0) * params.q);
          params.p -= delta_t*std::get<1>(leading_level_taylor).real();
          params.P -= delta_t*std::get<2>(leading_level_taylor);
          S -= delta_t*std::get<0>(leading_level_taylor);
    }
    
  };

  template<int N>
  struct Step2<N,1> {
    static void inhomogenous(int i, const InhomogenousMatrixPotential<N, 1> &V, HaWpParamSet<1> &params, const real_t &delta_t, complex_t &S) {
      const auto& leading_level_taylor = V.get_leading_level().taylor_at(complex_t(1,0) * params.q[0]); // inefficient since all levels are evaluated for each q
          params.p[0] -= delta_t*std::get<1>(leading_level_taylor)[i].real();
          params.P[0] -= delta_t*std::get<2>(leading_level_taylor)[i];
          S -= delta_t*std::get<0>(leading_level_taylor)[i];
    }
    static void homogenous(const HomogenousMatrixPotential<N, 1> &V, HaWpParamSet<1> &params, const real_t &delta_t, complex_t &S) {
      const auto& leading_level_taylor = V.get_leading_level().taylor_at(complex_t(1,0) * params.q[0]);
          params.p[0] -= delta_t*std::get<1>(leading_level_taylor).real();
          params.P[0] -= delta_t*std::get<2>(leading_level_taylor);
          S -= delta_t*std::get<0>(leading_level_taylor);
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

  template<class Packet, class Potential, class IP, int N, int D, class MultiIndex, class TQR>
  struct HelperF {
    static void build(CMatrix<Eigen::Dynamic, Eigen::Dynamic>& F, const Packet& packet, const IP &ip, const Potential& V) {
        int size = 0;
        for (auto& component : packet.components()) {
          size += component.coefficients().size();
        }
        F.resize(size,size);
        
        int i_offset = 0;
        for (int i = 0; i < N; ++i){
          int i_size = packet.component(i).coefficients().size();
          int j_offset = 0;

          for (int j = 0; j < N; ++j) {
            int j_size = packet.component(j).coefficients().size();
            // Set up operator
            auto op =
                [&V, i, j] (const CMatrix<D,Eigen::Dynamic> &nodes, const RMatrix<D,1> &pos)
            {
                const dim_t n_nodes = nodes.cols();
                CMatrix<1,Eigen::Dynamic> result(1, n_nodes);
                for(int l = 0; l < n_nodes; ++l) {
                  result(0, l) = V.evaluate_local_remainder_at(HelperArg<D>::first(nodes,l), HelperArg<D>::second(complex_t(1, 0) * pos))(i,j); // SUUUUPER INEFFICIENT. COMPUTE NxN Matrix when we only want one entry...
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

  template<class Packet, class Potential, class IP,int D, class MultiIndex, class TQR>
  struct HelperF<Packet,Potential,IP,1,D,MultiIndex,TQR> {
    static void build(CMatrix<Eigen::Dynamic, Eigen::Dynamic>& F, const Packet& packet, const IP &ip, const Potential &V) {
      auto op =
          [&V] (const CMatrix<D,Eigen::Dynamic> &nodes, const RMatrix<D,1> &pos) {
          const dim_t n_nodes = nodes.cols();
          CMatrix<1,Eigen::Dynamic> result(n_nodes);
          for(int l = 0; l < n_nodes; ++l) {
            result(0, l) = V.evaluate_local_remainder_at(HelperArg<D>::first(nodes,l), HelperArg<D>::second(complex_t(1,0)* pos));
          }
          return result;
      };

      // Build matrix
      F = ip.build_matrix(packet, packet, op);
    }
  };
    
  template<class Packet, class Potential, int N, int D, class MultiIndex, class TQR>
  struct Step3 {
    using IP = InhomogeneousInnerProduct<D, MultiIndex, TQR>;

    static void apply(Packet &packet, const Potential& V, const real_t& delta_t) {
      InhomogeneousInnerProduct<D, MultiIndex, TQR> ip;

      CMatrix<Eigen::Dynamic,Eigen::Dynamic> F;
      HelperF<Packet,Potential,InhomogeneousInnerProduct<D, MultiIndex, TQR>,N,D,MultiIndex,TQR>::build(F, packet, ip, V);

      auto M = -delta_t * ( complex_t( 0, 1 ) / packet.eps() ) * F;
      
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

  template<int N, int D>
  struct Step4 {
    static void apply(HaWpParamSet<D>& params,
                        const real_t &delta_t,
                        complex_t &S) {
      params.q += 0.5 * delta_t * params.p;
      params.Q +=  0.5 * delta_t * params.P;
      S += 0.25 * delta_t * params.p.dot(params.p);}
  };

  template<int N>
  struct Step4<N,1> {
    static void apply(HaWpParamSet<1>& params,
                        const real_t &delta_t,
                        complex_t &S){
      params.q[0] += 0.5 * delta_t * params.p[0];
      params.Q[0] +=  0.5 * delta_t * params.P[0];
      S += 0.25 * delta_t * params.p.dot(params.p);}
  };
    
  
  namespace propagators
  {
    template <int N, int D, class MultiIndex, class TQR>
    struct Hagedorn {
      static void propagate(InhomogeneousHaWp<D,MultiIndex> &packet,
                      const real_t &delta_t,
                      const InhomogenousMatrixPotential<N, D> &V,
                      CVector<N> &S
      ) {
        
        int i = 0;
        for (auto& component : packet.components()) {

          auto& params = component.parameters();
          Step1<N,D>::apply(params,delta_t,S[i]);
          Step2<N,D>::inhomogenous(i,V,params,delta_t,S[i]);
          i++;
        }
        Step3<InhomogeneousHaWp<D,MultiIndex>,InhomogenousMatrixPotential<N,D>,N,D,MultiIndex,TQR>::apply(packet, V, delta_t);
        
        i = 0;
        for (auto& component : packet.components()) {
          auto& params = component.parameters();
          Step4<N,D>::apply(params,delta_t,S[i]);
          i++;
        }
      }

      static void propagate( HomogeneousHaWp<D,MultiIndex> &packet,
                      const real_t &delta_t,
                      const HomogenousMatrixPotential<N, D> &V,
                      complex_t &S
      ) {
        auto& params = packet.parameters();
        Step1<N,D>::apply(params, delta_t, S);
        Step2<N,D>::homogenous(V,params,delta_t,S);
        Step3<HomogeneousHaWp<D,MultiIndex>,HomogenousMatrixPotential<N,D>,N,D,MultiIndex,TQR>::apply(packet, V, delta_t);
        Step4<N,D>::apply(params,delta_t,S);

        }
    };
    template<int D, class MultiIndex, class TQR>
    struct Hagedorn<1,D,MultiIndex,TQR> {
      static void propagate(ScalarHaWp<D, MultiIndex> &packet,
        const real_t &delta_t,
        const ScalarMatrixPotential<D> &V,
        complex_t &S) {
          auto& params = packet.parameters();
          Step1<1,D>::apply(params, delta_t, S);
          Step2<1,D>::homogenous(V,params,delta_t,S);
          Step3<ScalarHaWp<D,MultiIndex>,ScalarMatrixPotential<D>,1,D,MultiIndex,TQR>::apply(packet, V, delta_t);
          Step4<1,D>::apply(params,delta_t,S);
        }
    };
  }
}



