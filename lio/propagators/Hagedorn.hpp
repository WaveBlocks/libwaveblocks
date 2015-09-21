#include <vector>
#include <iostream>
#include "types.hpp"
#include "matrixPotentials/bases.hpp"
#include "matrixPotentials/potentials.hpp"
#include "libwaveblocks/waveblocks/hawp_commons.hpp"
#include "ben/libwaveblocks/waveblocks/inhomogeneous_inner_product.hpp"
#include "ben/libwaveblocks/waveblocks/tensor_product_qr.hpp"
#include "ben/libwaveblocks/waveblocks/gauss_hermite_qr.hpp"



namespace waveblocks
{
  namespace propagators
  {
    template <int N, int D, class MultiIndex>
    struct Hagedorn {
      static void propagate( waveblocks::InhomogeneousHaWp<D,MultiIndex> &packet,
                      const real_t &delta_t,
                      const HomogenousMatrixPotential<N, D> &V,
                      CVector<N> &S
      ) {
        using leading_level_type = typename HomogenousMatrixPotential<N,D>::leading_level_type;
        typename leading_level_type::potential_evaluation_type pot;
        typename leading_level_type::jacobian_evaluation_type jac;
        typename leading_level_type::hessian_evaluation_type hess;
        
        // 1.
        int i = 0;
        for (auto& component : packet.components()) {
          
          auto& params = component.parameters();
          params.q += 0.5 * delta_t * params.p;
          params.Q +=  0.5 * delta_t * params.P;
          S[i] += 0.25 * delta_t * params.p.dot(params.p);
        // 2.
          
          // leading levels
          std::tie (pot,jac,hess) = V.taylor_leading_at(params.q);
          params.p -= delta_t*jac;
          params.P -= delta_t*hess;
          S[i++] -= delta_t*pot;
        }
        
        // 3.
        // Set up quadrature
        using TQR = waveblocks::TensorProductQR < waveblocks::GaussHermiteQR<3>,
              waveblocks::GaussHermiteQR<4>,
              waveblocks::GaussHermiteQR<5> >;
        TQR::NodeMatrix nodes;
        TQR::WeightVector weights;
        waveblocks::InhomogeneousInnerProduct<D, MultiIndex, TQR> ip;

        
        CMatrix<Eigen::Dynamic,Eigen::Dynamic> F;
        int i_offset = 0;
        for (int i = 0; i < N; ++i){
          int i_size = packet.component(i).coefficients().size();
          int j_offset = 0;

          for (int j = 0; j < N; ++j) {
            int j_size = packet.component(j).coefficients().size();

            // Set up operator
            auto op =
                [&V, i, j] (CMatrix<D,Eigen::Dynamic> nodes, CMatrix<D,1> pos) -> CMatrix<1,Eigen::Dynamic>
            {
                const dim_t n_nodes = nodes.cols();
                CMatrix<1,Eigen::Dynamic> result(1, n_nodes);
                for(int l = 0; l < n_nodes; ++l) {
                  result(0, l) = V.evaluate_local_remainder_at(nodes.template block<D,1>(0,l),pos)(i,j); // SUUUUPER INEFFICIENT. COMPUTE NxN Matrix when we only want one entry...
                }
                return result;
            };
            // Build matrix
            F.block(i_offset,j_offset,i_size,j_size) = ip.build_matrix( packet.component(i), packet.component(j), op );
            j_offset += j_size;
          }
          i_offset += i_size;
        }

        
        auto M = -delta_t * ( complex_t( 0, 1 ) / packet.eps() ) * F;
        
        // Exponential
        CMatrix<Eigen::Dynamic,Eigen::Dynamic> expM;
        (Eigen::MatrixExponential<CMatrix<Eigen::Dynamic, Eigen::Dynamic> >( M )).compute( expM );
        int j_offset = 0;
        for(auto& component : packet.components()) {
          int j_size = component.coefficients().size();
          CVector<Eigen::Dynamic> c(component.coefficients().data());
          c = expM.block(0,j_offset,expM.rows(),j_size) * c;
          j_offset += j_size;
        }
        
        
        // 4.
        i = 0;
        for (auto& component : packet.components()) {
          
          auto& params = component.parameters();
          params.q += 0.5 * delta_t * params.p;
          params.Q +=  0.5 * delta_t * params.P;
          S[i++] += 0.25 * delta_t * params.p.dot(params.p);
        }
      }
    };
  }
}
