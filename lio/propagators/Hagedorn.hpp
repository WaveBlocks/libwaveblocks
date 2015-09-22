#include <vector>
#include <iostream>
#include "types.hpp"
#include "matrixPotentials/bases.hpp"
#include "matrixPotentials/potentials.hpp"
#include "hawp_commons.hpp"
#include "inhomogeneous_inner_product.hpp"
#include "tensor_product_qr.hpp"
#include "gauss_hermite_qr.hpp"



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
          std::tie (pot,jac,hess) = V.taylor_leading_at(complex_t(1,0) * params.q); // QUADRATIC REMAINDER
          params.p -= delta_t*jac.real();
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
                [&V, i, j] (const CMatrix<D,Eigen::Dynamic> &nodes, const RMatrix<D,1> &pos) // Do you want real or complex positions??? (You want real...)
            {
                const dim_t n_nodes = nodes.cols();
                CMatrix<1,Eigen::Dynamic> result(1, n_nodes);
                for(int l = 0; l < n_nodes; ++l) {
                  result(0, l) = V.evaluate_local_remainder_at(nodes.template block<D, 1>(0, l), complex_t(1, 0) * pos)(i,j); // SUUUUPER INEFFICIENT. COMPUTE NxN Matrix when we only want one entry...
                }
                return result;
            };

            // Build matrix
            auto comp = packet.component(i);
            auto compj = packet.component(j);
            auto resu = ip.build_matrix(comp,compj,op);
            F.block(i_offset,j_offset,i_size,j_size) = ip.build_matrix( packet.component(i), packet.component(j), op );
            ip.build_matrix( packet.component(i), packet.component(j), op);
            j_offset += j_size;
          }
          i_offset += i_size;
        }

        auto M = -delta_t * ( complex_t( 0, 1 ) / packet.eps() ) * F;
        
        // Exponential
        CMatrix<Eigen::Dynamic,Eigen::Dynamic> expM;
        (Eigen::MatrixExponential<CMatrix<Eigen::Dynamic, Eigen::Dynamic> >( M )).compute( expM );
        
        // Put all coefficients into a vector
        CVector<Eigen::Dynamic> coefficients;
        coefficients.resize(size);
        int j_offset = 0;
        for(auto& component : packet.components()) {
          int j_size = component.coefficients().size();
            for (int j = 0; j < j_size; ++j) {
              coefficients[j+j_offset] = component.coefficients()[j];
            }
          j_offset += j_size;
        }
        
        // Compute product
        coefficients = expM * coefficients;
        
        j_offset = 0;
        // Unpack coefficients
        for(auto& component: packet.components()) {
          int j_size = component.coefficients().size();
            for (int j = 0; j < j_size; ++j) {
              component.coefficients()[j] = coefficients[j+j_offset];
            }
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
