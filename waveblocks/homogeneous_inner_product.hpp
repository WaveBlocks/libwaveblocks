#pragma once

#include <cmath>
#include <iostream>
#include <tuple>
#include <vector>

#include <Eigen/Core>
#include <Eigen/Dense>
#include <unsupported/Eigen/MatrixFunctions>

#include "basic_types.hpp"
#include "hawp.hpp"
#include "hawp_commons.hpp"

namespace waveblocks {

template<dim_t D, class MultiIndex, class QR>
class HomogeneousInnerProduct
{
public:
    // TODO: make D-dimensional, replace Dynamic
    using CMatrixNN = CMatrix<Eigen::Dynamic, Eigen::Dynamic>;
    using CMatrix1N = CMatrix<1, Eigen::Dynamic>;
    using CMatrixDN = CMatrix<D, Eigen::Dynamic>;
    using NodeMatrix = typename QR::NodeMatrix;
    using WeightVector = typename QR::WeightVector;

    HomogeneousInnerProduct()
    {
    }

    CMatrixNN build_matrix(const AbstractScalarHaWpBasis<D, MultiIndex>& packet)
        const
    {
        const dim_t n_nodes = QR::number_nodes();
        const CMatrix<D,1>& q = complex_t(1, 0) * packet.parameters().q;
        const CMatrix<D,D>& Q = packet.parameters().Q;
        NodeMatrix nodes;
        WeightVector weights;
        std::tie(nodes, weights) = QR::nodes_and_weights();
        const CMatrixDN cnodes = complex_t(1, 0) * nodes;
        const CMatrix1N cweights = complex_t(1, 0) * weights;

        // Compute affine transformation.
        auto Qs = (Q * Q.adjoint()).inverse().sqrt().inverse();

        // Transform nodes.
        CMatrixDN transformed_nodes =
            q.replicate(1, n_nodes) + packet.eps() * (Qs * cnodes);

        // TODO: Apply operator.
        CMatrix1N values = CMatrix1N::Ones(1, n_nodes);

        Eigen::Array<complex_t, 1, Eigen::Dynamic> factor =
            std::pow(packet.eps(), D) * cweights.array() * values.array();
        //std::cout << "factor: " << factor << std::endl;

        HaWpBasisVector<Eigen::Dynamic> bases =
            packet.evaluate_basis(transformed_nodes);
        //std::cout << "bases(:,0):\n" << bases.col(0) << std::endl;

        //std::cout << "factor: " << factor.rows() << " x " << factor.cols() << "\n";
        //std::cout << "bases: " << bases.rows() << " x " << bases.cols() << "\n";

        // Build matrix.
        const dim_t N = bases.rows();
        CMatrixNN result = CMatrixNN::Zero(N, N);
        for(dim_t i = 0; i < N; ++i)
        {
            for(dim_t j = 0; j < N; ++j)
            {
                for(dim_t k = 0; k < n_nodes; ++k)
                {
                    result(i, j) += factor(k) * conj(bases(i, k)) * bases(j, k);
                }
            }
        }

        return result;
    }
};

}
