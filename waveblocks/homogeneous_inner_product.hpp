#pragma once

#include <cmath>
#include <functional>
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

/**
 * \brief Class providing calculation of homogeneous inner products.
 *
 * \tparam D dimensionality of processed wavepackets
 * \tparam MultiIndex multi-index type of processed wavepackets
 * \tparam QR quadrature rule to use, with N nodes
 */
template<dim_t D, class MultiIndex, class QR>
class HomogeneousInnerProduct
{
public:
    using CMatrixNN = CMatrix<Eigen::Dynamic, Eigen::Dynamic>;
    using CMatrix1N = CMatrix<1, Eigen::Dynamic>;
    using CMatrixN1 = CMatrix<Eigen::Dynamic, 1>;
    using CMatrixD1 = CMatrix<D, 1>;
    using CMatrixDD = CMatrix<D, D>;
    using CMatrixDN = CMatrix<D, Eigen::Dynamic>;
    using NodeMatrix = typename QR::NodeMatrix;
    using WeightVector = typename QR::WeightVector;
    using op_t = std::function<CMatrix1N(CMatrixDN,CMatrixD1)>;

    HomogeneousInnerProduct()
    {
    }

    /**
     * \brief Calculate the matrix of the inner product.
     *
     * Returns the matrix elements \f$\langle \phi | f | \phi \rangle\f$ with
     * an operator \f$f\f$.
     * The coefficients of the wavepacket are ignored.
     *
     * \param packet wavepacket \f$\phi\f$
     * \param op operator \f$f(x, q) : \mathbb{C}^{D \times N} \times
     *   \mathbb{C}^D \rightarrow \mathbb{C}^N\f$ which is evaluated at the
     *   nodal points \f$x\f$ and position \f$q\f$;
     *   default returns a vector of ones
     */
    CMatrixNN build_matrix(const AbstractScalarHaWpBasis<D, MultiIndex>& packet,
            const op_t& op=default_op)
        const
    {
        const dim_t n_nodes = QR::number_nodes();
        const CMatrixD1& q = complex_t(1, 0) * packet.parameters().q;
        const CMatrixDD& Q = packet.parameters().Q;
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

        // Apply operator.
        CMatrix1N values = op(transformed_nodes, q);

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

    complex_t quadrature(const AbstractScalarHaWp<D, MultiIndex>& packet,
            const op_t& op=default_op)
        const
    {
        const auto M = build_matrix(packet, op);
        // Quadrature with wavepacket coefficients, c^H M c.
        const CMatrixN1 coeffs = CMatrixN1::Map(
                packet.coefficients().data(), packet.coefficients().size());
        //std::cout << "\nM: " << M.rows() << " x " << M.cols() << "\n";
        //std::cout << "c: " << coeffs.rows() << " x " << coeffs.cols() << "\n";
        return coeffs.adjoint() * M * coeffs;
    }

private:
    static CMatrix1N default_op(const CMatrixDN& nodes, const CMatrixD1& pos)
    {
        return CMatrix1N::Ones(1, nodes.cols());
    }
};

}
