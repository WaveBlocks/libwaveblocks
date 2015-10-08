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
    using CDiagonalNN = Eigen::DiagonalMatrix<complex_t, Eigen::Dynamic>;
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
    CMatrixNN build_matrix(const AbstractScalarHaWp<D, MultiIndex>& packet,
                           const op_t& op=default_op) const {
        const dim_t n_nodes = QR::number_nodes();
        const CMatrixD1& q = packet.parameters().q.template cast<complex_t>();
        const CMatrixDD& Q = packet.parameters().Q;
        NodeMatrix nodes;
        WeightVector weights;
        std::tie(nodes, weights) = QR::nodes_and_weights();

        // Compute affine transformation
        const CMatrixDD Qs = (Q * Q.adjoint()).sqrt();

        // Transform nodes
        const CMatrixDN transformed_nodes = q.replicate(1, n_nodes) + packet.eps() * (Qs * nodes);

        // Apply operator
        const CMatrix1N values = op(transformed_nodes, q);

        // Prefactor
        const CMatrix1N factor =
            // std::conj(packet.prefactor()) * packet.prefactor() * Qs.determinant() = 1
            std::pow(packet.eps(), D) * weights.array() * values.array();

        // Evaluate basis
        const CMatrixNN basis = packet.evaluate_basis(transformed_nodes);

        // Build matrix
        const CDiagonalNN Dfactor(factor);
        const CMatrixNN result = basis.matrix().conjugate() * Dfactor * basis.matrix().transpose();

        #pragma omp parallel for schedule(guided)
                        complex_t resij = 0.0;
                        result(i, j) = resij;

        // Global phase cancels out
        return result;
    }

    complex_t quadrature(const AbstractScalarHaWp<D, MultiIndex>& packet,
                         const op_t& op=default_op) const {
        const auto M = build_matrix(packet, op);
        // Quadrature with wavepacket coefficients, c^H M c.
        const CMatrixN1 coeffs = CMatrixN1::Map(packet.coefficients().data(), packet.coefficients().size());
        return coeffs.adjoint() * M * coeffs;
    }

private:
    static CMatrix1N default_op(const CMatrixDN& nodes, const CMatrixD1& pos)
    {
        (void)pos;
        return CMatrix1N::Ones(1, nodes.cols());
    }
};

}
