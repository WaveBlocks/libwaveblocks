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

template<dim_t D, class MultiIndex, class QR>
class InhomogeneousInnerProduct
{
public:
    using CMatrixNN = CMatrix<Eigen::Dynamic, Eigen::Dynamic>;
    using CMatrix1N = CMatrix<1, Eigen::Dynamic>;
    using CMatrixN1 = CMatrix<Eigen::Dynamic, 1>;
    using CMatrixD1 = CMatrix<D, 1>;
    using CMatrixDD = CMatrix<D, D>;
    using CMatrixDN = CMatrix<D, Eigen::Dynamic>;
    using RMatrixDD = RMatrix<D, D>;
    using RMatrixD1 = RMatrix<D, 1>;
    using CDiagonalNN = Eigen::DiagonalMatrix<complex_t, Eigen::Dynamic>;
    using NodeMatrix = typename QR::NodeMatrix;
    using WeightVector = typename QR::WeightVector;
    using op_t = std::function<CMatrix1N(CMatrixDN,RMatrixD1)>;

    InhomogeneousInnerProduct()
    {
    }

    CMatrixNN build_matrix(const AbstractScalarHaWp<D, MultiIndex>& pacbra,
                           const AbstractScalarHaWp<D, MultiIndex>& packet,
                           const op_t& op=default_op) const {
        const dim_t n_nodes = QR::number_nodes();
        const complex_t S_bra = pacbra.parameters().S;
        const complex_t S_ket = packet.parameters().S;
        NodeMatrix nodes;
        WeightVector weights;
        std::tie(nodes, weights) = QR::nodes_and_weights();

        // Mix parameters and compute affine transformation
        std::pair<RMatrixD1, RMatrixDD> PImix = pacbra.parameters().mix(packet.parameters());
        const RMatrixD1& q0 = std::get<0>(PImix);
        const RMatrixDD& Qs = std::get<1>(PImix);

        // Transform nodes
        const CMatrixDN transformed_nodes = q0.template cast<complex_t>().replicate(1, n_nodes) + packet.eps() * (Qs.template cast<complex_t>() * nodes);

        // Apply operator
        const CMatrix1N values = op(transformed_nodes, q0);

        // Prefactor
        const CMatrix1N factor =
            std::conj(pacbra.prefactor()) * packet.prefactor() * Qs.determinant() *
            std::pow(packet.eps(), D) * weights.array() * values.array();

        // Evaluate basis
        const CMatrixNN basisr = pacbra.evaluate_basis(transformed_nodes);
        const CMatrixNN basisc = packet.evaluate_basis(transformed_nodes);

        // Build matrix
        const CDiagonalNN Dfactor(factor);
        const CMatrixNN result = basisr.matrix().conjugate() * Dfactor * basisc.matrix().transpose();

        // Global phase
        const complex_t phase = std::exp(complex_t(0,1) * (S_ket - std::conj(S_bra)) / std::pow(packet.eps(),2));
        return phase * result;
    }

    complex_t quadrature(const AbstractScalarHaWp<D, MultiIndex>& pacbra,
                         const AbstractScalarHaWp<D, MultiIndex>& packet,
                         const op_t& op=default_op) const {
        const auto M = build_matrix(pacbra, packet, op);
        // Quadrature with wavepacket coefficients, c^H M c.
        const CMatrixN1 coeffs_bra = CMatrixN1::Map(pacbra.coefficients().data(), pacbra.coefficients().size());
        const CMatrixN1 coeffs_ket = CMatrixN1::Map(packet.coefficients().data(), packet.coefficients().size());
        return coeffs_bra.adjoint() * M * coeffs_ket;
    }

private:
    static CMatrix1N default_op(const CMatrixDN& nodes, const RMatrixD1& pos)
    {
        (void)pos;
        return CMatrix1N::Ones(1, nodes.cols());
    }
};

}
