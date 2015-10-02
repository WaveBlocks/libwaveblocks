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
    using RMatrixD1 = RMatrix<D, 1>;
    using NodeMatrix = typename QR::NodeMatrix;
    using WeightVector = typename QR::WeightVector;
    using op_t = std::function<CMatrix1N(CMatrixDN,RMatrixD1)>;

    InhomogeneousInnerProduct()
    {
    }

    CMatrixNN build_matrix(const AbstractScalarHaWp<D, MultiIndex>& pacbra,
            const AbstractScalarHaWp<D, MultiIndex>& packet,
            const op_t& op=default_op)
        const
    {
        const dim_t n_nodes = QR::number_nodes();
        const CMatrixD1& qr = complex_t(1, 0) * pacbra.parameters().q;
        const CMatrixD1& qc = complex_t(1, 0) * packet.parameters().q;
        const CMatrixDD& Qr = pacbra.parameters().Q;
        const CMatrixDD& Qc = packet.parameters().Q;
        const CMatrixDD& Pr = pacbra.parameters().P;
        const CMatrixDD& Pc = packet.parameters().P;
        const complex_t S_bra = pacbra.parameters().S;
        const complex_t S_ket = packet.parameters().S;
        NodeMatrix nodes;
        WeightVector weights;
        std::tie(nodes, weights) = QR::nodes_and_weights();
        const CMatrixDN cnodes = complex_t(1, 0) * nodes;
        const CMatrix1N cweights = complex_t(1, 0) * weights;

        // Mix parameters, compute affine transformation.
        CMatrix<D,D> Gr = Pr * Qr.inverse();
        CMatrix<D,D> Gc = Pc * Qc.inverse();
        RMatrix<D,D> r = (Gc - Gr.adjoint()).imag();
        RMatrix<D,1> s = ((Gc * qc) - (Gr.adjoint() * qr)).imag();
        RMatrix<D,1> q0 = r.inverse() * s;
        RMatrix<D,D> Q0 = 0.5 * r;
        RMatrix<D,D> Qs = Q0.sqrt().inverse();

        
        //~ auto PIket = packet.parameters();
        //~ auto PIbra = pacbra.parameters();
//~ 
        //~ auto PImix = PIbra.mix(PIket);
//~ 
        //~ std::cout << "qmix: " << std::get<0>(PImix) << std::endl;
        //~ std::cout << "Qmix: " << std::get<1>(PImix) << std::endl;
//~ 
        //~ std::cout << "q0? [" << q0 << "]"<<std::endl << std::endl;
        //~ std::cout << "Qs? [" << Qs << "]"<<std::endl << std::endl;
//~ 
        //~ std::cout << "Qr " << Qr << std::endl;
        // Transform nodes.
        CMatrixDN transformed_nodes = complex_t(1, 0) *
            q0.replicate(1, n_nodes) + packet.eps() * (Qs * cnodes);

        // Apply operator.
        CMatrix1N values = op(transformed_nodes, q0);

        Eigen::Array<complex_t, 1, Eigen::Dynamic> factor =
            std::pow(packet.eps(), D) * cweights.array() * values.array() *
              Qs.determinant() * pacbra.prefactor() * packet.prefactor();
              
        HaWpBasisVector<Eigen::Dynamic> basisr =
            pacbra.evaluate_basis(transformed_nodes);
        HaWpBasisVector<Eigen::Dynamic> basisc =
            packet.evaluate_basis(transformed_nodes);
        //std::cout << "bases(:,0):\n" << bases.col(0) << std::endl;

        //std::cout << "factor: " << factor.rows() << " x " << factor.cols() << "\n";
        //std::cout << "bases: " << bases.rows() << " x " << bases.cols() << "\n";

        // Build matrix.
        CMatrixNN result = CMatrixNN::Zero(basisr.rows(), basisc.rows());
        for(dim_t i = 0; i < basisr.rows(); ++i)
        {
            for(dim_t j = 0; j < basisc.rows(); ++j)
            {
                for(dim_t k = 0; k < n_nodes; ++k)
                {
                    result(i, j) += factor(k) * conj(basisr(i, k)) * basisc(j, k);
                }
            }
        }

        // TODO: Phase calculation ("S" parameter?)
        auto phase = std::exp(
          complex_t(0,1) * (S_ket - std::conj(S_bra))/std::pow(packet.eps(),2));
          
        return phase * result;
    }

    complex_t quadrature(const AbstractScalarHaWp<D, MultiIndex>& pacbra,
            const AbstractScalarHaWp<D, MultiIndex>& packet,
            const op_t& op=default_op)
        const
    {
        const auto M = build_matrix(pacbra, packet, op);
        // Quadrature with wavepacket coefficients, c^H M c.
        const CMatrixN1 coeffs_bra = CMatrixN1::Map(
                pacbra.coefficients().data(), pacbra.coefficients().size());
        const CMatrixN1 coeffs_ket = CMatrixN1::Map(
                packet.coefficients().data(), packet.coefficients().size());
        //std::cout << "\nM: " << M.rows() << " x " << M.cols() << "\n";
        //std::cout << "c: " << coeffs.rows() << " x " << coeffs.cols() << "\n";
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
