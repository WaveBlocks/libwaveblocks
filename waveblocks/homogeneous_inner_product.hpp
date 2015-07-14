#pragma once

#include <iostream>
#include <vector>

#include <Eigen/Core>
#include <Eigen/Dense>
#include <unsupported/Eigen/MatrixFunctions>

#include "basic_types.hpp"
#include "hawp.hpp"
#include "basic_types.hpp"

namespace waveblocks {

template<dim_t D, class MultiIndex, class QR>
class HomogeneousInnerProduct
{
public:
    // TODO: make D-dimensional, replace Dynamic
    using CMatrixDD = CMatrix<Eigen::Dynamic, Eigen::Dynamic>;
    using CMatrix1D = CMatrix<1, Eigen::Dynamic>;
    using RMatrix1D = RMatrix<1, Eigen::Dynamic>;

    HomogeneousInnerProduct()
    {
    }

    CMatrixDD build_matrix(const HaWp<D, MultiIndex>& packet, const QR& qr)
        const
    {
        const size_t order = qr.order;
        const real_t eps = packet.basis.eps;
        const CMatrix<D,1>& q = complex_t(1, 0) * packet.basis.parameters->q;
        const CMatrix<D,D>& Q = packet.basis.parameters->Q;
        const CMatrix1D nodes = complex_t(1, 0) *
            RMatrix1D::Map(qr.nodes.data(), qr.nodes.size());
        const CMatrix1D weights = complex_t(1, 0) *
            RMatrix1D::Map(qr.weights.data(), qr.weights.size());

        // Compute affine transformation.
        auto Qs = (Q * Q.adjoint()).inverse().sqrt().inverse();

        // Transform nodes.
        CMatrix1D transformed_nodes =
            q.replicate(1, order) + eps * (Qs * nodes);

        // TODO: Apply operator.
        CMatrix1D values = CMatrix1D::Ones(1, order);

        // TODO: eps^DIM
        Eigen::Array<complex_t, 1, Eigen::Dynamic>
            factor = eps * weights.array() * values.array();
        std::cout << "factor: " << factor << std::endl;

        HaWpBasisVector<Eigen::Dynamic> bases = packet.basis.
            at(transformed_nodes).all();
        //auto bases = packet.basis.
        //    at(transformed_nodes).reduce(
        //            std::vector<complex_t>(order, 0));
        std::cout << "bases:\n" << bases << std::endl;

        // TODO: Placeholder
        return CMatrixDD::Ones(order, order);
    }
};

}
