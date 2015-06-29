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
    using CMatrixDD = CMatrix<Eigen::Dynamic, Eigen::Dynamic>;

    HomogeneousInnerProduct()
    {
    }

    CMatrixDD build_matrix(const HaWp<D, MultiIndex>& packet, const QR& qr)
        const
    {
        const RMatrix<D,1>& q = packet.basis.parameters->q;
        const CMatrix<D,D>& Q = packet.basis.parameters->Q;

        // Transform nodes.
        std::cout << "Q: " << Q << std::endl;
        auto Q0 = (Q * Q.adjoint()).inverse();
        std::cout << "Q0: " << Q0 << std::endl;
        auto Qs = Q0.sqrt().inverse();
        std::cout << "Qs: " << Qs << std::endl;

        // Placeholder
        return CMatrixDD::Ones(qr.order, qr.order);
    }
};

}
