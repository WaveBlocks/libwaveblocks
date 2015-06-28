#pragma once

#include <iostream>
#include <vector>

#include <Eigen/Core>
#include <Eigen/Dense>

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
        // Placeholder
        std::cout << (*packet.coefficients)[0] << std::endl;
        return CMatrixDD::Ones(qr.order, qr.order);
    }
};

}
