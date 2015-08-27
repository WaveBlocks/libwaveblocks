#pragma once

#include <iostream>
#include <vector>

#include <Eigen/Core>
#include <Eigen/Dense>

#include "basic_types.hpp"
#include "tables_gausshermite.hpp"

namespace waveblocks {

template <class... RULES>
struct TensorProductQR
{
    static const dim_t D = sizeof...(RULES);

    using NodeMatrix = Eigen::Matrix<real_t,D,Eigen::Dynamic>;
    using WeightVector = Eigen::Matrix<real_t,1,Eigen::Dynamic>;

    static const NodeMatrix nodes()
    {
        std::cout << "Dimensions: " << D << "\n";
        size_t sizes[D] = { RULES::nodes().size() ... };
        std::cout << "size[0]: " << sizes[0] << "\n";

        const dim_t dim = D;
        return NodeMatrix(dim, 1);
    }

    /*
    static const WeightVector weights()
    {
    }
    */
};

}
