#pragma once

#include <vector>

#include "basic_types.hpp"
#include "tables_gausshermite.hpp"

namespace waveblocks {

template <dim_t ORDER>
struct GaussHermiteQR
{
    static const dim_t D = 1;
    static const dim_t order = ORDER;

    using NodeMatrix = Eigen::Matrix<real_t,1,Eigen::Dynamic>;
    using WeightVector = Eigen::Matrix<real_t,1,Eigen::Dynamic>;

    static dim_t number_nodes()
    {
        return order;
    }

    static const std::vector<real_t>& nodes()
    {
        return gauss_hermite_rules[ORDER-1].nodes;
    }

    static const std::vector<real_t>& weights()
    {
        return gauss_hermite_rules[ORDER-1].weights;
    }

    static std::tuple<NodeMatrix,WeightVector> nodes_and_weights()
    {
        return std::make_tuple(
                NodeMatrix::Map(nodes().data(), nodes().size()),
                WeightVector::Map(weights().data(), weights().size()));
    }
};

}
