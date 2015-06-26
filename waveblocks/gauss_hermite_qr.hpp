#pragma once

#include <iostream>
#include <vector>

#include <Eigen/Core>
#include <Eigen/Dense>

#include "basic_types.hpp"

namespace waveblocks {

struct GaussHermiteQR
{
public:
    const dim_t order;
    const std::vector<real_t> nodes, weights;

    GaussHermiteQR(dim_t order)
        : order(order)
        , nodes(calculate_nodes())
        , weights(calculate_weights())
    {
    }

private:
    std::vector<real_t> calculate_nodes() const
    {
        // Placeholder
        return std::vector<real_t>(order, 0);
    }

    std::vector<real_t> calculate_weights() const
    {
        // Placeholder
        return std::vector<real_t>(order, 0);
    }
};

}
