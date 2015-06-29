#pragma once

#include <iostream>
#include <vector>

#include <Eigen/Core>
#include <Eigen/Dense>

#include "basic_types.hpp"
#include "tables_gausshermite.hpp"

namespace waveblocks {

template <dim_t ORDER>
struct GaussHermiteQR
{
    static const dim_t order = ORDER;
    const std::vector<real_t> nodes, weights;

    GaussHermiteQR()
        : nodes(gauss_hermite_rules[ORDER-1].nodes)
        , weights(gauss_hermite_rules[ORDER-1].weights)
    {
    }
};

}
