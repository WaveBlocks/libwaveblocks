#pragma once

#include <iostream>
#include <tuple>
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

    static dim_t number_nodes()
    {
        dim_t result = 1;
        for(auto n_nodes : { RULES::number_nodes() ... })
        {
            result *= n_nodes;
        }
        return result;
    }

    static std::tuple<const NodeMatrix,const WeightVector> nodes_and_weights()
    {
        const dim_t dim = D;
        const dim_t n_nodes = number_nodes();
        NodeMatrix nodes(dim, n_nodes);
        WeightVector weights = WeightVector::Ones(1, n_nodes);

        const dim_t sizes[D] = { RULES::number_nodes() ...};
        dim_t cycles[D];
        cycles[D-1] = 1;
        for(dim_t d = D-2; d >= 0; --d)
        {
            cycles[d] = cycles[d+1] * sizes[d+1];
        }

        // Fill in nodes and weights.
        {
            dim_t d = 0;
            for(auto nw : { RULES::nodes_and_weights() ... })
            {
                for(dim_t n = 0; n < n_nodes; ++n)
                {
                    const dim_t base_idx = (n / cycles[d]) % sizes[d];
                    nodes(d, n) = std::get<0>(nw)(base_idx);
                    weights(n) *= std::get<1>(nw)(base_idx);
                }

                ++d;
            }
        }

        return std::make_tuple(nodes, weights);
    }
};

}
