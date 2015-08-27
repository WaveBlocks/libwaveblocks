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
        WeightVector weights(1, n_nodes);

        return std::make_tuple(nodes, weights);
    }

private:
    // returns next node index
    /*
    static dim_t fill_nodes(dim_t depth, dim_t idx)
    {
        for()
        {
        }
    }
    */
};

}
