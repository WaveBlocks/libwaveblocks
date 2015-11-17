#pragma once

#include <iostream>
#include <tuple>
#include <vector>

#include <Eigen/Core>
#include <Eigen/Dense>

#include "basic_types.hpp"
#include "tables_gausshermite.hpp"
#include "util/combinatorics.hpp"

namespace waveblocks {


/**
 * \brief Structure providing weighted nodes for Genz-Keister quadrature.
 */
template <dim_t D, dim_t LEVEL>
struct GenzKeisterQR
{
    using NodeMatrix = Eigen::Matrix<real_t,D,Eigen::Dynamic>;
    using WeightVector = Eigen::Matrix<real_t,1,Eigen::Dynamic>;


    /**
     * \brief Return the number of nodes for the given order.
     */
    static dim_t number_nodes()
    {
        if (!cached) calculate_nodes_and_weights();
        return cached_nodes.cols();
    }

    /**
     * \brief Return the quadrature nodes.
     */
    static NodeMatrix nodes()
    {
        if (!cached) calculate_nodes_and_weights();
        return cached_nodes;
    }


    /**
     * \brief Return the quadrature weights.
     */
    static WeightVector weights()
    {
        if (!cached) calculate_nodes_and_weights();
        return cached_weights;
    }

    /**
     * \brief Return the quadrature nodes and weights.
     */
    static std::tuple<NodeMatrix,WeightVector>
        nodes_and_weights()
    {
        if (!cached) calculate_nodes_and_weights();
        return std::make_tuple(cached_nodes, cached_weights);
    }


    /**
     * \brief Free the precalculated nodes and weights.
     */
    static void clear_cache()
    {
        cached = false;
        cached_nodes.resize(D, 0);
        cached_weights.resize(1, 0);
    }

private:
    static bool cached;
    static NodeMatrix cached_nodes;
    static WeightVector cached_weights;

    /**
     * \brief Precalculate nodes and weights.
     */
    static void calculate_nodes_and_weights()
    {
        clear_cache();

        const dim_t K = LEVEL - 1;
        const auto parts = partitions<D>(K);
        std::cout << "Partitions:\n";
        for (const auto& part : parts)
        {
            std::cout << part << "\n";
        }

        cached = true;
    }
};

// Initialize static members.
template <dim_t D, dim_t LEVEL>
bool GenzKeisterQR<D, LEVEL>::cached = false;

template <dim_t D, dim_t LEVEL>
typename GenzKeisterQR<D, LEVEL>::NodeMatrix
    GenzKeisterQR<D, LEVEL>::cached_nodes;

template <dim_t D, dim_t LEVEL>
typename GenzKeisterQR<D, LEVEL>::WeightVector
    GenzKeisterQR<D, LEVEL>::cached_weights;

}
