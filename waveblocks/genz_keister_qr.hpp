#pragma once

#include <iostream>
#include <tuple>
#include <vector>

#include <Eigen/Core>
#include <Eigen/Dense>

#include "basic_types.hpp"
#include "tables_genzkeister.hpp"
#include "util/combinatorics.hpp"

namespace waveblocks {


/**
 * \brief Structure providing weighted nodes for Genz-Keister quadrature.
 *
 * \tparam D dimensionality of the Genz-Keister rule
 * \tparam LEVEL the level of the Genz-Keister rule, must be between 1 and 30
 *   inclusive
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

        // Check LEVEL validity.
        if (LEVEL < 1 || LEVEL > 30)
        {
            std::cerr << "GenzKeisterQR: LEVEL must be between 1 and 30!\n";
        }

        const dim_t K = LEVEL - 1;
        const partitions_t<D> parts = partitions<D>(K);
        dim_t node_offset = 0;
        for (const auto& part : parts)
        {
            int sum = 0;
            for (int p : part) sum += p + genz_keister_zs(p);
            if (sum <= K)
            {
                // Calculate nodes and weights for this partition.
                const NodeMatrix part_nodes = nodes_for_partition(part);
                const dim_t n_nodes = part_nodes.cols();
                const real_t part_weight = weight_for_partition(part);

                // DEBUG
                std::cout << "\nPart:";
                for (int p : part) std::cout << " " << p;
                std::cout << "\nNodes:\n" << part_nodes << "\n";
                std::cout << "\nWeight:\n" << part_weight << "\n";

                // Append them to the result.
                cached_nodes.conservativeResize(Eigen::NoChange, node_offset + n_nodes);
                cached_weights.conservativeResize(Eigen::NoChange, node_offset + n_nodes);

                cached_nodes.rightCols(n_nodes) = part_nodes;
                cached_weights.tail(n_nodes).fill(part_weight);

                node_offset += n_nodes;
            }
        }

        // Transform weights.
        //cached_weights /= exp(-

        cached = true;
    }

    /**
     * \brief Return fully symmetric quadrature nodes for given partition.
     *
     * \param[in] part partition
     */
    static NodeMatrix nodes_for_partition(const partition_t<D>& part)
    {
        const dim_t xi = nz<D>(part);
        const permutations_t<D> perms = permutations<D>(part);
        NodeMatrix nodes = NodeMatrix::Zero(D, perms.size() * (1 << (D-xi)));

        dim_t node_offset = 0;
        for (const auto& perm : perms)
        {
            // Copy nodes with flipped signs from non-zero generators.
            int nnz_row = 0;
            for (int row = 0; row < D; ++row)
            {
                if (perm[row] == 0) continue;

                for (int col = 0; col < (1 << (D-xi)); ++col)
                {
                    // Generators corresponding to current permutation.
                    const int sign = -2 * ((col >> nnz_row) & 1) + 1;
                    nodes(row, node_offset + col) = sign *
                        genz_keister_generators(perm[row]);
                }

                // Count non-zero rows.
                ++ nnz_row;
            }

            node_offset += (1 << (D-xi));
        }

        return nodes;
    }

    /**
     * \brief Return quadrature weight for given partition.
     *
     * \param[in] part partition
     */
    static real_t weight_for_partition(const partition_t<D>& part)
    {
        // TODO
        return 0;
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
