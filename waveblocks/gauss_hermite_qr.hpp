#pragma once

#include <tuple>
#include <vector>

#include "basic_types.hpp"
#include "tables_gausshermite.hpp"




namespace waveblocks {

    /**
     * \brief Structure providing weighted nodes for Gauss Hermite quadrature.
     *
     * \tparam ORDER requested order of the quadrature rule
     */
    template <dim_t ORDER>
    struct GaussHermiteQR
    {
        static const dim_t D = 1;
        static const dim_t order = ORDER;

        using NodeMatrix = Eigen::Matrix<real_t,1,Eigen::Dynamic>;
        using WeightVector = Eigen::Matrix<real_t,1,Eigen::Dynamic>;

        /**
         * \brief Returns the number of nodes for the given order.
         *
         * In the case of GaussHermiteQR, this is the same as the order.
         */
        static dim_t number_nodes()
        {
            return order;
        }

        /**
         * \brief Returns the quadrature nodes.
         */
        static const std::vector<real_t>& nodes()
        {
            return gauss_hermite_rules[ORDER-1].nodes;
        }


        /**
         * \brief Returns the quadrature weights.
         */
        static const std::vector<real_t>& weights()
        {
            return gauss_hermite_rules[ORDER-1].weights;
        }

        /**
         * \brief Returns the quadrature nodes and weights.
         */
        static std::tuple<NodeMatrix,WeightVector> nodes_and_weights()
        {
            return std::make_tuple(
                                   NodeMatrix::Map(nodes().data(), nodes().size()),
                                   WeightVector::Map(weights().data(), weights().size()));
        }
    };

}
