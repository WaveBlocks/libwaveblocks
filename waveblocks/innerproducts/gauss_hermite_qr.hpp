#pragma once

#include <tuple>
#include <vector>

#include "../types.hpp"

#include "tables_gausshermite.hpp"


namespace waveblocks {
    namespace innerproducts {
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
             * \brief Return the number of nodes for the given order.
             *
             * In the case of GaussHermiteQR, this is the same as the order.
             */
            static dim_t number_nodes()
            {
                return order;
            }

            /**
             * \brief Return the quadrature nodes.
             */
            static NodeMatrix nodes()
            {
                const auto& nodes = gauss_hermite_rules[ORDER-1].nodes;
                return NodeMatrix::Map(nodes.data(), nodes.size());
            }

            /**
             * \brief Return the quadrature weights.
             */
            static WeightVector weights()
            {
                const auto& weights = gauss_hermite_rules[ORDER-1].weights;
                return WeightVector::Map(weights.data(), weights.size());
            }

            /**
             * \brief Return the quadrature nodes and weights.
             */
            static std::tuple<NodeMatrix,WeightVector> nodes_and_weights()
            {
                return std::make_tuple(nodes(), weights());
            }

            /**
             * \brief Free the precalculated nodes and weights.
             *
             * Does nothing, only for consistency with TensorProductQR.
             */
            static void clear_cache() {}
        };
    }
}
