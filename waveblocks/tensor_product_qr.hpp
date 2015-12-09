#pragma once

#include <iostream>
#include <tuple>
#include <vector>

#include <Eigen/Core>
#include <Eigen/Dense>

#include "basic_types.hpp"
#include "tables_gausshermite.hpp"


namespace waveblocks {
    namespace innerproducts {

        /**
         * \brief Structure providing weighted nodes for Tensor Product quadrature.
         *
         * \tparam RULES list of other quadrature rules to use as components of the
         *   tensor product.
         */
        template <class... RULES>
        struct TensorProductQR
        {
            static const dim_t D = sizeof...(RULES);

            using NodeMatrix = Eigen::Matrix<real_t,D,Eigen::Dynamic>;
            using WeightVector = Eigen::Matrix<real_t,1,Eigen::Dynamic>;

            /**
             * \brief Return the number of nodes for the given order.
             *
             * In the case of TensorProductQR, this is the product of the numbers of
             * nodes of the constituent quadrature rules (RULES).
             */
            static dim_t number_nodes()
            {
                dim_t result = 1;
                for(auto n_nodes : { RULES::number_nodes() ... })
                    {
                        result *= n_nodes;
                    }
                return result;
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
                const dim_t dim = D;
                const dim_t n_nodes = number_nodes();
                cached_nodes.resize(dim, n_nodes);
                cached_weights.resize(1, n_nodes);
                cached_weights.setOnes();

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
                                    cached_nodes(d, n) = std::get<0>(nw)(base_idx);
                                    cached_weights(n) *= std::get<1>(nw)(base_idx);
                                }

                            ++d;
                        }
                }
                cached = true;
            }
        };

        // Initialize static members.
        template <class... RULES>
        bool TensorProductQR<RULES...>::cached = false;

        template <class... RULES>
        typename TensorProductQR<RULES...>::NodeMatrix
        TensorProductQR<RULES...>::cached_nodes;

        template <class... RULES>
        typename TensorProductQR<RULES...>::WeightVector
        TensorProductQR<RULES...>::cached_weights;
    }
}
