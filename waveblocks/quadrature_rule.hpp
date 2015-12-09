#pragma once

#include <vector>

#include "basic_types.hpp"

namespace waveblocks {
    namespace innerproducts {

        /**
         * \brief Node and weight values for a 1D quadrature rule.
         *
         * Should only be used internally.
         */
        struct QuadratureRule
        {
            const size_t order;
            const std::vector<real_t> nodes;
            const std::vector<real_t> weights;

            QuadratureRule(const size_t order,
                           const std::vector<real_t>& nodes,
                           const std::vector<real_t>& weights)
                : order(order), nodes(nodes), weights(weights)
            {}
        };
    }
}
