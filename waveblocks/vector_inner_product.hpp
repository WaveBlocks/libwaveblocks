#pragma once

#include <cmath>
#include <functional>
#include <iostream>
#include <tuple>
#include <vector>

#include <Eigen/Core>
#include <Eigen/Dense>

#include "basic_types.hpp"
#include "hawp_commons.hpp"

namespace waveblocks {

/**
 * \brief Class providing inner product calculation of multi-component
 *   wavepackets.
 *
 * \tparam D dimensionality of processed wavepackets
 * \tparam MultiIndex multi-index type of processed wavepackets
 * \tparam QR quadrature rule to use, with N nodes
 */
template<dim_t D, class MultiIndex, class QR>
class VectorInnerProduct
{
public:
};

}
