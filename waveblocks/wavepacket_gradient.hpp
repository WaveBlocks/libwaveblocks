#ifndef WAVEBLOCKS_WAVEPACKET_GRADIENT
#define WAVEBLOCKS_WAVEPACKET_GRADIENT

#include "basic_types.hpp"
#include "multi_index.hpp"

#include "hagedorn_parameter_set.hpp"
#include "sliced_shape_enumeration.hpp"

#include <Eigen/Core>

#include <map>
#include <vector>

namespace waveblocks {

/**
 * Computes gradient by scatter-type stencil application.
 * See Chapter 3.8.1 for details.
 */
template<dim_t D, class S>
void evaluateWavepacketGradient(const std::vector<complex_t> &coefficients, 
                                const HagedornParameterSet<D> &parameters, 
                                const SlicedShapeEnumeration<D,S> &slices,
                                const Eigen::Matrix<real_t,D,1> &x,
                                dim_t axis,
                                std::vector<complex_t> &result,
                                std::map<MultiIndex<D>, complex_t> &extension)
{
    //iterate over each slice [i = index of current slice]
    for (std::size_t i = 0; i < slices.count(); i++) {
        //loop over all multi-indices within current slice [j = position of multi-index within current slice]
        for (std::size_t j = 0; j < slices[i].size(); j++) {
            //central node
            
        }
    }
}

}

#endif