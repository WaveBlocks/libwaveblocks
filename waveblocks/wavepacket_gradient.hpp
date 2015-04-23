#ifndef WAVEBLOCKS_WAVEPACKET_GRADIENT
#define WAVEBLOCKS_WAVEPACKET_GRADIENT

#include "basic_types.hpp"
#include "multi_index.hpp"

#include "hagedorn_parameter_set.hpp"
#include "sliced_shape_enumeration.hpp"
#include "shape_extension.hpp"

#include <Eigen/Core>

#include <map>
#include <vector>

namespace waveblocks {

/**
 * Computes gradient by gather-type stencil application.
 * See master thesis chapter 3.8 for details.
 */
template<dim_t D, class S>
void evaluateWavepacketGradientGatherType(const std::vector<complex_t> &coefficients, 
                                const HagedornParameterSet<D> &parameters, 
                                const SlicedShapeEnumeration<D,S> &enumeration,
                                const ShapeExtensionEnumeration<D,S> &extension,
                                const Eigen::Matrix<real_t,D,1> &x,
                                dim_t axis,
                                std::vector<complex_t> &result)
{
    auto &eps = parameters.eps;
    auto &P = parameters.P;
    auto &p = parameters.p;
    
    Eigen::Matrix<complex_t,D,D> Pbar = P * P.adjoint() * P.inverse().transpose();
    
    std::size_t shape_size = enumeration.size();
    
    result.clear();
    result.resize(shape_size + extension.size());
    
    //iterate over each slice [i = index of current slice]
    auto slices = enumeration.slices();
    for (std::size_t i = 0; i < slices.count(); i++) {
        //loop over all multi-indices within current slice [j = position of multi-index within current slice]
        for (std::size_t j = 0; j < slices[i].size(); j++) {
            MultiIndex<D> curr_index = slices[i][j];
            std::size_t curr_ordinal = slices[i].offset() + j;
            
            //central node
            complex_t cc = coefficients[curr_ordinal]*p(axis,0);
            
            //backward neighbours
            Eigen::Matrix<complex_t,D,1> cb;
            for (dim_t d = 0; d < D; d++) {
                if (curr_index[d] != 0) {
                    MultiIndex<D> prev_index = curr_index; prev_index[d] -= 1;
                    std::size_t prev_ordinal = slices[i-1].find(prev_index);
                    assert (prev_ordinal < slices[i-1].size()); //assert that backward neighbour exist
                    prev_ordinal += slices[i-1].offset();
                    
                    cb[d] = std::sqrt(real_t(curr_index[d])) * coefficients[prev_ordinal];
                }
            }
            
            //forward neighbours
            Eigen::Matrix<complex_t,D,1> cf;
            for (dim_t d = 0; d < D; d++) {
                MultiIndex<D> next_index = curr_index; next_index[d] += 1;
                
                bool is_ext = true;
                std::size_t next_ordinal;
                if (i+1 != slices.count()) {
                    next_ordinal = slices[i+1].find(next_index);
                    if (next_ordinal < slices[i+1].size()) {
                        is_ext = false;
                        next_ordinal += slices[i+1].offset();
                    }
                }
                if (is_ext) {
                    next_ordinal = extension.find(next_index);
                    assert (next_ordinal < extension.size()); //assert that forward neighbour exists
                    next_ordinal += enumeration.size();
                }
                
                cf[d] = std::sqrt(real_t(curr_index[d]+1)) * coefficients[next_ordinal];
            }
            
            result[curr_ordinal] = eps/std::sqrt(real_t(2))*(Pbar.row(axis)*cf + P.row(axis)*cb) + cc*p(axis,0);
        }
    }
    
    std::size_t curr_ordinal = enumeration.size();
    for (auto curr_index : extension) {
        //backward neighbours
        Eigen::Matrix<complex_t,D,1> cb;
        for (dim_t d = 0; d < D; d++) {
            if (curr_index[d] != 0) {
                MultiIndex<D> prev_index = curr_index; prev_index[d] -= 1;
                std::size_t prev_ordinal = extension.find(prev_index);
                if (prev_ordinal < extension.size())
                    cb[d] = std::sqrt(real_t(prev_index[d] + 1))*coefficients[prev_ordinal];
            }
        }
        
        result[curr_ordinal] = eps/std::sqrt(real_t(2))*(P.row(axis)*cb);
        
        ++curr_ordinal;
    }
}

/**
 * Computes gradient by scatter-type stencil application.
 * See master thesis chapter 3.8 for details.
 */
template<dim_t D, class S>
void evaluateWavepacketGradient(const std::vector<complex_t> &coefficients, 
                                const HagedornParameterSet<D> &parameters, 
                                const SlicedShapeEnumeration<D,S> &enumeration,
                                const ShapeExtensionEnumeration<D,S> &extension,
                                const Eigen::Matrix<real_t,D,1> &x,
                                dim_t axis,
                                std::vector<complex_t> &result)
{
    auto &eps = parameters.eps;
    auto &P = parameters.P;
    auto &p = parameters.p;
    
    Eigen::Matrix<complex_t,D,D> Pbar = P * P.adjoint() * P.inverse().transpose();
    
    std::size_t shape_size = enumeration.size();
    
    result.clear();
    result.resize(shape_size + extension.size());
    
    //iterate over each slice [i = index of current slice]
    auto slices = enumeration.slices();
    for (std::size_t i = 0; i < slices.count(); i++) {
        //loop over all multi-indices within current slice [j = position of multi-index within current slice]
        for (std::size_t j = 0; j < slices[i].size(); j++) {
            MultiIndex<D> curr_index = slices[i][j];
            std::size_t curr_ordinal = slices[i].offset() + j;
            complex_t coefficient = coefficients[curr_ordinal];
            
            //central node
            result[curr_ordinal] += coefficient*p(axis,0);
            
            //backward neighbours
            for (dim_t d = 0; d < D; d++) {
                if (curr_index[d] != 0) {
                    MultiIndex<D> prev_index = curr_index; prev_index[d] -= 1;
                    std::size_t prev_ordinal = slices[i-1].find(prev_index);
                    assert (prev_ordinal < slices[i-1].size()); //assert that backward neighbour exist
                    prev_ordinal += slices[i-1].offset();
                    
                    result[prev_ordinal] += 
                            eps/std::sqrt(real_t(2)) *
                            std::sqrt(real_t(curr_index[d])) *
                            coefficient *
                            Pbar(axis,d);
                }
            }
            
            //forward neighbours
            for (dim_t d = 0; d < D; d++) {
                MultiIndex<D> next_index = curr_index; next_index[d] += 1;
                
                bool is_ext = true;
                std::size_t next_ordinal;
                if (i+1 != slices.count()) {
                    next_ordinal = slices[i+1].find(next_index);
                    if (next_ordinal < slices[i+1].size()) {
                        is_ext = false;
                        next_ordinal += slices[i+1].offset();
                    }
                }
                if (is_ext) {
                    next_ordinal = extension.find(next_index);
                    assert (next_ordinal < extension.size()); //assert that forward neighbour exists
                    next_ordinal += enumeration.size();
                }
                
                result[next_ordinal] += 
                        eps/std::sqrt(real_t(2)) *
                        std::sqrt(real_t(curr_index[d]+1)) *
                        coefficient * 
                        P(axis,d);
            }
        }
    }
}

}

#endif