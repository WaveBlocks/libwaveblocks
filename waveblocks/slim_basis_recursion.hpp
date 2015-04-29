#ifndef WAVEBLOCKS_SLIM_BASIS_RECURSION
#define WAVEBLOCKS_SLIM_BASIS_RECURSION

#include <vector>
#include <map>

#include "basic_types.hpp"

#include "hagedorn_parameter_set.hpp"
#include "lexical_shape_enumerator.hpp"
#include "sliced_shape_enumeration.hpp"

namespace waveblocks {

template<dim_t D>
complex_t evaluateGroundState(const HagedornParameterSet<D> &parameters, 
                              const Eigen::Matrix<real_t,D,1> &x)
{
    const real_t pi = 3.14159265359;
    
    Eigen::Matrix<real_t,D,1> dx = x - parameters.q;
    
    complex_t pr1 = dx.transpose()*parameters.P*parameters.Q.inverse()*dx;
    real_t pr2 = parameters.p.transpose()*dx;
    
    complex_t exponent = complex_t(0.0, 1.0)/(parameters.eps*parameters.eps) * (0.5*pr1 + pr2);
    
    return 1.0/std::pow(pi*parameters.eps*parameters.eps, D/4.0) * std::exp(exponent);
}

template<dim_t D>
inline complex_t evaluateBasis(const HagedornParameterSet<D> &parameters,
                               dim_t axis, 
                               MultiIndex<D> k, 
                               complex_t curr_basis, 
                               const Eigen::Matrix<complex_t,D,1> &prev_bases,
                               const Eigen::Matrix<real_t,D,1> &x)
{
    Eigen::Matrix<complex_t,D,D> Qinv = parameters.Q.inverse();
    Eigen::Matrix<complex_t,D,D> QhQinvt = parameters.Q.adjoint()*Qinv.transpose();
    
    //compute {sqrt(k[i])*phi[k-e[i]]}
    //  e[i]: unit vector aligned to i-th axis
    Eigen::Matrix<complex_t,D,1> prev_bases_scaled = prev_bases;
    for (dim_t d = 0; d < D; d++)
        prev_bases_scaled(d,0) *= std::sqrt( real_t(k[d]) );
    
    Eigen::Matrix<real_t,D,1> dx = x - parameters.q;
    
    complex_t pr1 = std::sqrt(2.0)/parameters.eps * complex_t(Qinv.row(axis)*dx) * curr_basis;
    complex_t pr2 = QhQinvt.row(axis)*prev_bases_scaled;
    
    return (pr1 - pr2) / std::sqrt( real_t(k[axis])+1.0);
}

template<dim_t D, class S>
complex_t evaluateWavepacket(const std::vector<complex_t> &coefficients, 
                             const HagedornParameterSet<D> &parameters, 
                             const SlicedShapeEnumeration<D,S> &enumeration,
                             const Eigen::Matrix<real_t,D,1> &x)
{
    auto slices = enumeration.slices();
    
    std::vector<complex_t> curr_slice_values;
    std::vector<complex_t> prev_slice_values;
    std::vector<complex_t> next_slice_values;
    
    complex_t phi0 = evaluateGroundState(parameters, x);
    next_slice_values.push_back(phi0);
    complex_t result = phi0*coefficients[0];
    
    //loop over all slices [i = index of next slice]
    for (std::size_t i = 1; i < slices.count(); i++) {
        //exchange slices
        prev_slice_values = curr_slice_values;
        curr_slice_values = next_slice_values;
        next_slice_values.resize(slices[i].size());
        
        //loop over all multi-indices within next slice [j = position of multi-index within next slice]
        for (std::size_t j = 0; j < slices[i].size(); j++) {
            MultiIndex<D> next_index = slices[i][j];
            //find valid precursor: find first non-zero entry
            dim_t axis = D;
            for (dim_t d = 0; d < D; d++) {
                if (next_index[d] != 0) {
                    axis = d;
                    break;
                }
            }
            
            assert(axis != D); //assert that multi-index contains some non-zero entries
            
            //retrieve the basis value within current slice
            MultiIndex<D> curr_index = next_index; curr_index[axis] -= 1; //get backward neighbour
            std::size_t curr_ordinal = slices[i-1].find(curr_index);
            
            assert(curr_ordinal < slices[i-1].size()); //assert that multi-index has been found within current slice
            complex_t curr_basis = curr_slice_values[curr_ordinal];
            
            //retrieve the basis values within previous slice
            Eigen::Matrix<complex_t,D,1> prev_bases;
            for (dim_t d = 0; d < D; d++) {
                if (curr_index[d] == 0) {
                    //precursor is out of shape therefore set this precursor value to zero
                    prev_bases[d] = complex_t(0.0, 0.0);
                }
                else {
                    MultiIndex<D> prev_index = curr_index; prev_index[d] -= 1; //get backward neighbour
                    std::size_t prev_ordinal = slices[i-2].find(prev_index);
                    assert (prev_ordinal < slices[i-2].size()); //assert that multi-index has been found within previous slice
                    prev_bases[d] = prev_slice_values[prev_ordinal];
                }
            }
            
            //compute basis value within next slice
            complex_t next_basis = evaluateBasis(parameters, axis, curr_index, curr_basis, prev_bases, x);
            next_slice_values[j] = next_basis;
            result += next_basis*coefficients[slices[i].offset() + j];
        }
    }
    
    return result;
}

}

#endif