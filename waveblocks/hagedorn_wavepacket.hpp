#ifndef WAVEBLOCKS_HAGEDORN_WAVEPACKET
#define WAVEBLOCKS_HAGEDORN_WAVEPACKET

#include <Eigen/Core>
#include <Eigen/Dense>

#include <vector>
#include <memory>

#include "basic_types.hpp"
#include "hagedorn_parameter_set.hpp"
#include "sliced_shape_enumeration.hpp"
#include "hagedorn_coefficient_vector.hpp"

#include "kahan_sum.hpp"

namespace waveblocks {

template<dim_t D>
class HagedornWavepacket
{
private:
    std::shared_ptr< const HagedornParameterSet<D> > parameters_;
    std::shared_ptr< const CoefficientVector<D> > coefficients_;
    
    complex_t evaluateBasis(const Eigen::Matrix<complex_t,D,D> &Qinv,
                                const Eigen::Matrix<complex_t,D,D> &QhQinvt,
                                const HagedornParameterSet<D> &parameters,
                                dim_t axis, 
                                MultiIndex<D> k, 
                                complex_t curr_basis, 
                                const Eigen::Matrix<complex_t,D,1> &prev_bases,
                                const Eigen::Matrix<real_t,D,1> &x) const
    {
        //compute {sqrt(k[i])*phi[k-e[i]]}
        //  e[i]: unit vector aligned to i-th axis
        Eigen::Matrix<complex_t,D,1> prev_bases_scaled = prev_bases;
        for (dim_t d = 0; d < D; d++)
            prev_bases_scaled(d,0) *= std::sqrt( real_t(k[d]) );
        
        Eigen::Matrix<real_t,D,1> dx = x - parameters.q;
        
        complex_t temp = Qinv.row(axis)*dx;
        
        complex_t pr1 = std::sqrt(2.0)/parameters.eps * temp * curr_basis;
        complex_t pr2 = QhQinvt.row(axis)*prev_bases_scaled;
        
        return (pr1 - pr2) / std::sqrt( real_t(k[axis])+1.0);
    }
    
public:
    HagedornWavepacket(std::shared_ptr< const HagedornParameterSet<D> > parameters, 
                       std::shared_ptr< const CoefficientVector<D> > coefficients)
        : parameters_(parameters)
        , coefficients_(coefficients)
    { }
    
    HagedornWavepacket(const HagedornWavepacket<D> &other)
        : parameters_(other.parameters_)
        , coefficients_(other.coefficients_)
    { }
    
    HagedornWavepacket &operator=(const HagedornWavepacket<D> &other)
    {
        parameters_ = other.parameters_;
        coefficients_ = other.coefficients_;
        
        return *this;
    }
    
    std::shared_ptr< const HagedornParameterSet<D> > parameters() const
    {
        return parameters_;
    }
    
    std::shared_ptr< const CoefficientVector<D> > coefficients() const
    {
        return coefficients_;
    }
    
    std::shared_ptr< const ShapeEnumeration<D> > enumeration() const
    {
        return coefficients_->enumeration();
    }
    
    complex_t operator()(const Eigen::Matrix<real_t,D,1> &x) const
    {
        auto & parameters = *parameters_;
        auto & coefficients = *coefficients_;
        auto & enumeration = *coefficients_->enumeration();
        
        Eigen::Matrix<complex_t,D,D> Qinv = parameters.Q.inverse();
        Eigen::Matrix<complex_t,D,D> QhQinvt = parameters.Q.adjoint()*Qinv.transpose();
        
        std::vector<complex_t> curr_slice_values;
        std::vector<complex_t> prev_slice_values;
        std::vector<complex_t> next_slice_values;
        
        //Use Kahan's algorithm to accumulate bases with O(1) numerical error instead of O(Sqrt(N))
        KahanSum<complex_t> psi;
        
        complex_t phi0 = parameters_->evaluateGroundState(x);
        next_slice_values.push_back(phi0);
        
        psi += phi0*coefficients[0];
        
        //loop over all slices [i = index of next slice]
        for (std::size_t i = 1; i < enumeration.count_slices(); i++) {
            //exchange slices
            std::swap(prev_slice_values, curr_slice_values);
            std::swap(curr_slice_values, next_slice_values);
            next_slice_values = std::vector<complex_t>(enumeration.slice(i).size());
            
            //loop over all multi-indices within next slice [j = position of multi-index within next slice]
            for (std::size_t j = 0; j < enumeration.slice(i).size(); j++) {
                MultiIndex<D> next_index = enumeration.slice(i)[j];
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
                std::size_t curr_ordinal = enumeration.slice(i-1).find(curr_index);
                
                assert(curr_ordinal < enumeration.slice(i-1).size()); //assert that multi-index has been found within current slice
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
                        std::size_t prev_ordinal = enumeration.slice(i-2).find(prev_index);
                        assert (prev_ordinal < enumeration.slice(i-2).size()); //assert that multi-index has been found within previous slice
                        prev_bases[d] = prev_slice_values[prev_ordinal];
                    }
                }
                
                //compute basis value within next slice
                complex_t next_basis = evaluateBasis(Qinv, QhQinvt, parameters, axis, curr_index, curr_basis, prev_bases, x);
                next_slice_values[j] = next_basis;
                psi += next_basis*coefficients[enumeration.slice(i).offset() + j];
            }
        }
        
        return psi();
    }
};

}

#endif