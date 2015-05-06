#ifndef WAVEBLOCKS_GRADIENT_OPERATOR_HPP
#define WAVEBLOCKS_GRADIENT_OPERATOR_HPP

#include "basic_types.hpp"
#include "multi_index.hpp"

#include "hagedorn_parameter_set.hpp"
#include "sliced_shape_enumeration.hpp"
#include "shape_extension.hpp"
#include "hagedorn_wavepacket.hpp"
#include "kahan_sum.hpp"

#include <Eigen/Core>

#include <vector>
#include <array>
#include <valarray>

namespace waveblocks {

template<dim_t D>
class GradientOperator
{
public:
    std::shared_ptr< ShapeEnumeration<D> > base_enum_;
    std::shared_ptr< ShapeEnumeration<D> > grad_enum_;
    
public:
    GradientOperator(const std::shared_ptr< ShapeEnumeration<D> > &base_enum,
                     const std::shared_ptr< ShapeEnumeration<D> > &grad_enum_)
        : base_enum_(base_enum)
        , grad_enum_(grad_enum_)
    { }
    
    GradientOperator(const GradientOperator<D> &other)
        : base_enum_(other.base_enum_)
        , grad_enum_(other.grad_enum_)
    { }
    
    GradientOperator &operator=(const GradientOperator<D> &other)
    {
        base_enum_ = other.base_enum_;
        grad_enum_ = other.grad_enum_;
        
        return *this;
    }
    
    HagedornWavepacket<D,D> operator()(const HagedornWavepacket<D,1> &wavepacket) const
    {
        if (wavepacket.enumeration() != base_enum_)
            throw "incompatible shapes";
        
        const auto & q = wavepacket.parameters()->q;
        const auto & p = wavepacket.parameters()->p;
        const auto & Q = wavepacket.parameters()->Q;
        const auto & P = wavepacket.parameters()->P;
        
        Eigen::Matrix<complex_t,D,D> Pbar = P.conjugate();
        
        HagedornWavepacket<D,D> gradpacket(wavepacket.eps(), wavepacket.parameters(), grad_enum_);
        
        //iterate over each slice [i = index of current slice]
        for (std::size_t i = 0; i < grad_enum_->count_slices(); i++) {
            //loop over all multi-indices within current slice [j = position of multi-index within current slice]
            for (std::size_t j = 0; j < grad_enum_->slice(i).size(); j++) {
                MultiIndex<D> curr_index = grad_enum_->slice(i)[j];
                
                //central node
                complex_t cc;
                if (base_enum_->contains(curr_index)) {
                    std::size_t curr_ordinal = base_enum_->slice(i).find(curr_index);
                    
                    assert (curr_ordinal < base_enum_->slice(i).size());
                    
                    curr_ordinal += base_enum_->slice(i).offset();
                    
                    cc = wavepacket.coefficient(curr_ordinal)(0,0);
                }
                
                //backward neighbours
                Eigen::Matrix<complex_t,D,1> cb;
                for (dim_t d = 0; d < D; d++) {
                    if (curr_index[d] != 0) {
                        MultiIndex<D> prev_index = curr_index; prev_index[d] -= 1;
                        if (base_enum_->contains(prev_index)) {
                            std::size_t prev_ordinal = base_enum_->slice(i-1).find(prev_index);
                            
                            assert (prev_ordinal < base_enum_->slice(i-1).size());
                            
                            prev_ordinal += base_enum_->slice(i-1).offset();
                            
                            cb[d] = wavepacket.coefficient(prev_ordinal)(0,0) * std::sqrt(real_t(curr_index[d]));
                        }
                    }
                }
                
                //forward neighbours
                Eigen::Matrix<complex_t,D,1> cf;
                if (i+1 < base_enum_->count_slices() && base_enum_->contains(curr_index)) {
                    for (dim_t d = 0; d < D; d++) {
                        MultiIndex<D> next_index = curr_index; next_index[d] += 1;
                        
                        if (base_enum_->contains(next_index)) {
                            std::size_t next_ordinal = base_enum_->slice(i+1).find(next_index);
                            
                            assert (next_ordinal < base_enum_->slice(i+1).size());
                            
                            next_ordinal += base_enum_->slice(i+1).offset();
                            
                            cf[d] = wavepacket.coefficient(next_ordinal)(0,0) * std::sqrt(real_t(curr_index[d]+1));
                        }
                    }
                }
                
                Eigen::Matrix<complex_t,D,1> Pcc = p*cc;
                Eigen::Matrix<complex_t,D,1> Pcf = Pbar*cf;
                Eigen::Matrix<complex_t,D,1> Pcb = P*cb;
                
                Eigen::Matrix<complex_t,D,1> bgrad = (Pcf + Pcb)*wavepacket.eps()/std::sqrt(real_t(2)) + Pcc;
                
                for (dim_t d = 0; d < D; d++) {
                    const_cast< std::valarray<complex_t>& >(* (gradpacket.coefficients()[d]))[grad_enum_->slice(i).offset() + j] = bgrad(d,0);
                }
            }
        }
        
        return gradpacket;
    }
};

}

#endif