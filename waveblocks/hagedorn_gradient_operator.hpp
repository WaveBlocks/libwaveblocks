#ifndef WAVEBLOCKS_HAGEDORN_GRADIENT_OPERATOR_HPP
#define WAVEBLOCKS_HAGEDORN_GRADIENT_OPERATOR_HPP

#include "basic_types.hpp"

#include "hagedorn_parameter_set.hpp"
#include "shape_enumeration_default.hpp"
#include "shape_extended.hpp"
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
private:
    std::shared_ptr< ShapeEnumeration<D> > grad_enum_;
    
    /**
     * \param[in] wavepacket a hagedorn wavepacket
     * \return
     * \parblock
     * A tuple that contains for each dimension of the gradient its coefficient vector.
     * \endparblock
     */
    std::array< std::vector<complex_t>, D > apply_(const HagedornWavepacket<D>& wavepacket) const
    {
        const auto eps = wavepacket.eps();
        const auto & p = wavepacket.parameters().p;
        const auto & P = wavepacket.parameters().P;
        const auto base_enum = wavepacket.enumeration();
        const auto & base_coeffs = wavepacket.coefficients();
        
        Eigen::Matrix<complex_t,D,D> Pbar = P.conjugate();
        
        std::array< std::vector<complex_t>, D > grad_coeffs;
        for (dim_t d = 0; d < D; d++) {
            grad_coeffs[d] = std::vector<complex_t>(grad_enum_.size());
        }
        
        //iterate over each slice [i = index of current slice]
        for (std::size_t i = 0; i < grad_enum_->count_slices(); i++) {
            //loop over all multi-indices within current slice [j = position of multi-index within current slice]
            for (std::size_t j = 0; j < grad_enum_->slice(i).size(); j++) {
                std::array<int,D> curr_index = grad_enum_->slice(i)[j];
                
                //central node
                complex_t cc;
                if (base_enum->contains(curr_index)) {
                    std::size_t curr_ordinal = base_enum->slice(i).find(curr_index);
                    
                    assert (curr_ordinal < base_enum->slice(i).size());
                    
                    curr_ordinal += base_enum->slice(i).offset();
                    
                    cc = base_coeffs[curr_ordinal];
                }
                
                //backward neighbours
                Eigen::Matrix<complex_t,D,1> cb;
                for (dim_t d = 0; d < D; d++) {
                    if (curr_index[d] != 0) {
                        std::array<int,D> prev_index = curr_index; prev_index[d] -= 1;
                        if (base_enum->contains(prev_index)) {
                            std::size_t prev_ordinal = base_enum->slice(i-1).find(prev_index);
                            
                            assert (prev_ordinal < base_enum->slice(i-1).size());
                            
                            prev_ordinal += base_enum->slice(i-1).offset();
                            
                            cb[d] = base_coeffs[prev_ordinal] * std::sqrt(real_t(curr_index[d]));
                        }
                    }
                }
                
                //forward neighbours
                Eigen::Matrix<complex_t,D,1> cf;
                if (i+1 < base_enum->count_slices() && base_enum->contains(curr_index)) {
                    for (dim_t d = 0; d < D; d++) {
                        std::array<int,D> next_index = curr_index; next_index[d] += 1;
                        
                        if (base_enum->contains(next_index)) {
                            std::size_t next_ordinal = base_enum->slice(i+1).find(next_index);
                            
                            assert (next_ordinal < base_enum->slice(i+1).size());
                            
                            next_ordinal += base_enum->slice(i+1).offset();
                            
                            cf[d] = base_coeffs[next_ordinal] * std::sqrt(real_t(curr_index[d]+1));
                        }
                    }
                }
                
                Eigen::Matrix<complex_t,D,1> Pcc = p*cc;
                Eigen::Matrix<complex_t,D,1> Pcf = Pbar*cf;
                Eigen::Matrix<complex_t,D,1> Pcb = P*cb;
                
                Eigen::Matrix<complex_t,D,1> bgrad = (Pcf + Pcb)*eps/std::sqrt(real_t(2)) + Pcc;
                
                for (dim_t d = 0; d < D; d++) {
                    grad_coeffs[d][grad_enum_->slice(i).offset() + j] = bgrad(d,0);
                }
            }
        }
        
        return grad_coeffs;
    }
    
public:
    GradientOperator() = delete;
    
    GradientOperator(const std::shared_ptr< ShapeEnumeration<D> > &grad_enum)
        : grad_enum_(grad_enum)
    { }
    
    GradientOperator(const GradientOperator<D> &other) = default;
    
    GradientOperator &operator=(const GradientOperator<D> &other) = default;
    
    std::array< HagedornWavepacket<D>, D > operator()(const HagedornWavepacket<D>& wavepacket) const
    {
        auto grad_coeffs = apply_(wavepacket);
        
        std::array< HagedornWavepacket<D>, D > tuple;
        for (dim_t d = 0; d < D; d++) {
            tuple[d] = HagedornWavepacket<D>{wavepacket.eps(), wavepacket.parameters(), grad_enum_, std::move(grad_coeffs[d])};
        }
        
        return tuple;
    }
};

}

#endif