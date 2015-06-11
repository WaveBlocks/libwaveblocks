#ifndef WAVEBLOCKS_HAGEDORN_GRADIENT_OPERATOR_HPP
#define WAVEBLOCKS_HAGEDORN_GRADIENT_OPERATOR_HPP

#include <vector>
#include <array>

#include <Eigen/Core>

#include "hawp_paramset.hpp"
#include "shape_enum.hpp"

namespace waveblocks {

template<dim_t D, class MultiIndex>
class GradientOperator
{
private:
    double eps_;
    
    const HaWpParamSet<D>* parameters_;
    
    /**
     * \brief enumeration of basic shape
     */
    const ShapeEnum<D,MultiIndex>* base_enum_;
    
    /**
     * \brief enumeration of extended shape
     */
    const ShapeEnum<D,MultiIndex>* grad_enum_;
    
public:
    GradientOperator() = default;
    GradientOperator(const GradientOperator& that) = default;
    GradientOperator(GradientOperator&& that) = default;
    
    GradientOperator &operator=(const GradientOperator& that) = default;
    GradientOperator &operator=(GradientOperator&& that) = default;
    
    GradientOperator(double eps,
                     const HaWpParamSet<D>* parameters,
                     const ShapeEnum<D,MultiIndex>* base_enum,
                     const ShapeEnum<D,MultiIndex>* grad_enum)
        : eps_(eps)
        , parameters_(parameters)
        , base_enum_(base_enum)
        , grad_enum_(grad_enum)
    { }
    
    /**
     * \param[in] wavepacket a hagedorn wavepacket
     * \return
     * \parblock
     * A tuple that contains for each dimension of the gradient its coefficient vector.
     * \endparblock
     */
    std::array< std::vector<complex_t>, std::size_t(D) > apply(const std::vector<complex_t>& base_coeffs) const
    {
        const auto & p = parameters_->p;
        const auto & P = parameters_->P;
        
        Eigen::Matrix<complex_t,D,D> Pbar = P.conjugate();
        
        std::array< std::vector<complex_t>, std::size_t(D) > grad_coeffs;
        for (dim_t d = 0; d < D; d++) {
            grad_coeffs[d] = std::vector<complex_t>(grad_enum_->n_entries());
        }
        
        //iterate over each slice [i = index of current slice]
        for (int i = 0; i < grad_enum_->n_slices(); i++) {
            //loop over all multi-indices within current slice [j = position of multi-index within current slice]
            for (std::size_t j = 0; j < grad_enum_->slice(i).size(); j++) {
                MultiIndex curr_index = grad_enum_->slice(i)[j];
                
                //central node
                complex_t cc;
                std::size_t curr_ordinal;
                bool central_node_exists;
                
                if ( (central_node_exists = base_enum_->slice(i).try_find(curr_index, curr_ordinal)) ) {
                    
                    curr_ordinal += base_enum_->slice(i).offset();
                    
                    cc = base_coeffs[curr_ordinal];
                }
                
                //backward neighbours
                Eigen::Matrix<complex_t,D,1> cb;
                for (dim_t d = 0; d < D; d++) {
                    if (curr_index[d] != 0) {
                        MultiIndex prev_index = curr_index; prev_index[d] -= 1;
                        
                        std::size_t prev_ordinal;
                        
                        if (base_enum_->slice(i-1).try_find(prev_index, prev_ordinal)) {
                            
                            prev_ordinal += base_enum_->slice(i-1).offset();
                            
                            cb[d] = base_coeffs[prev_ordinal] * std::sqrt(real_t(curr_index[d]));
                        }
                    }
                }
                
                //forward neighbours
                Eigen::Matrix<complex_t,D,1> cf;
                if (central_node_exists && i+1 < base_enum_->n_slices()) {
                    for (dim_t d = 0; d < D; d++) {
                        MultiIndex next_index = curr_index; next_index[d] += 1;
                        
                        std::size_t next_ordinal;
                        
                        if (base_enum_->slice(i+1).try_find(next_index, next_ordinal)) {
                            
                            next_ordinal += base_enum_->slice(i+1).offset();
                            
                            cf[d] = base_coeffs[next_ordinal] * std::sqrt(real_t(curr_index[d]+1));
                        }
                    }
                }
                
                Eigen::Matrix<complex_t,D,1> Pcc = p*cc;
                Eigen::Matrix<complex_t,D,1> Pcf = Pbar*cf;
                Eigen::Matrix<complex_t,D,1> Pcb = P*cb;
                
                Eigen::Matrix<complex_t,D,1> bgrad = (Pcf + Pcb)*eps_/std::sqrt(real_t(2)) + Pcc;
                
                for (dim_t d = 0; d < D; d++) {
                    grad_coeffs[d][grad_enum_->slice(i).offset() + j] = bgrad(d,0);
                }
            }
        }
        
        return grad_coeffs;
    }
};

} // namespace waveblocks

#endif