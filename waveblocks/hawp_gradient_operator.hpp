#ifndef WAVEBLOCKS_GRADIENT_OPERATOR_HPP
#define WAVEBLOCKS_GRADIENT_OPERATOR_HPP

#include "hawp_commons.hpp"
#include "hawp_gradient.hpp"

namespace waveblocks
{

template<dim_t D, class MultiIndex>
class HaWpGradient
{
public:
    class Component : public AbstractScalarHaWp<D,MultiIndex>
    {
        Component(HaWpGradient const * const owner)
            : owner_(owner)
        { }
        
        
    private:
        HaWpGradient const * const owner_;
        
        
    };
    
    HaWpGradient(double eps, HaWpParamSet<D> const& params, 
                 ShapeEnumSharedPtr<D,MultiIndex> shape,
                 std::array< std::vector<complex_t>, std::size_t(D) > const& coefficients)
    {
        
    }
    
private:
    
}; // class HaWpGradient

template<dim_t D, class MultiIndex>
class HaWpGradientOperator
{
public:
    HaWpGradient<D,MultiIndex> operator()(AbstractScalarHaWp<D,MultiIndex> const& wp) const
    {
        wp.compute_extended_shape();
        
        ShapeEnumSharedPtr<D,MultiIndex> extended_shape = wp.extended_shape();
        
        HaWpGradientEvaluator<D,MultiIndex> evaluator(wp.eps(), wp.parameters().get(), wp.shape().get(), extended_shape.get());
        std::array< std::vector<complex_t>, std::size_t(D) > coeffs_result = evaluator.apply(wp.coefficients());
        
        return {wp.eps(), wp.parameters(), extended_shape, coeffs_result};
    }
    
}; // class HaWpGradientOperator

} // namespace waveblocks

#endif