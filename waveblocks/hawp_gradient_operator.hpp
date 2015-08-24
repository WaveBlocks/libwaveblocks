#ifndef WAVEBLOCKS_GRADIENT_OPERATOR_HPP
#define WAVEBLOCKS_GRADIENT_OPERATOR_HPP

#include "hawp_commons.hpp"
#include "hawp_gradient.hpp"

namespace waveblocks
{

/**
 * \brief This class represents the gradient of a (scalar) wavepacket.
 * 
 * \tparam D dimensionality and number of components
 * \tparam MultiIndex
 */
template<dim_t D, class MultiIndex>
class HaWpGradient
{
public:
    class Component : public AbstractScalarHaWp<D,MultiIndex>
    {
        Component(HaWpGradient const * const owner)
            : owner_(owner)
        { }
        
        double eps() const override
        {
            return owner_->eps();
        }
        
        HaWpParamSet<D> const& parameters() const override
        {
            return owner_->parameters();
        }
        
        ShapeEnumSharedPtr<D,MultiIndex> const& shape() const override
        {
            return owner_->shape();
        }
        
        std::vector<complex_t> & coefficients()
        {
            return coefficients_;
        }
        
        std::vector<complex_t> const& coefficients() const
        {
            return coefficients_;
        }
        
    private:
        HaWpGradient const * const owner_;
        
        std::vector<complex_t> coefficients_;
    };
    
    HaWpGradient()
        : components_(Component(this))
    { }
    
    double & eps()
    {
        return eps_;
    }
    
    double const& eps() const
    {
        return eps_;
    }
    
    HaWpParamSet<D> & parameters()
    {
        return parameters_;
    }
    
    HaWpParamSet<D> const& parameters() const
    {
        return parameters_;
    }
    
    ShapeEnumSharedPtr<D,MultiIndex> & shape()
    {
        return shape_;
    }
    
    ShapeEnumSharedPtr<D,MultiIndex> const& shape() const
    {
        return shape_;
    }
    
    
    Component & component(std::size_t n)
    {
        return components_[n];
    }
    
    Component const& component(std::size_t n) const
    {
        return components_[n];
    }
    
    Component & operator[](std::size_t n)
    {
        return component(n);
    }
    
    Component const& operator[](std::size_t n) const
    {
        return component(n);
    }
    
    std::size_t n_components() const
    {
        return components_.size();
    }
    
private:
    double eps_;
    HaWpParamSet<D> parameters_;
    ShapeEnumSharedPtr<D,MultiIndex> shape_;
    std::array<Component,D> components_;
    
}; // class HaWpGradient

/**
 * \brief Instances of this class take a scalar Hagedorn wavepacket, apply the
 * gradient operator on it and return a D-component Hagedorn wavepacket.
 * 
 * All you have to do is:
 * \code{.cpp}
 * ScalarHaWp<D,MultiIndex> wavepacket;
 * 
 * // ... alter wavepacket ...
 * 
 * HaWpGradientOperator<D,MultiIndex> nabla;
 * HaWpGradient<D, MultiIndex> gradient = nabla(wavepacket);
 * \endcode
 * 
 * This class uses the HaWpGradientEvaluator to evaluate the coefficients of
 * the resulting wavepacket.
 * This class simplifies taking gradients since it assembles the resulting 
 * wavepacket (HaWpGradientEvaluator just returns the new coefficients).
 * 
 * You cannot apply the gradient to multi-component wavepackets (yet).
 * If you want the gradient of multi-component wavepackets you
 * will have to loop over all components and apply the gradient operator
 * on each component.
 * 
 * \tparam D dimensionality of the processed wavepackets
 * \tparam MultiIndex multi-index type of the processed wavepackets
 */
template<dim_t D, class MultiIndex>
class HaWpGradientOperator
{
public:
    HaWpGradient<D,MultiIndex> operator()(AbstractScalarHaWp<D,MultiIndex> const& wp) const
    {
        wp.compute_extended_shape();
        
        HaWpGradient<D,MultiIndex> gradwp;
        gradwp.eps() = wp.eps();
        gradwp.parameters() = wp.parameters();
        gradwp.shape() = wp.extended_shape();
        
        HaWpGradientEvaluator<D,MultiIndex> evaluator(wp.eps(), wp.parameters().get(), wp.shape().get(), gradwp.shape().get());
        std::array< std::vector<complex_t>, std::size_t(D) > coeffs_result = evaluator.apply(wp.coefficients());
        
        for (dim_t c = 0; c < D; c++) {
            gradwp[c].coefficients() = std::move(coeffs_result[c]);
        }
        
        return gradwp;
    }
    
}; // class HaWpGradientOperator

} // namespace waveblocks

#endif