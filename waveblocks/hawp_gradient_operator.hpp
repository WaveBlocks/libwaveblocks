#ifndef WAVEBLOCKS_GRADIENT_OPERATOR_HPP
#define WAVEBLOCKS_GRADIENT_OPERATOR_HPP

#include "hawp_commons.hpp"
#include "hawp_gradient_evaluator.hpp"

namespace waveblocks
{

/**
 * \brief This class represents the gradient of a (scalar) wavepacket.
 * 
 * \tparam D dimensionality and number of components
 * \tparam MultiIndex
 */
template<dim_t D, class MultiIndex>
class HaWpGradient : public AbstractScalarHaWpBasis<D,MultiIndex>
{
public:
    class Component : public AbstractScalarHaWp<D,MultiIndex>
    {
    public:
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
        
        ShapeEnumSharedPtr<D,MultiIndex> shape() const override
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
        : components_(D, Component(this))
    { }
    
    double & eps()
    {
        return eps_;
    }
    
    double eps() const override
    {
        return eps_;
    }
    
    HaWpParamSet<D> & parameters()
    {
        return parameters_;
    }
    
    HaWpParamSet<D> const& parameters() const override
    {
        return parameters_;
    }
    
    ShapeEnumSharedPtr<D,MultiIndex> & shape()
    {
        return shape_;
    }
    
    ShapeEnumSharedPtr<D,MultiIndex> shape() const override
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
    std::vector<Component> components_;
    
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
    /**
     * \e Thread-Safety: Computing the gradient involves creating a shape extension.
     * Since computing a shape extension is very expensive, shape extensions are cached.
     * Concurrently applying any gradient operator to the same wavepacket is unsafe (and is pointless anyway)
     * since cached shape extensions are stored inside the wavepacket objects without mutex guard.
     * Till now applying the same gradient operator to different wavepacket objects in parallel
     * is completely safe. But to ensure future compatibility, each thread should use its 
     * own gradient operator instance.
     * 
     * \param wp scalar wavepacket
     * \return D-dimensional wavepacket
     */
    HaWpGradient<D,MultiIndex> operator()(AbstractScalarHaWp<D,MultiIndex> const& wp) const
    {
        HaWpGradient<D,MultiIndex> gradwp;
        gradwp.eps() = wp.eps();
        gradwp.parameters() = wp.parameters();
        gradwp.shape() = wp.extended_shape();
        
        HaWpGradientEvaluator<D,MultiIndex> evaluator(wp.eps(), &wp.parameters(), wp.shape().get(), gradwp.shape().get());
        std::array< std::vector<complex_t>, std::size_t(D) > coeffs_result = evaluator.apply(wp.coefficients());
        
        for (dim_t c = 0; c < D; c++) {
            gradwp[c].coefficients() = std::move(coeffs_result[c]);
        }
        
        return gradwp;
    }
    
}; // class HaWpGradientOperator

} // namespace waveblocks

#endif