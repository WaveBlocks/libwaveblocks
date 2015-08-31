#ifndef WAVEBLOCKS_GRADIENT_OPERATOR_HPP
#define WAVEBLOCKS_GRADIENT_OPERATOR_HPP

#include "hawp_commons.hpp"
#include "hawp_gradient_evaluator.hpp"

namespace waveblocks
{

/**
 * \brief This class represents the gradient \f$ \nabla \Phi \f$ of a (scalar) Hagedorn wavepacket \f$ \Phi \f$.
 * 
 * The gradient of a \f$ D \f$-dimensional wavepacket has \f$ D \f$ components.
 * 
 * Creation
 * --------
 * Apply HaWpGradientOperator to an AbstractScalarHaWp to create an instance:
 * \code
 * HaWpGradientOperator<D,MultiIndex> nabla;
 * HaWpGradient<D,MultiIndex> gradwp = nabla(wavepacket);
 * \endcode
 * 
 * \tparam D dimensionality and number of components
 * \tparam MultiIndex
 */
template<dim_t D, class MultiIndex>
class HaWpGradient : public AbstractScalarHaWpBasis<D,MultiIndex>
{
public:
    /**
     * \brief This class is component of a Hagedorn wavepacket gradient.
     * 
     * Such a component is a full-fledged scalar Hagedorn wavepacket
     * that shares the scaling parameter \f$ \varepsilon \f$,
     * parameter set \f$ \Pi \f$ and basis shape extension \f$ \mathfrak{K}_{ext} \f$
     * with other components.
     */
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
        
        /**
         * \brief Grants writeable access to the coefficients \f$ c \f$
         * of the wavepacket.
         */
        std::vector<complex_t> & coefficients()
        {
            return coefficients_;
        }
        
        /**
         * \brief Grants read-only access to the coefficients \f$ c \f$
         * of the wavepacket.
         */
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
    
    /**
     * \brief Grants writeable access to the semi-classical scaling parameter 
     * \f$ \varepsilon \f$ of the wavepacket.
     */
    double & eps()
    {
        return eps_;
    }
    
    double eps() const override
    {
        return eps_;
    }
    
    /**
     * \brief Grants writeable access to the Hagedorn parameter set 
     * \f$ \Pi \f$ of the wavepacket.
     */
    HaWpParamSet<D> & parameters()
    {
        return parameters_;
    }
    
    HaWpParamSet<D> const& parameters() const override
    {
        return parameters_;
    }
    
    /**
     * \brief Grants access to the basis shape 
     * \f$ \mathfrak{K} \f$ of the wavepacket.
     * 
     * \return 
     * Reference to the shape enumeration pointer. 
     * You can assign a new pointer to it!
     */
    ShapeEnumSharedPtr<D,MultiIndex> & shape()
    {
        return shape_;
    }
    
    ShapeEnumSharedPtr<D,MultiIndex> shape() const override
    {
        return shape_;
    }
    
    /**
     * \brief Grants writeable access to the \f$ n \f$-th component
     * \f$ \Phi_n \f$.
     * 
     * \param n The index \f$ n \f$ of the requested component.
     * \return Reference to the requested component.
     */
    Component & component(std::size_t n)
    {
        return components_[n];
    }
    
    /**
     * \brief Grants read-only access to the \f$ n \f$-th component
     * \f$ \Phi_n \f$.
     * 
     * \param n The index \f$ n \f$ of the requested component.
     * \return Reference to the requested component.
     */
    Component const& component(std::size_t n) const
    {
        return components_[n];
    }
    
    /**
     * \brief Grants writeable access to the \f$ n \f$-th component
     * \f$ \Phi_n \f$.
     * 
     * \param n The index \f$ n \f$ of the requested component.
     * \return Reference to the requested component.
     */
    Component & operator[](std::size_t n)
    {
        return component(n);
    }
    
    /**
     * \brief Grants read-only access to the \f$ n \f$-th component
     * \f$ \Phi_n \f$.
     * 
     * \param n The index \f$ n \f$ of the requested component.
     * \return Reference to the requested component.
     */
    Component const& operator[](std::size_t n) const
    {
        return component(n);
    }
    
    /**
     * \brief Returns the number of components.
     */
    std::size_t n_components() const
    {
        return components_.size();
    }
    
    /**
     * \brief Evaluate the value of all components at once.
     * 
     * Evaluates \f$ \Psi(x) = \{\Phi_i(x)\} \f$, 
     * where \f$ x \f$ is is a complex quadrature point.
     * 
     * \param grid 
     * Complex quadrature points.
     * Complex matrix of shape (dimensionality, number of quadrature points)
     * \return
     * Complex matrix of shape (number of components, number of quadrature points)
     * 
     * \tparam N Number of quadrature points or Eigen::Dynamic if not known at compile-time.
     */
    template<int N>
    CMatrix<Eigen::Dynamic,N> evaluate(CMatrix<D,N> const& grid) const
    {
        CMatrix<Eigen::Dynamic,N> result(n_components(),grid.cols());
        
        for (std::size_t c = 0; c < n_components(); c++) {
            result.row(c) = component(c).evaluate(grid);
        }
        
        return result;
    }
    
    /**
     * \brief Evaluates the value of all components at once.
     * 
     * Evaluates \f$ \Psi(x) = \{\Phi_i(x)\} \f$, 
     * where \f$ x \f$ is is a real quadrature point.
     * 
     * \param rgrid
     * Real quadrature points.
     * Real matrix of shape (dimensionality, number of quadrature points)
     * \return
     * Complex matrix of shape (number of components, number of quadrature points)
     * 
     * \tparam N Number of quadrature points or Eigen::Dynamic if not known at compile-time.
     */
    template<int N>
    CMatrix<Eigen::Dynamic,N> evaluate(RMatrix<D,N> const& rgrid) const
    {
        CMatrix<D,N> cgrid = rgrid.template cast<complex_t>();
        return evaluate(cgrid);
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
     * \brief Applies this gradient operator to a \e scalar Hagedorn wavepacket.
     * 
     * \e Vectorial \e wavepackets: 
     * You cannot apply this function directly to vectorial wavepackets \f$ \Psi \f$. You have to
     * apply the gradient to each component \f$ \Phi_n \f$ (which is scalar) of the vectorial wavepacket.
     * \f$ \nabla \Psi = \left( \nabla \Phi_1, \dots, \nabla \Phi_N \right)^T \f$
     * 
     * \e Thread-Safety: Computing the gradient involves creating a shape extension.
     * Since computing a shape extension is very expensive, shape extensions are cached.
     * Concurrently applying any gradient operator to the same wavepacket is unsafe (and is pointless anyway)
     * since cached shape extensions are stored inside the wavepacket objects without mutex guard.
     * Till now applying the same gradient operator to different wavepacket objects in parallel
     * is completely safe. But to ensure future compatibility, each thread should use its 
     * own gradient operator instance.
     * 
     * \param wp scalar Hagedorn wavepacket
     * \return wavepacket gradient
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