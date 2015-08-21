#ifndef WAVEBLOCKS_HAGEDORN_WAVEPACKET
#define WAVEBLOCKS_HAGEDORN_WAVEPACKET

#include <Eigen/Core>
#include <Eigen/Dense>

#include <vector>
#include <array>
#include <memory>
#include <initializer_list>
#include <memory>

#include "basic_types.hpp"
#include "math_util.hpp"
#include "hawp_paramset.hpp"
#include "shape_enum.hpp"

#include "hawp_evaluator.hpp"
#include "hawp_gradient.hpp"

#include "kahan_sum.hpp"

namespace waveblocks {

/**
 * \brief represents a hagedorn wavepacket basis
 * 
 */
template<dim_t D, class MultiIndex>
struct HaWpBasis
{
public:
    double eps;
    const HaWpParamSet<D>* parameters;
    const ShapeEnum<D,MultiIndex>* enumeration;
    
    HaWpBasis(double eps, 
              const HaWpParamSet<D>* parameters,
              const ShapeEnum<D,MultiIndex>* enumeration)
        : eps(eps)
        , parameters(parameters)
        , enumeration(enumeration)
    { }
    
    template<int N>
    HaWpEvaluator<D, MultiIndex,N> at(const Eigen::Matrix<complex_t,D,N> &x) const
    {
        return {eps, parameters, enumeration, x};
    }
    
    template<int N>
    HaWpEvaluator<D, MultiIndex,N> at(const Eigen::Matrix<real_t,D,N> &x) const
    {
        Eigen::Matrix<complex_t,D,N> xtmp = x.template cast<complex_t>();
        return {eps, parameters, enumeration, xtmp};
    }
};

/**
 * \brief represents a hagedorn wavepacket
 * 
 */
template<dim_t D, class MultiIndex>
struct HaWp
{
public:
    HaWpBasis<D, MultiIndex> basis;
    const std::vector<complex_t>* coefficients;
    
    HaWp(double eps, 
         const HaWpParamSet<D>* parameters,
         const ShapeEnum<D,MultiIndex>* enumeration,
         const std::vector<complex_t>* coefficients)
        : basis(eps, parameters, enumeration)
        , coefficients(coefficients)
    { }
    
    HaWp(const HaWpBasis<D,MultiIndex>& basis,
         const std::vector<complex_t>* coefficients)
        : basis(basis)
        , coefficients(coefficients)
    { }
};

template<dim_t D, class MultiIndex>
using ShapeEnumPtr = std::shared_ptr< ShapeEnum<D, MultiIndex> >;

template<dim_t D>
class AbstractScalarWavepacket
{
public:
    virtual Eigen::Array<complex_t,1,Eigen::Dynamic> evaluate(ComplexGrid<D,Eigen::Dynamic> const& grid) const = 0;
    virtual HaWpBasisVector<Eigen::Dynamic> evaluate_basis(ComplexGrid<D,Eigen::Dynamic> const& grid) const = 0;
};

template<dim_t D, class MultiIndex>
class AbstractScalarHaWp : public AbstractScalarWavepacket<D>
{
public:
    virtual double eps() const = 0;
    virtual HaWpBasis<D, MultiIndex> const& paramaters() const = 0;
    virtual ShapeEnumPtr<D, MultiIndex> shape() const = 0;
    virtual std::vector<complex_t> const& coefficients() const = 0;
    
    template<int N>
    HaWpEvaluator<D,MultiIndex,N> create_evaluator(ComplexGrid<D,N> const& grid) const
    {
        return {eps(), &paramaters(), shape().get(), grid};
    }
    
//     template<int N>
//     Eigen::Array<complex_t,1,N> evaluate(ComplexGrid<D,N> const& grid) const
//     {
//         return create_evaluator<N>(grid).reduce(coefficients());
//     }
    
    virtual HaWpBasisVector<Eigen::Dynamic> evaluate(ComplexGrid<D,Eigen::Dynamic> const& grid) const
    {
        return create_evaluator<Eigen::Dynamic>().reduce(coefficients());
    }
    
//     template<int N>
//     HaWpBasisVector<N> evaluate_basis(ComplexGrid<D,N> const& grid) const
//     {
//         return create_evaluator(grid).all();
//     }
    
    virtual HaWpBasisVector<Eigen::Dynamic> evaluate_basis(ComplexGrid<D,Eigen::Dynamic> const& grid) const
    {
        return create_evaluator<Eigen::Dynamic>(grid).all();
    }
};

template<dim_t D, class MultiIndex>
class ScalarHaWp : public AbstractScalarHaWp<D, MultiIndex>
{
public:
    double & eps()
    {
        return eps_;
    }
    
    double eps() const override
    {
        return eps_;
    }
    
    HaWpBasis<D, MultiIndex> & paramaters()
    {
        return parameters_;
    }
    
    HaWpBasis<D, MultiIndex> const& paramaters() const override
    {
        return parameters_;
    }
    
    ShapeEnumPtr<D, MultiIndex> & shape()
    {
        return shape_;
    }
    
    ShapeEnumPtr<D, MultiIndex> shape() const override
    {
        return shape_;
    }
    
    std::vector<complex_t> & coefficients()
    {
        return coefficients_;
    }
    
    std::vector<complex_t> const& coefficients() const override
    {
        return coefficients_;
    }
    
private:
    double eps_;
    HaWpBasis<D, MultiIndex> parameters_;
    ShapeEnumPtr<D, MultiIndex> shape_;
    std::vector<complex_t> coefficients_;
};

template<dim_t D, class MultiIndex>
class HomogeneousHaWp
{
public:
    class Component : public AbstractScalarHaWp<D,MultiIndex>
    {
    public:
        Component(HomogeneousHaWp const* const owner)
            : owner_(owner)
            , coefficients_()
        { }
        
        Component(Component&& that)
            : shape_(that.shape_)
            , coefficients_(std::move(that.coefficients_))
        { }
        
        Component(Component const& that)
            : shape_(that.shape_)
            , coefficients_(that.coefficients_)
        { }
        
        Component & operator=(Component&& that)
        {
            shape_ = that.shape_;
            coefficients_ = std::move(that.coefficients_);
            return *this;
        }
        
        Component & operator=(Component const& that)
        {
            shape_ = that.shape_;
            coefficients_ = that.coefficients_;
            return *this;
        }
        
        double eps() const override
        {
            return owner_->eps();
        }
        
        HaWpBasis<D, MultiIndex> const& paramaters() const override
        {
            return owner_->paramaters();
        }
        
        
        
        ShapeEnumPtr<D, MultiIndex> & shape()
        {
            return shape_;
        }
        
        ShapeEnumPtr<D, MultiIndex> shape() const override
        {
            return shape_;
        }
        
        std::vector<complex_t> & coefficients()
        {
            return coefficients_;
        }
        
        std::vector<complex_t> const& coefficients() const override
        {
            return coefficients_;
        }
        
    private:
        HomogeneousHaWp const* const owner_;
        
        ShapeEnumPtr<D, MultiIndex> shape_;
        std::vector<complex_t> coefficients_;
    };
    
    HomogeneousHaWp(std::size_t n)
        : components_(n, Component(this))
    { }
    
    double & eps()
    {
        return eps_;
    }
    
    double eps() const
    {
        return eps_;
    }
    
    HaWpBasis<D, MultiIndex> & paramaters()
    {
        return parameters_;
    }
    
    HaWpBasis<D, MultiIndex> const& paramaters() const
    {
        return parameters_;
    }
    
    std::vector<Component> & components()
    {
        return components_;
    }
    
    std::vector<Component> const& components() const
    {
        return components_;
    }
    
    Component & operator[](std::size_t n)
    {
        return components_[n];
    }
    
    Component const& operator[](std::size_t n) const
    {
        return components_[n];
    }
    
    std::size_t n_components() const
    {
        return components_.size();
    }
    
private:
    double eps_;
    HaWpBasis<D, MultiIndex> parameters_;
    std::vector<Component> components_;
};

template<dim_t D, class MultiIndex>
class InhomogeneousHaWp
{
public:
    class Component : public AbstractScalarHaWp<D,MultiIndex>
    {
    public:
        Component(InhomogeneousHaWp const* const owner)
            : owner_(owner)
            , shape_()
            , coefficients_()
        { }
        
        Component(Component&& that)
            : parameters_(std::move(that.paramaters_))
            , shape_(that.shape_)
            , coefficients_(std::move(that.coefficients_))
        { }
        
        Component(Component const& that)
            : parameters_(that.paramaters_)
            , shape_(that.shape_)
            , coefficients_(that.coefficients_)
        { }
        
        Component & operator=(Component&& that)
        {
            parameters_ = std::move(that.paramaters_);
            shape_ = that.shape_;
            coefficients_ = std::move(that.coefficients_);
            return *this;
        }
        
        Component & operator=(Component const& that)
        {
            parameters_ = that.parameters_;
            shape_ = that.shape_;
            coefficients_ = that.coefficients_;
            return *this;
        }
        
        double eps() const override
        {
            return owner_->eps();
        }
        
        
        
        HaWpBasis<D, MultiIndex> & paramaters()
        {
            return parameters_;
        }
        
        HaWpBasis<D, MultiIndex> const& paramaters() const override
        {
            return parameters_;
        }
        
        ShapeEnumPtr<D, MultiIndex> & shape()
        {
            return shape_;
        }
        
        ShapeEnumPtr<D, MultiIndex> shape() const override
        {
            return shape_;
        }
        
        std::vector<complex_t> & coefficients()
        {
            return coefficients_;
        }
        
        std::vector<complex_t> const& coefficients() const override
        {
            return coefficients_;
        }
        
    private:
        InhomogeneousHaWp const* const owner_;
        
        HaWpBasis<D, MultiIndex> parameters_;
        ShapeEnumPtr<D, MultiIndex> shape_;
        std::vector<complex_t> coefficients_;
    };
    
    InhomogeneousHaWp(std::size_t n)
        : eps_()
        , components_(n, Component(this))
    { }
    
    double & eps()
    {
        return eps_;
    }
    
    double eps() const
    {
        return eps_;
    }
    
    std::vector<Component> & components()
    {
        return components_;
    }
    
    std::vector<Component> const& components() const
    {
        return components_;
    }
    
    Component & operator[](std::size_t n)
    {
        return components_[n];
    }
    
    Component const& operator[](std::size_t n) const
    {
        return components_[n];
    }
    
    std::size_t n_components() const
    {
        return components_.size();
    }
    
private:
    double eps_;
    std::vector<Component> components_;
};

namespace hawp {

template<dim_t D, class MultiIndex>
HaWpBasis<D,MultiIndex> basis(double eps, 
                              const HaWpParamSet<D>* parameters,
                              const ShapeEnum<D,MultiIndex>* enumeration)
{
    return {eps, parameters, enumeration};
}

template<dim_t D>
complex_t prefactor(const HaWpParamSet<D>& paramaters)
{
    return real_t(1)/paramaters.sqrt_detQ();;
}

template<dim_t D, class MultiIndex>
HaWpGradientOperator<D,MultiIndex> nabla(double eps, 
                                     const HaWpParamSet<D>* parameters,
                                     const ShapeEnum<D,MultiIndex>* base_enum,
                                     const ShapeEnum<D,MultiIndex>* grad_enum)
{
    return {eps, parameters, base_enum, grad_enum};
}

template<dim_t D, class MultiIndex>
HaWpGradientOperator<D,MultiIndex> nabla(const HaWpBasis<D,MultiIndex>& basis,
                                     const ShapeEnum<D,MultiIndex>* grad_enum)
{
    return {basis.eps, basis.parameters, basis.enumeration, grad_enum};
}

} // namespace hawp

}

#endif