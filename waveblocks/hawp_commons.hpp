#ifndef WAVEBLOCKS_HAWP_COMMONS_HPP
#define WAVEBLOCKS_HAWP_COMMONS_HPP

#include "hawp.hpp"

#include "shape_enum_union.hpp"
#include "shape_enum_extended.hpp"

#include <stdexcept>

namespace waveblocks
{

template<dim_t D, class MultiIndex>
using ShapeEnumSharedPtr = std::shared_ptr< ShapeEnum<D, MultiIndex> >;

// template<dim_t D>
// class AbstractScalarWavepacket
// {
// public:
//     virtual Eigen::Array<complex_t,1,Eigen::Dynamic> evaluate(ComplexGrid<D,Eigen::Dynamic> const& grid) const = 0;
//     virtual HaWpBasisVector<Eigen::Dynamic> evaluate_basis(ComplexGrid<D,Eigen::Dynamic> const& grid) const = 0;
// };

template<dim_t D, class MultiIndex>
class ShapeExtensionCache
{
public:
    /**
     * \brief shape reference to shape to check actuality of cache
     * \return 
     */
    ShapeEnumSharedPtr<D,MultiIndex> get_extended_shape(std::shared_ptr< ShapeEnum<D,MultiIndex> > shape) const
    {
        if (shape.get() != cached_extended_shape_source_)
            update_extended_shape(shape); // this function is not thread-safe
        
        return cached_extended_shape_;
    }
    
    /**
     * \brief Manually sets extended shape
     * 
     * \param shape source of new extended shape
     * \param extension new extended shape
     */
    void set_extended_shape(std::shared_ptr< ShapeEnum<D,MultiIndex> > shape, 
                            std::shared_ptr< ShapeEnum<D,MultiIndex> > extension)
    {
        cached_extended_shape_source_ = shape.get();
        cached_extended_shape_ = extension;
    }
    
    /**
     * \brief Recomputes extended shape if source shape changed.
     * 
     * \param shape new source shape
     */
    void update_extended_shape(std::shared_ptr< ShapeEnum<D,MultiIndex> > shape) const
    {
        if (shape.get() != cached_extended_shape_source_) {
            cached_extended_shape_source_ = shape.get();
            cached_extended_shape_ = std::make_shared< ShapeEnum<D,MultiIndex> >(shape_enum::extend(shape.get()));
        }
    }
    
private:
    mutable ShapeEnumSharedPtr<D,MultiIndex> cached_extended_shape_;
    mutable ShapeEnum<D,MultiIndex> * cached_extended_shape_source_; // source of the cached extended shape
};

template<dim_t D, class MultiIndex>
class AbstractScalarHaWpBasis
{
public:
    virtual double eps() const = 0;
    virtual HaWpParamSet<D> const& parameters() const = 0;
    virtual ShapeEnumSharedPtr<D, MultiIndex> shape() const = 0;
    
    template<int N>
    HaWpEvaluator<D,MultiIndex,N> create_evaluator(CMatrix<D,N> const& grid) const
    {
        return {eps(), &parameters(), shape().get(), grid};
    }
    
    template<int N>
    HaWpBasisVector<N> evaluate_basis(CMatrix<D,N> const& grid) const
    {
        return create_evaluator(grid).all();
    }
    
    template<int N>
    HaWpBasisVector<N> evaluate_basis(RMatrix<D,N> const& rgrid) const
    {
        CMatrix<D,N> cgrid = rgrid.template cast <complex_t>();
        return evaluate_basis(cgrid);
    }
    
    //     virtual HaWpBasisVector<Eigen::Dynamic> evaluate_basis(ComplexGrid<D,Eigen::Dynamic> const& grid) const
    //     {
    //         return create_evaluator<Eigen::Dynamic>(grid).all();
    //     }
    
    /**
     * \brief Computes the extended shape if necessary (caching!) and returns it.
     * 
     * Computing an extended shape is expensive. Therefore this function
     * manually managed cache.
     * 
     * Thus you need to call update_extended_shape() if you change
     * the shape of this wavepacket.
     * 
     * \e Thread-Safety: The stored pointer to the cached shape extension is not guarded by a mutex.
     * Therefore race conditions may occur when calling this function concurrently.
     * 
     * \return shared pointer to extended shape
     */
    ShapeEnumSharedPtr<D,MultiIndex> extended_shape() const
    {
        return shape_extension_cache_.get_extended_shape( this->shape() );
    }
    
    /**
     * \brief Manually updates the stored shape extension if necessary.
     * 
     * Be aware that this operation is EXPENSIVE.
     * Currently, generating a share extension scales O(D*log(D)*N) 
     * where D is the wavepacket dimensionality and N is the number of shape lattice points.
     */
    void update_extended_shape()
    {
        shape_extension_cache_.update_extended_shape( this->shape() );
    };
    
private:
    ShapeExtensionCache<D,MultiIndex> shape_extension_cache_;
};

/**
 * \brief Abstract basis class that represents a scalar (1-component) hagedorn wavepacket.
 * 
 * It provides read-only access to epsilon, shape, parameters and coefficients.
 * Furthermore it is able to evaluate itself on quadrature points.
 * 
 * \tparam D wavepacket dimensionality
 * \tparam MultiIndex
 */
template<dim_t D, class MultiIndex>
class AbstractScalarHaWp : public AbstractScalarHaWpBasis<D,MultiIndex>
{
public:
    virtual double eps() const = 0;
    virtual HaWpParamSet<D> const& parameters() const = 0;
    virtual ShapeEnumSharedPtr<D, MultiIndex> shape() const = 0;
    virtual std::vector<complex_t> const& coefficients() const = 0;
    
    template<int N>
    CMatrix<1,N> evaluate(CMatrix<D,N> const& grid) const
    {
        if (this->shape()->n_entries() != coefficients().size())
            throw std::runtime_error("shape.size() != coefficients.size()");
        
        return this->template create_evaluator<N>(grid).reduce( coefficients() );
    }
    
    template<int N>
    CMatrix<1,N> evaluate(RMatrix<D,N> const& rgrid) const
    {
        CMatrix<D,N> cgrid = rgrid.template cast <complex_t>();
        return evaluate(cgrid);
    }
    
//     virtual HaWpBasisVector<Eigen::Dynamic> evaluate(ComplexGrid<D,Eigen::Dynamic> const& grid) const
//     {
//         return create_evaluator<Eigen::Dynamic>().reduce(coefficients());
//     }
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
    
    HaWpParamSet<D> & parameters()
    {
        return parameters_;
    }
    
    HaWpParamSet<D> const& parameters() const override
    {
        return parameters_;
    }
    
    ShapeEnumSharedPtr<D, MultiIndex> & shape()
    {
        return shape_;
    }
    
    ShapeEnumSharedPtr<D, MultiIndex> shape() const override
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
    HaWpParamSet<D> parameters_;
    ShapeEnumSharedPtr<D, MultiIndex> shape_;
    std::vector<complex_t> coefficients_;
    
}; // class ScalarHaWp

/**
 * \brief Represents a Hagedorn wavepacket with C components.
 * All components share the same Hagedorn parameterset.
 * 
 * The number of components C is determined at runtime.
 * 
 * \tparam D wavepacket dimensionality
 * \tparam MultiIndex
 */
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
            : owner_(that.owner_)
            , shape_(that.shape_)
            , coefficients_(std::move(that.coefficients_))
        { }
        
        Component(Component const& that)
            : owner_(that.owner_)
            , shape_(that.shape_)
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
        
        HaWpParamSet<D> const& parameters() const override
        {
            return owner_->parameters();
        }
        
        
        
        ShapeEnumSharedPtr<D, MultiIndex> & shape()
        {
            return shape_;
        }
        
        ShapeEnumSharedPtr<D, MultiIndex> shape() const override
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
        
        ShapeEnumSharedPtr<D, MultiIndex> shape_;
        std::vector<complex_t> coefficients_;
    };
    
    HomogeneousHaWp(std::size_t n)
        : eps_(0.0)
        , parameters_()
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
    
    HaWpParamSet<D> & parameters()
    {
        return parameters_;
    }
    
    HaWpParamSet<D> const& parameters() const
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
    
    ShapeEnumSharedPtr<D,MultiIndex> compute_union_shape() const
    {
        std::vector< ShapeEnum<D,MultiIndex>* > list(n_components());
        for (std::size_t c = 0; c < n_components(); c++) {
            list[c] = components()[c].shape().get();
        }
        return std::make_shared< ShapeEnum<D,MultiIndex> >(shape_enum::strict_union(list));
    }
    
    template<int N>
    CMatrix<Eigen::Dynamic,N> evaluate(CMatrix<D,N> const& grid) const
    {
        CMatrix<Eigen::Dynamic,N> result(n_components(),grid.cols());
        
        for (std::size_t c = 0; c < n_components(); c++) {
            result.row(c) = component(c).evaluate(grid);
        }
        
        return result;
    }
    
    template<int N>
    CMatrix<Eigen::Dynamic,N> evaluate(RMatrix<D,N> const& rgrid) const
    {
        CMatrix<D,N> cgrid = rgrid.template cast<complex_t>();
        return evaluate(cgrid);
    }
    
private:
    double eps_;
    HaWpParamSet<D> parameters_;
    std::vector<Component> components_;
    
}; // class HomogeneousHaWp

/**
 * \brief Represents a Hagedorn wavepacket with C components.
 * All components have a different set of Hagedorn parameters, shapes and coefficients.
 * 
 * The number of components C is determined at runtime.
 * 
 * \tparam D wavepacket dimensionality
 * \tparam MultiIndex
 */
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
        
        
        
        HaWpParamSet<D> & parameters()
        {
            return parameters_;
        }
        
        HaWpParamSet<D> const& parameters() const override
        {
            return parameters_;
        }
        
        ShapeEnumSharedPtr<D, MultiIndex> & shape()
        {
            return shape_;
        }
        
        ShapeEnumSharedPtr<D, MultiIndex> shape() const override
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
        
        HaWpParamSet<D> parameters_;
        ShapeEnumSharedPtr<D, MultiIndex> shape_;
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
    
    template<int N>
    CMatrix<Eigen::Dynamic,N> evaluate(CMatrix<D,N> const& grid) const
    {
        CMatrix<Eigen::Dynamic,N> result(n_components(),grid.cols());
        
        for (std::size_t c = 0; c < n_components(); c++) {
            result.row(c) = component(c).evaluate(grid);
        }
        
        return result;
    }
    
    template<int N>
    CMatrix<Eigen::Dynamic,N> evaluate(RMatrix<D,N> const& rgrid) const
    {
        CMatrix<D,N> cgrid = rgrid.template cast<complex_t>();
        return evaluate(cgrid);
    }
    
    ShapeEnumSharedPtr<D,MultiIndex> compute_union_shape() const
    {
        std::vector< ShapeEnum<D,MultiIndex>* > list(n_components());
        for (std::size_t c = 0; c < n_components(); c++) {
            list[c] = components()[c].shape().get();
        }
        return std::make_shared< ShapeEnum<D,MultiIndex> >(shape_enum::strict_union(list));
    }
    
private:
    double eps_;
    std::vector<Component> components_;
}; // class InhomogeneousHaWp

} // namespace waveblocks

#endif