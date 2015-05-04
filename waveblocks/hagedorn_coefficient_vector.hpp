#ifndef WAVEBLOCKS_HAGEDORN_COEFFIENT_VECTOR_HPP
#define WAVEBLOCKS_HAGEDORN_COEFFIENT_VECTOR_HPP

#include <memory>

#include "basic_types.hpp"
#include "multi_index.hpp"
#include "sliced_shape_enumeration.hpp"

namespace waveblocks {

/**
 * A coefficent-vector is a 2-tuple that contains:
 *    [1] a pointer to the describing shape enumeration
 *    [2] a field containing the coefficient values
 *        the size of the field matches the size of the enumeration
 */
template<dim_t D>
class CoefficientVector
{
private:
    std::shared_ptr< const ShapeEnumeration<D> > enumeration_;
    std::vector<complex_t> values_;
    
public:
    CoefficientVector(const std::shared_ptr< const ShapeEnumeration<D> > &enumeration)
        : enumeration_(enumeration)
        , values_(enumeration->size())
    { }
    
    CoefficientVector(const std::shared_ptr< const ShapeEnumeration<D> > &enumeration,
                      const std::vector<complex_t> &values)
        : enumeration_(enumeration)
        , values_(values)
    { 
        assert (enumeration.size() == values.size());
    }
    
    /**
     * Copy constructor
     * 
     * Actions:
     *    [1] copy reference to source enumeration
     *    [2] copy all coefficent values from the source field
     */
    CoefficientVector(const CoefficientVector &other)
        : enumeration_(other.enumeration_)
        , values_(other.values_)
    { }
    
    /**
     * Conversion constructor
     * 
     * Actions:
     *    [A] if both vectors share same enumeration:
     *        [1] copy reference
     *        [2] copy values
     * 
     *    [B] enumeration is different
     *        @TODO transformation
     */
    CoefficientVector(const std::shared_ptr< const ShapeEnumeration<D> > &enumeration, 
                      const CoefficientVector<D> &other)
        : enumeration_(enumeration)
        , values_(enumeration.size())
    {
        if (enumeration_.get() == other.enumeration_.get()) {
            values_ = other.values_;
        } else {
            throw "not implemented yet";
        }
    }
    
    /**
     * Assignment Operator
     * 
     * equivalent to the copy constructor
     */
    CoefficientVector<D> &operator=(const CoefficientVector<D> &other)
    {
        enumeration_ = other.enumeration_;
        values_ = other.values_;
        
        return *this;
    }
    
    const std::shared_ptr< const ShapeEnumeration<D> > &enumeration() const
    {
        return enumeration_;
    }
    
    
    std::size_t size() const
    {
        return values_.size();
    }
    
    
    std::vector<complex_t>::iterator begin()
    {
        return values_.begin();
    }
    
    std::vector<complex_t>::const_iterator begin() const
    {
        return values_.begin();
    }
    
    
    std::vector<complex_t>::iterator end()
    {
        return values_.end();
    }
    
    std::vector<complex_t>::const_iterator end() const
    {
        return values_.end();
    }
    
    
    complex_t &operator[](std::size_t index)
    {
        return values_[index];
    }
    
    complex_t operator[](std::size_t index) const
    {
        return values_[index];
    }
};

}

#endif