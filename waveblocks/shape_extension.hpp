#ifndef WAVEBLOCKS_SHAPE_EXTENSION_HPP
#define WAVEBLOCKS_SHAPE_EXTENSION_HPP

#include <vector>
#include <memory>
#include <algorithm>

#include "basic_types.hpp"
#include "multi_index.hpp"

#include "hypercubic_shape.hpp"

namespace waveblocks {

template<dim_t D, class S>
class ExtendedShape
{
private:
    S shape_;
    
public:
    ExtendedShape(S shape)
        : shape_(shape)
    { }
    
    int getSurface(dim_t axis, const MultiIndex<D> &index) const
    {
        int value = shape_.getSurface(axis, index);
        
        //extend only when there actually are some nodes
        if (value >= 0)
            value += 1;
        
        for (dim_t d = 0; d < D; d++) {
            if (d != axis && index[d] != 0) {
                //backward neighbour index
                MultiIndex<D> prev_index = index; prev_index[d] -= 1;
                
                value = std::max(value, shape_.getSurface(axis, prev_index));
            }
        }
        
        return value;
    }
};

template<dim_t D>
class ExtendedShape<D,HyperCubicShape<D>>
{
private:
    HyperCubicShape<D> expansion_;
    
public:
    ExtendedShape(HyperCubicShape<D> shape)
        : expansion_(shape)
    {
        MultiIndex<D> limits = shape.limits();
        for (dim_t d = 0; d < D; d++)
            limits[d] += 1;
        expansion_ = HyperCubicShape<D>(limits);
    }
    
    int getSurface(dim_t axis, const MultiIndex<D> &index) const
    {
        return expansion_.getSurface(axis,index);
    }
};

/**
 * Extended shapes are needed for the gradient computation of the wavepacket.
 * 
 * This class is not used anymore. The original idea was that the enumeration
 * of the extended shape consists of two parts: The first part is equal to the
 * enumeration of the original shape. This class forms the second part. It stores
 * the additional nodes that are part of the shape extension 'shell'. This saves
 * some memory as the first part is identical to the enumeration of the corresponding
 * basic shape.
 * Since the gradient computation involves the evaluation of a wavepacket, this
 * two-part idea implicates that two implementations of the evaluation are needed.
 * One implementation for a wavepacket using a basic shape and one for using 
 * an extended shape. This two-part idea has been abandoned for the sake of
 * simplicity.
 * 
 * This class has been fully tested and works.
 */
template<dim_t D, class S>
class ShapeExtensionEnumeration
{
private:
    S shape_;
    ExtendedShape<D,S> extension_;
    std::shared_ptr<std::vector<MultiIndex<D>>> table_;
    
public:
    ShapeExtensionEnumeration(S shape)
        : shape_(shape)
        , extension_(shape)
        , table_(std::make_shared<std::vector<MultiIndex<D>>>())
    {
        MultiIndex<D> index = {};
        
        while (true) {
            int imin = std::max(0, 1+shape_.getSurface(D-1,index));
            int imax = extension_.getSurface(D-1,index);
            for (dim_t i = imin; i <= imax; i++)
            {
                index[D-1] = i;
                table_->push_back(index);
            }
            index[D-1] = 0;
            
            if (D == 1) {
                return;
            }
            else {
                dim_t j = D-2;
                while (index[j] == extension_.getSurface(j, index)) {
                    index[j] = 0;
                    if (j == 0)
                        return;
                    else
                        j = j-1;
                }
                index[j] += 1;
            }
        }
    }
    
    ShapeExtensionEnumeration(const ShapeExtensionEnumeration &that)
        : shape_(that.shape_)
        , extension_(that.extension_)
        , table_(that.table_)
    { }
    
    ShapeExtensionEnumeration &operator=(const ShapeExtensionEnumeration &that)
    {
        shape_ = that.shape_;
        extension_ = that.extension_;
        table_ = that.table_;
        return *this;
    }
    
    std::size_t size() const
    {
        return table_->size();
    }
    
    std::size_t find(MultiIndex<D> index) const
    {
        std::less< MultiIndex<D> > comp;
        
        auto it = std::lower_bound(table_->begin(), table_->end(), index, comp);
        
        if (*it == index)
            return (it - table_->begin());
        else
            return table_->size(); //not found
    }
    
    MultiIndex<D> operator[](std::size_t ordinal) const
    {
        assert (ordinal < table_.size());
        
        return table_->at(ordinal);
    }
    
    typedef typename std::vector<MultiIndex<D>>::const_iterator Iterator;
    
    Iterator begin() const
    {
        return table_->begin();
    }
    
    Iterator end() const
    {
        return table_->end();
    }
};

}

#endif