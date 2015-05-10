#ifndef WAVEBLOCKS_SHAPE_EXTENSION_HPP
#define WAVEBLOCKS_SHAPE_EXTENSION_HPP

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
    
    std::string description() const
    {
        std::stringstream out;
        out << "ExtendedShape<?>[" << shape_.description() << "]";
        return out.str();
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
    
    std::string description() const
    {
        return "ExtendedHyperCubicShape";
    }
};

}

#endif