#ifndef WAVEBLOCKS_SHAPE_EXTENSION_HPP
#define WAVEBLOCKS_SHAPE_EXTENSION_HPP

#include "basic_types.hpp"

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
    
    int bbox(dim_t axis) const
    {
        return shape_.bbox(axis)+1;
    }
    
    template<class MultiIndex>
    int limit(const MultiIndex &index, dim_t axis) const
    {
        int value = shape_.limit(index, axis);
        
        //extend only when there actually are some nodes
        if (value >= 0)
            value += 1;
        
        for (dim_t d = 0; d < D; d++) {
            if (d != axis && index[d] != 0) {
                //backward neighbour index
                MultiIndex prev_index = index; prev_index[d] -= 1;
                
                value = std::max(value, shape_.limit(prev_index, axis) );
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
        std::array<int,D> limits2;
        for (dim_t d = 0; d < D; d++)
            limits2[d] = shape.bbox(d)+2;
        expansion_ = HyperCubicShape<D>(limits2);
    }
    
    int bbox(dim_t axis) const
    {
        return expansion_.bbox(axis);
    }
    
    template<class MultiIndex>
    int limit(const MultiIndex &index, dim_t axis) const
    {
        return expansion_.template limit< MultiIndex >(index, axis);
    }
    
    std::string description() const
    {
        return "ExtendedHyperCubicShape";
    }
};

}

#endif