#ifndef WAVEBLOCKS_ABSTRACT_SHAPE_HPP
#define WAVEBLOCKS_ABSTRACT_SHAPE_HPP

#include <iostream>

namespace waveblocks
{
    template<dim_t D>
    class AbstractShape
    {
    public:
        virtual ~AbstractShape() { }
        
        virtual int bbox(dim_t axis) const = 0;
        virtual int limit(int const* base_node, dim_t axis) const = 0;
        
        virtual void print(std::ostream & out) const = 0;
    };
    
    template<dim_t D>
    std::ostream & operator<<(std::ostream & out, AbstractShape<D> const& shape)
    {
        shape.print(out);
        return out;
    }
}



#endif