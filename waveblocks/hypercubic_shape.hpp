#ifndef WAVEBLOCKS_HYPERCUBIC_SHAPE
#define WAVEBLOCKS_HYPERCUBIC_SHAPE

#include "multi_index.hpp"

namespace waveblocks {

template<std::size_t D>
class HyperCubicShape
{
private:
    MultiIndex<D> limits_;
    
public:
    HyperCubicShape(MultiIndex<D> limits) 
        : limits_(limits)
    { }
    
    HyperCubicShape(const HyperCubicShape &that)
        : limits_(that.limits_)
    { }
    
    HyperCubicShape &operator=(const HyperCubicShape &that)
    {
        limits_ = that.limits_;
        return *this;
    }
    
    int getSurface(std::size_t axis, MultiIndex<D> position) const
    {
        { (void)(position); } //disable unused-parameter warning
        return limits_[axis]-1;
    }
};

}

#endif