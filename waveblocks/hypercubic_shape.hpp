#ifndef WAVEBLOCKS_HYPERCUBIC_SHAPE
#define WAVEBLOCKS_HYPERCUBIC_SHAPE

#include <string>
#include <sstream>

#include "multi_index.hpp"

namespace waveblocks {

template<dim_t D>
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
    
    int getSurface(dim_t axis, MultiIndex<D> position) const
    {
        { (void)(position); } //disable unused-parameter warning
        
        for (dim_t d = 0; d < D; d++) {
            if (d != axis && position[d] >= limits_[d])
                return -1;
        }
        return limits_[axis]-1;
    }
    
    MultiIndex<D> limits() const
    {
        return limits_;
    }
    
    std::string description() const
    {
        std::stringstream out;
        out << "HyperCubicShape{";
        out << "dimension: " << D << ", ";
        out << "limits (exclusive): " << limits_ << "}";
        return out.str();
    }
};

}

#endif