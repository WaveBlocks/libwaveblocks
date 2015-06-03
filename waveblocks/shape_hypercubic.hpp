#ifndef WAVEBLOCKS_HYPERCUBIC_SHAPE
#define WAVEBLOCKS_HYPERCUBIC_SHAPE

#include <string>
#include <sstream>
#include <array>
#include <stdexcept>
#include <initializer_list>

#include "basic_types.hpp"

namespace waveblocks {

template<dim_t D>
class HyperCubicShape
{
private:
    std::array<int,D> limits_;
    
public:
    /**
     * \param[in] limits (exclusive)
     */
    HyperCubicShape(const std::array<int,D> &limits) 
        : limits_(limits)
    { }
    
    HyperCubicShape(int size)
    {
        for (std::size_t d = 0; d < D; d++)
            limits_[d] = size;
    }
    
    HyperCubicShape(std::initializer_list<int> list)
    {
        int deflt = 0;
        std::size_t i = 0;
        for (int e : list) {
            limits_[i++] = deflt = e;
        }
        //fill remaining elements with last value of initializer list
        while (i < D) {
            limits_[i++] = deflt;
        }
    }
    
    HyperCubicShape(const HyperCubicShape &that)
        : limits_(that.limits_)
    { }
    
    HyperCubicShape &operator=(const HyperCubicShape &that)
    {
        limits_ = that.limits_;
        return *this;
    }
    
    template<class MultiIndex>
    int limit(MultiIndex position, dim_t axis) const
    {
        { (void)(position); } //disable unused-parameter warning
        
        for (dim_t d = 0; d < D; d++) {
            if (d != axis && position[d] >= limits_[d])
                return -1;
        }
        return limits_[axis]-1;
    }
    
    /**
     * \return bbox inclusive
     */
    int bbox(dim_t axis) const
    {
        return limits_[axis]-1;
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