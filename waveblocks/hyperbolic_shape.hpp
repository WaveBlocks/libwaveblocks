#ifndef HYPERBOLIC_SHAPE_HPP
#define HYPERBOLIC_SHAPE_HPP

#include <cmath>

#include "multi_index.hpp"

namespace waveblocks {

template<dim_t D>
class HyperbolicCutShape
{
private:
    double K_;
    
public:
    HyperbolicCutShape(double K) : K_(K) {}
    
    HyperbolicCutShape(const HyperbolicCutShape &that) : K_(that.K_) {}
    
    HyperbolicCutShape &operator=(const HyperbolicCutShape &that)
    {
        K_ = that.K_;
        return *this;
    }
    
    MultiIndex<D> getLimits() const
    {
        MultiIndex<D> bbox;
        for (dim_t i = 0; i < D; i++)
            bbox[i] = std::floor(K_-1);
        return bbox;
    }
    
    int getSurface(dim_t axis, const MultiIndex<D> &coordinate) const
    {
        //choose large enough integer type to prevent overflow bugs
        long long product = 1;
        for (dim_t i = 0; i < D; i++)
            if (i != axis)
                product *= (1 + coordinate[i]);
        
        //std::floor is necessary since C++ rounds negative numbers up (so -0.8 becomes 0)
        return (int)std::floor(K_/product - 1);
    }
};

template<dim_t D>
class LimitedHyperbolicCutShape
{
private:
    double K_;
    MultiIndex<D> limits_;
    
public:
    LimitedHyperbolicCutShape(double K, MultiIndex<D> limits) : K_(K), limits_(limits) {}
    
    LimitedHyperbolicCutShape(const LimitedHyperbolicCutShape &that) : K_(that.K_), limits_(that.limits_) {}
    
    LimitedHyperbolicCutShape &operator=(const LimitedHyperbolicCutShape &that)
    {
        K_ = that.K_;
        limits_ = that.limits_;
        return *this;
    }
    
    MultiIndex<D> getLimits() const
    {
        MultiIndex<D> bbox;
        for (dim_t i = 0; i < D; i++) {
            bbox[i] = std::min((int)limits_[i]-1,(int)std::floor(K_-1));
        }
        return bbox;
    }
    
    int getSurface(std::size_t axis, const MultiIndex<D> &coordinate) const
    {
        //choose large enough integer type to prevent overflow bugs
        long long product = 1;
        for (dim_t i = 0; i < D; i++) {
            if (i != axis) {
                product *= (1 + coordinate[i]);
            }
        }
        
        //std::floor is necessary since C++ rounds negative numbers up (so -0.5 becomes 0)
        return std::min((int)limits_[axis]-1, (int)std::floor(K_/product - 1));
    }
};

}

#endif