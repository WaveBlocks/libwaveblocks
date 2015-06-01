#ifndef HYPERBOLIC_SHAPE_HPP
#define HYPERBOLIC_SHAPE_HPP

#include <cmath>
#include <string>
#include <sstream>

#include "basic_types.hpp"

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
    
    int bbox(dim_t axis) const
    {
        (void)(axis); //suppress -Wunused-parameter
        return std::floor(K_-1.0);
    }
    
    template<class MultiIndex>
    int limit(const MultiIndex &coordinate, dim_t axis) const
    {
        //choose large enough integer type to prevent overflow
        long long product = 1;
        for (dim_t i = 0; i < D; i++)
            if (i != axis)
                product *= (1 + coordinate[i]);
        
        //std::floor is necessary since C++ rounds negative numbers up (so -0.8 becomes 0)
        return (int)std::floor(K_/product - 1);
    }
    
    std::string description() const
    {
        std::stringstream out;
        out << "HyperbolicCutShape{";
        out << "dimension: " << D << ", ";
        out << "sparsity: " << K_ << "}";
        return out.str();
    }
};

template<dim_t D>
class LimitedHyperbolicCutShape
{
private:
    double K_;
    std::array<int,D> limits_;
    
public:
    LimitedHyperbolicCutShape(double K, const std::array<int,D> &limits)
        : K_(K)
        , limits_(limits) 
    { }
    
    LimitedHyperbolicCutShape(double K, int size)
        : K_(K)
    {
        for (std::size_t d = 0; d < D; d++)
            limits_[d] = size;
    }
    
    LimitedHyperbolicCutShape(double K, std::initializer_list<int> list)
        : K_(K)
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
    
    LimitedHyperbolicCutShape(const LimitedHyperbolicCutShape &that) : K_(that.K_), limits_(that.limits_) {}
    
    LimitedHyperbolicCutShape &operator=(const LimitedHyperbolicCutShape &that)
    {
        K_ = that.K_;
        limits_ = that.limits_;
        return *this;
    }
    
    int bbox(dim_t axis) const
    {
        return std::min( limits_[axis]-1, (int)std::floor(K_-1.0));
    }
    
    template<class MultiIndex>
    int limit(const MultiIndex &coordinate, dim_t axis) const
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
    
    std::string description() const
    {
        std::stringstream out;
        out << "LimitedHyperbolicCutShape{";
        out << "dimension: " << D << ", ";
        out << "sparsity: " << K_ << ", ";
        out << "limits (exclusive): " << limits_ << "}";
        return out.str();
    }
};

}

#endif