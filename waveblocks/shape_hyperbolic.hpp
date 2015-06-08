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
    double S_;
    
public:
    HyperbolicCutShape(double K) : S_(K) {}
    
    HyperbolicCutShape(const HyperbolicCutShape &that) : S_(that.S_) {}
    
    HyperbolicCutShape &operator=(const HyperbolicCutShape &that)
    {
        S_ = that.S_;
        return *this;
    }
    
    int bbox(dim_t axis) const
    {
        (void)(axis); //suppress -Wunused-parameter
        return std::floor(S_-1.0);
    }
    
    template<class MultiIndex>
    int limit(const MultiIndex &coordinate, dim_t axis) const
    {
        double s = S_;
        
        for (dim_t i = 0; i < D; i++) {
            if (i != axis) {
                s /= 1 + coordinate[i];
            }
        }
        
        return (int)s - 1;
    }
    
    std::string description() const
    {
        std::stringstream out;
        out << "HyperbolicCutShape{";
        out << "dimension: " << D << ", ";
        out << "sparsity: " << S_ << "}";
        return out.str();
    }
};

/**
 * \f[
 * \mathfrak{K}(D,S,K) := {(k_1,\dots,k_D) | 
 *      0 \leq k_d < K_d \forall d \in [1,\dots,D] \land
 *      \displaystyle\prod_{d=1}^{D} (1+k_d) \leq S}
 * \f]
 * \tparam D number of multi-index dimensions
 */
template<dim_t D>
class LimitedHyperbolicCutShape
{
private:
    double S_;
    std::array<int,D> limits_;
    
public:
    /**
     * \param S sparsity parameter S
     * \param limits list of all limits
     */
    LimitedHyperbolicCutShape(double S, const std::array<int,D> &limits)
        : S_(S)
        , limits_(limits) 
    { }
    
    /**
     * \brief S sparsity parameter S
     * \brief 
     */
    LimitedHyperbolicCutShape(double S, int size)
        : S_(S)
    {
        for (std::size_t d = 0; d < D; d++)
            limits_[d] = size;
    }
    
    /**
     * \brief S sparsity parame S
     * \brief list of all limits
     */
    LimitedHyperbolicCutShape(double S, std::initializer_list<int> list)
        : S_(S)
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
    
    /**
     * 
     */
    int bbox(dim_t axis) const
    {
        return std::min( limits_[axis]-1, (int)S_ - 1);
    }
    
    /**
     * 
     */
    template<class MultiIndex>
    int limit(const MultiIndex &coordinate, dim_t axis) const
    {
        double s = S_;
        
        for (dim_t i = 0; i < D; i++) {
            if (i != axis) {
                if (coordinate[i] >= limits_[i])
                    return -1;
                else
                    s /= 1 + coordinate[i];
            }
        }
        
        return std::min((int)limits_[axis]-1, (int)s - 1);
    }
    
    std::string description() const
    {
        std::stringstream out;
        out << "LimitedHyperbolicCutShape{";
        out << "dimension: " << D << ", ";
        out << "sparsity: " << S_ << ", ";
        out << "limits (exclusive): " << limits_ << "}";
        return out.str();
    }
};

}

#endif