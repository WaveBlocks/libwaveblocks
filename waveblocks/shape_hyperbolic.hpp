#ifndef HYPERBOLIC_SHAPE_HPP
#define HYPERBOLIC_SHAPE_HPP

#include <cmath>
#include <string>
#include <sstream>

#include "basic_types.hpp"

namespace waveblocks {

/**
 * \ingroup BasicShape
 * 
 * This class implements the hyperbolic cut basis shape which is a special 
 * type of sparse basis set. A basis shape is essentially all information 
 * and operations related to the set \f$ \mathfrak{K} \f$ of multi-indices \f$ \boldsymbol{k} \f$. 
 * The hyperbolic cut shape in \f$ D \f$ dimensions with \e sparsity \f$S\f$ and 
 * \e limits \f$ \boldsymbol{K} = (K_1,\ldots,K_D) \f$ is defined as the set
 * 
 * \f[
 * \mathfrak{K}(D,S,\boldsymbol{K}) := \left\{(k_1,\dots,k_D) \mid 
 *      \displaystyle\prod_{d=1}^{D} (1+k_d) \leq S \right\}
 * \f]
 * 
 * \tparam D number of dimensions
 */
template<dim_t D>
class HyperbolicCutShape
{
private:
    int S_;
    
public:
    /**
     * \brief General constructor to set sparsity parameter \f$ S \f$
     * 
     * \param S sparsity parameter \f$ S \f$
     */
    HyperbolicCutShape(int S) : S_(S) {}
    
    /**
     * \sa HyperCubicShape#bbox
     */
    int bbox(dim_t axis) const
    {
        (void)(axis); //suppress -Wunused-parameter
        return S_ - 1;
    }
    
    
    
    /**
     * \sa HyperCubicShape#limit
     */
    template<class MultiIndex>
    int limit(const MultiIndex &base_node, dim_t axis) const
    {
        double s = S_;
        
        for (dim_t i = 0; i < D; i++) {
            if (i != axis) {
                s /= 1 + base_node[i];
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
 * \ingroup BasicShape
 * 
 * This class implements the hyperbolic cut basis shape which is a special 
 * type of sparse basis set. A basis shape is essentially all information 
 * and operations related to the set \f$ \mathfrak{K} \f$ of multi-indices \f$ \boldsymbol{k} \f$. 
 * The limited hyperbolic cut shape in \f$ D \f$ dimensions with \e sparsity \f$S\f$ and 
 * \e limits \f$ \boldsymbol{K} = (K_1,\ldots,K_D) \f$ is defined as the set
 * 
 * \f[
 * \mathfrak{K}(D,S,\boldsymbol{K}) := \left\{(k_1,\dots,k_D) \mid 
 *      0 \leq k_d < K_d \; \forall d \in \{ 1,\ldots,D \} \land
 *      \displaystyle\prod_{d=1}^{D} (1+k_d) \leq S \right\}
 * \f]
 * 
 * \tparam D number of dimensions
 */
template<dim_t D>
class LimitedHyperbolicCutShape
{
private:
    int S_;
    std::array<int,D> limits_;
    
public:
    /**
     * \brief General constructor to define sparsity parameter and limits.
     * 
     * \param S sparsity parameter \f$ S \f$ 
     * \param limits tuple of all limits \f$ \boldsymbol{K} \f$ 
     */
    LimitedHyperbolicCutShape(int S, const std::array<int,D> &limits)
        : S_(S)
        , limits_(limits) 
    { }
    
    /**
     * \brief Specialized constructor to set all limits \f$ K_d \f$ to the same value \f$ K^\star \f$.
     * 
     * \param S sparsity parameter \f$ S \f$ 
     * \param size \f$ K^\star \f$
     */
    LimitedHyperbolicCutShape(int S, int size)
        : S_(S)
    {
        for (std::size_t d = 0; d < D; d++)
            limits_[d] = size;
    }
    
    /**
     * \brief General constructor to define sparsity parameter and limits.
     * 
     * \param S sparsity parameter \f$ S \f$ 
     * \param list list of all limits \f$ \boldsymbol{K} \f$ 
     */
    LimitedHyperbolicCutShape(int S, std::initializer_list<int> list)
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
     * \sa HyperCubicShape#bbox
     */
    int bbox(dim_t axis) const
    {
        return std::min( limits_[axis]-1, S_ - 1);
    }
    
    /**
     * \sa HyperCubicShape#limit
     */
    template<class MultiIndex>
    int limit(const MultiIndex &base_node, dim_t axis) const
    {
        double s = S_;
        
        for (dim_t i = 0; i < D; i++) {
            if (i != axis) {
                if (base_node[i] >= limits_[i])
                    return -1;
                else
                    s /= 1 + base_node[i];
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