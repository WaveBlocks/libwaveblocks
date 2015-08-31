#ifndef WAVEBLOCKS_HYPERCUBIC_SHAPE
#define WAVEBLOCKS_HYPERCUBIC_SHAPE

#include <string>
#include <array>
#include <stdexcept>
#include <initializer_list>

#include "basic_types.hpp"
#include "shape_base.hpp"

namespace waveblocks {

/**
 * \ingroup BasicShape
 * 
 * \brief This class implements the hypercubic basis shape.
 * 
 * A \f$ D \f$-dimensional hypercubic shape with limits \f$ \mathbf{K}=\{K_1,\dots,K_D\} \f$ is defined as the set
 * 
 * \f[
 * \mathfrak{K}(D,\mathbf{K}) := \left\{ (k_1, \dots, k_D) \in \mathbb{N}_0^D | k_i < K_i \forall i \right\}
 * \f]
 * 
 */
template<dim_t D>
class HyperCubicShape : public AbstractShape<D>
{
private:
    std::array<int,D> limits_;
    
public:
    /**
     * \param[in] limits array of all limits \f$ \{K_i\} \f$
     */
    HyperCubicShape(const std::array<int,D> &limits) 
        : limits_(limits)
    { }
    
    /**
     * \brief Set limits to \f$ K_d := K \; \forall d \f$.
     * 
     * \param limit Limit \f$ K \f$.
     */
    HyperCubicShape(int limit)
    {
        for (std::size_t d = 0; d < D; d++)
            limits_[d] = limit;
    }
    
    /**
     * \param list list of all limits \f$ \{K_i\} \f$
     */
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
    
    virtual int limit(int const* base_node, dim_t axis) const override
    {
        { (void)(base_node); } //disable unused-parameter warning
        
        for (dim_t d = 0; d < D; d++) {
            if (d != axis && base_node[d] >= limits_[d])
                return -1;
        }
        return limits_[axis]-1;
    }
    
    virtual int bbox(dim_t axis) const override
    {
        return limits_[axis]-1;
    }
    
    virtual void print(std::ostream & out) const override
    {
        out << "HyperCubicShape{ limits (exclusive): [";
        for (dim_t i = 0; i < D-1; i++) {
            out << limits_[i] << ",";
        }
        out << limits_[D-1] << "]";
        out << "}";
    }
};

}

#endif