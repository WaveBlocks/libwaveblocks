#ifndef WAVEBLOCKS_HYPERCUBIC_SHAPE
#define WAVEBLOCKS_HYPERCUBIC_SHAPE

#include <string>
#include <sstream>
#include <array>
#include <stdexcept>
#include <initializer_list>

#include "basic_types.hpp"

namespace waveblocks {

/**
 * \ingroup BasicShape
 */
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
    
    /**
     * \brief Returns shape limit in one direction starting from a base node.
     * 
     * Returns for a given \e axis \f$ j \f$ and a given \e base \e node \f$ \boldsymbol{n} \f$
     * the largest element \f$ k^\star \f$ that satisfies: 
     * \f[ \boldsymbol{k} \in \mathfrak{K}, \;
     * k_i =
     *    \begin{cases}
     *       n_i,& i \neq j\\
     *       k^\star, & i = j
     *    \end{cases}
     * \f]
     * 
     * A \e base \e node is a node whose \f$ j \f$ -th entry is zero.
     * 
     * Notice that \f$ n_j \f$ does not influence return value \f$ k^\star \f$. 
     * Therefore \f$ n_j \f$ can be any value since it is simply ignored.
     * 
     * \param base_node base point \f$ \boldsymbol{n} \f$. \f$ n_j \f$ can be arbitrary.
     * \param axis \f$ j \f$
     * \return \f$ k^\star \f$
     */
    template<class MultiIndex>
    int limit(MultiIndex base_node, dim_t axis) const
    {
        { (void)(base_node); } //disable unused-parameter warning
        
        for (dim_t d = 0; d < D; d++) {
            if (d != axis && base_node[d] >= limits_[d])
                return -1;
        }
        return limits_[axis]-1;
    }
    
    /**
     * \brief Returns upper shape limit in one direction.
     * 
     * Returns for a given \e axis \f$ j \f$ a (preferably) as small as possible \e limit \f$ L_j \f$
     * such that:
     * \f[
     * \forall \boldsymbol{k} \in \mathfrak{K} \,\colon\; k_j \leq L_j
     * \f]
     * \param axis \f$ j \f$
     * \return \f$ L_j \f$
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