#ifndef WAVEBLOCKS_SHAPE_ENUMERATION_OPERATIONS_HPP
#define WAVEBLOCKS_SHAPE_ENUMERATION_OPERATIONS_HPP

#include <vector>
#include <algorithm>

#include "waveblocks/shape_enum.hpp"

namespace waveblocks {

namespace shape_enum {

/**
 * \brief constructs union of two slices
 * 
 * \f[ \forall k \colon k \in union(\mathfrak{K}_1,\mathfrak{K}_2) \iff k \in \mathfrak{K}_1 \land k \in \mathfrak{K}_2 \f]
 * 
 * \param slice1 first slice
 * \param slice2 second slice
 * \return union of both slices
 */
template<class MultiIndex>
std::vector<MultiIndex> strict_union(const std::vector<MultiIndex>& slice1, 
                                     const std::vector<MultiIndex>& slice2)
{
    std::less<MultiIndex> less;
    
    auto it1 = slice1.cbegin();
    auto it2 = slice2.cbegin();
    
    std::vector<MultiIndex> slice12(slice1.size() + slice2.size());
    
    auto it12 = slice12.begin();
    
    while (it1 != slice1.end() && it2 != slice2.end()) {
        if (*it1 == *it2) {
            // lhs = rhs
            *it12 = *it1;
            ++it1;
            ++it2;
            ++it12;
        }
        else if ( less(*it1, *it2) ) {
            // lhs < rhs
            *it12 = *it1;
            ++it1;
            ++it12;
        }
        else {
            // lhs > rhs
            *it12 = *it2;
            ++it2;
            ++it12;
        }
    }
    
    it12 = std::copy(it1, slice1.end(), it12);
    it12 = std::copy(it2, slice2.end(), it12);
    
    slice12.resize(it12 - slice12.begin());
    slice12.shrink_to_fit();
    
    return slice12;
}

}

}

#endif