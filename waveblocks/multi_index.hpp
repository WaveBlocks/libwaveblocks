#ifndef WAVEBLOCKS_MULTI_INDEX_HPP
#define WAVEBLOCKS_MULTI_INDEX_HPP

#include <iostream>
#include <array>

#include "basic_types.hpp"

namespace waveblocks {

template<dim_t D>
using MultiIndex = std::array<unsigned,std::size_t(D)>;

template<dim_t D>
class LexicalMultiIndexCompare
{
public:
    /**
     * yields true if the first argument appears before the second
     */
    bool operator()(const MultiIndex<D> &a, const MultiIndex<D> &b)
    {
        for (dim_t i = 0; i < D; i++) {
            if (a[i] < b[i])
                return true;
            if (a[i] > b[i])
                return false;
        }
        return false; //strict weak ordering: a = b => comp(a,b) = false
    }
};
    
// This code does not work on GNU-Compiler
//
// template<dim_t D>
// std::ostream &operator<<(std::ostream &out, typename const MultiIndex<D> &index)
// {
//     std::cout << "(";
//     for (dim_t i = 0; i < D-1; i++)
//         std::cout << index[i] << ", ";
//     if (D != 0)
//         std::cout << index[D-1];
//     std::cout << ")";
//     return out;
// }

template<std::size_t D>
std::ostream &operator<<(std::ostream &out, const std::array<unsigned,D> &index)
{
    std::cout << "(";
    for (std::size_t i = 0; i < D-1; i++)
        std::cout << index[i] << ", ";
    if (D != 0)
        std::cout << index[D-1];
    std::cout << ")";
    return out;
}

}



#endif