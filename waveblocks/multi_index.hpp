#ifndef WAVEBLOCKS_MULTI_INDEX_HPP
#define WAVEBLOCKS_MULTI_INDEX_HPP

#include <iostream>
#include <array>

#include "basic_types.hpp"

#include "tiny_multi_index.hpp"

namespace waveblocks {

template<dim_t D>
using MultiIndex = TinyMultiIndex<std::size_t,D>;

template<dim_t D>
std::ostream &operator<<(std::ostream &out, const MultiIndex<D> &index)
{
    std::cout << "(";
    for (dim_t i = 0; i < D-1; i++)
        std::cout << index[i] << ", ";
    if (D != 0)
        std::cout << index[D-1];
    std::cout << ")";
    return out;
}

}

#endif