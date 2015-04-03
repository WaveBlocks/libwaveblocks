#ifndef WAVEBLOCKS_MULTI_INDEX_HPP
#define WAVEBLOCKS_MULTI_INDEX_HPP

#include <iostream>
#include <array>

namespace waveblocks {

template<std::size_t D>
using MultiIndex = std::array<unsigned,D>;

template<std::size_t D>
std::ostream &operator<<(std::ostream &out, const MultiIndex<D> &index)
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