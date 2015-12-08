#pragma once

#include <iostream>
#include <array>


namespace std {
    template<class T, std::size_t D>
    std::ostream &operator<<(std::ostream &out, const std::array<T,D> &value)
    {
        out << '{';
        for (std::size_t i = 0; i < D; i++) {
            if (i != 0)
                out << ", ";
            out << value[i];
        }
        out << '}';
        return out;
    }
}
