#ifndef WAVEBLOCKS_TYPES_HPP
#define WAVEBLOCKS_TYPES_HPP

#include <array>
#include <complex>

namespace waveblocks {
typedef std::size_t dim_t;
typedef double real_t;
typedef std::complex<real_t> complex_t;
}

#include "MultiIndexUint.hpp"

namespace waveblocks {

template<dim_t D>
using MultiIndex = IntegerMultiIndex<std::size_t, D, 3>;

}

#include "PolymorphousShape.hpp"

namespace waveblocks {

template<dim_t D>
using Shape = BasicShape<D>;

}

#endif