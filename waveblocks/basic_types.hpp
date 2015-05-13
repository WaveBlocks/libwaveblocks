#ifndef WAVEBLOCKS_BASIC_TYPES_HPP
#define WAVEBLOCKS_BASIC_TYPES_HPP

#include <complex>
#include <Eigen/Core>

namespace waveblocks {

typedef double real_t;
typedef std::complex<real_t> complex_t;
typedef int dim_t;

template<int R, int C>
using CMatrix = Eigen::Matrix<complex_t,R,C>;

template<int R, int C>
using RMatrix = Eigen::Matrix<real_t,R,C>;

}

#endif