#ifndef WAVEBLOCKS_HAGEDORN_WAVEPACKET
#define WAVEBLOCKS_HAGEDORN_WAVEPACKET

namespace waveblocks {

template<std::size_t D>
struct HagedornParameters
{
    Eigen::Matrix<real_t, D,1> q, p;
    Eigen::Matrix<complex_t, D,D> Q, P;
};



}

#endif