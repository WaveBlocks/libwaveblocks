#ifndef WAVEBLOCKS_TEST_CHECK_WAVEPACHET_HPP
#define WAVEBLOCKS_TEST_CHECK_WAVEPACHET_HPP

#include <iostream>
#include <fstream>

#include "waveblocks/basic_types.hpp"
#include "waveblocks/multi_index.hpp"
#include "waveblocks/hagedorn_wavepacket.hpp"

namespace waveblocks {

/**
 * Compares wavepacket to reference value contained in a space delimited csv files
 * with columns: x; real part of psi; imaginary part of psi.
 * \param wavepacket
 * \param filename
 */
template<dim_t D, int V>
void compareWavepacketToReferenceFile(const HagedornWavepacket<D,V> &wavepacket, const char *filename)
{
    std::cout << "compare wavepacket to reference files {" << std::endl;
    std::cout << "   [FILE] " << filename << std::endl;
    std::ifstream in(filename);

    if (!in.good())
        std::cout << "   [ERROR] File not found!" << std::endl;

    std::size_t lines = 0;
    while (in.good()) {
        //read position
        Eigen::Matrix<real_t,D,1> x;
        for (dim_t d = 0; d < D; d++)
            in >> x(d,0);

        //read reference values
        Eigen::Matrix<real_t,V,1> ref_real;
        for (int v = 0; v < V; v++) {
            real_t val;
            in >> val;
            ref_real(v,0) = val;
        }

        Eigen::Matrix<real_t,V,1> ref_imag;
        for (int v = 0; v < V; v++) {
            real_t val;
            in >> val;
            ref_imag(v,0) = val;
        }

        Eigen::Matrix<complex_t,V,1> ref;
        for (int v = 0; v < V; v++) {
            ref(v,0) = complex_t(ref_real(v,0), ref_imag(v,0));
        }

        //compute wavepacket value
        if (in.good()) {
            ++lines;

            Eigen::Matrix<complex_t,V,1> psi = wavepacket(x);

            auto error = (psi - ref).norm()/ref.norm();

            if (error > 10e-10) {
                std::cout << "      [FAILURE] mismatch at line " << lines << ". error = " << error << std::endl;
            }
        }
    }

    std::cout << "      [INFO] processed " << lines << " lines" << std::endl;
    std::cout << "}" << std::endl;
}

}

#endif
