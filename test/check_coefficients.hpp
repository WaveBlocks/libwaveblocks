#ifndef WAVEBLOCKS_TEST_CHECK_COEFFICIENTS_HPP
#define WAVEBLOCKS_TEST_CHECK_COEFFICIENTS_HPP

#include <iostream>
#include <fstream>

#include "waveblocks/basic_types.hpp"
#include "waveblocks/multi_index.hpp"
#include "waveblocks/hagedorn_wavepacket.hpp"

namespace waveblocks {

template<dim_t D, int C>
void compareCoefficientsToReferenceFile(const HagedornWavepacket<D,C> &wavepacket, const char *filename)
{
    auto enumeration = wavepacket.enumeration();

    std::cout << "compare wavepacket coefficients to reference file {" << std::endl;
    std::ifstream csv(filename);
    if (csv.fail())
        std::cout << "   [ERROR] File not found!" << std::endl;

    std::size_t lines = 0;
    while (csv.good()) {
        real_t tol = 1e-10;

        Eigen::Matrix<complex_t,C,1> ref;

        //read index
        MultiIndex<C> index;
        for (dim_t c = 0; c < C; c++) {
            int entry;
            csv >> entry;
            index[c] = entry;
        }

        //read real gradient part
        real_t temp[C];
        for (dim_t c = 0; c < C; c++) {
            real_t value;
            csv >> value;
            temp[c] = value;
        }

        //read imaginary gradient part
        for (dim_t c = 0; c < C; c++) {
            real_t value;
            csv >> value;
            ref(c,0) = complex_t(temp[c],value);
        }

        if (csv.good()) {
            ++lines;

            if (!enumeration->contains(index) && ref.norm() > tol) {
                std::cout << "   [FAILURE] missing a node in our shape enumeration: " << index << std::endl;
                continue;
            }

            Eigen::Matrix<complex_t,C,1> coeff = wavepacket.coefficient( enumeration->find(index) );

            real_t error = (coeff - ref).norm()/ref.norm();

            if ( error > tol) {
                std::cout << "   [FAILURE] coefficient at node " << index << " does not match entry in csv file. error = " << error << std::endl;
                std::cout << "     computed value: " << coeff.transpose() << std::endl;
                std::cout << "     reference value: " << ref.transpose() << std::endl;
            }
        }
    }
    std::cout << "   [INFO] lines processed: " << lines << std::endl;
    std::cout << "}" << std::endl;
}

}

#endif
