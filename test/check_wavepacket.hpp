#ifndef WAVEBLOCKS_TEST_CHECK_WAVEPACHET_HPP
#define WAVEBLOCKS_TEST_CHECK_WAVEPACHET_HPP

#include <iostream>
#include <fstream>
#include <array>

#include "waveblocks/basic_types.hpp"
#include "waveblocks/hawp.hpp"

namespace waveblocks {

/**
 * Compares wavepacket to reference value contained in a space delimited csv files
 * with columns: x; real part of psi; imaginary part of psi.
 * \param wavepacket
 * \param filename
 */
template<dim_t D, class MultiIndex>
void compareWavepacketToReferenceFile(double eps,
                                      const HaWpParamSet<D>& parameters,
                                      const ShapeEnum<D,MultiIndex>& enumeration,
                                      const std::vector<complex_t>& coefficients,
                                      const char *filename)
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
        real_t ref_real, ref_imag;
        in >> ref_real;
        in >> ref_imag;
        
        complex_t ref = complex_t(ref_real, ref_imag);
        
        //compute wavepacket value
        if (in.good()) {
            ++lines;
            
            complex_t psi = hawp::basis(eps,&parameters,&enumeration).at(x).reduce(coefficients)(0,0);
            
            auto error = std::norm(psi - ref)/std::norm(ref);
            
            if (error > 10e-10) {
                std::cout << "      [FAILURE] mismatch at line " << lines << ". error = " << error << std::endl;
            }
        }
    }
    
    std::cout << "      [INFO] processed " << lines << " lines" << std::endl;
    std::cout << "}" << std::endl;
}

// template<dim_t D, int C>
// void compareWavepacketToReferenceFile(const HagedornWavepacket<D, C> &wavepacket, const char *filename)
// {
//     std::cout << "compare wavepacket to reference files {" << std::endl;
//     std::cout << "   [FILE] " << filename << std::endl;
//     std::ifstream in(filename);
// 
//     if (!in.good())
//         std::cout << "   [ERROR] File not found!" << std::endl;
// 
//     std::size_t lines = 0;
//     while (in.good()) {
//         //read position
//         Eigen::Matrix<real_t,D,1> x;
//         for (dim_t d = 0; d < D; d++)
//             in >> x(d,0);
// 
//         //read reference values
//         Eigen::Matrix<real_t,C,1> ref_real;
//         for (int c = 0; c < C; c++) {
//             real_t val;
//             in >> val;
//             ref_real(c,0) = val;
//         }
// 
//         Eigen::Matrix<real_t,C,1> ref_imag;
//         for (int c = 0; c < C; c++) {
//             real_t val;
//             in >> val;
//             ref_imag(c,0) = val;
//         }
// 
//         Eigen::Matrix<complex_t,C,1> ref;
//         for (int c = 0; c < C; c++) {
//             ref(c,0) = complex_t(ref_real(c,0), ref_imag(c,0));
//         }
// 
//         //compute wavepacket value
//         if (in.good()) {
//             ++lines;
// 
//             Eigen::Matrix<complex_t,C,1> psi = wavepacket(x);
// 
//             auto error = (psi - ref).norm()/ref.norm();
//             
//             if (error > 10e-10) {
//                 std::cout << "      [FAILURE] mismatch at line " << lines << ". error = " << error << std::endl;
//             }
//         }
//     }
// 
//     std::cout << "      [INFO] processed " << lines << " lines" << std::endl;
//     std::cout << "}" << std::endl;
// }

}

#endif
