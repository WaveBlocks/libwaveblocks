#ifndef WAVEBLOCKS_TEST_CHECK_COEFFICIENTS_HPP
#define WAVEBLOCKS_TEST_CHECK_COEFFICIENTS_HPP

#include <iostream>
#include <array>
#include <fstream>

#include "waveblocks/basic_types.hpp"
#include "waveblocks/hagedorn_wavepacket.hpp"

namespace waveblocks {

template<dim_t D>
void compareCoefficientsToReferenceFile(const std::array< HagedornWavepacket<D>, D> &gradient, const char *filename)
{
    auto enumeration = gradient[0].enumeration();
    
    std::cout << "compare wavepacket coefficients to reference file {" << std::endl;
    std::ifstream csv(filename);
    if (csv.fail())
        std::cout << "   [ERROR] File not found!" << std::endl;
    
    std::size_t lines = 0;
    while (csv.good()) {
        real_t tol = 1e-10;
        
        Eigen::Matrix<complex_t,D,1> ref;
        
        // read multi-index
        std::array<int,D> index;
        for (dim_t d = 0; d < D; d++) {
            int entry;
            csv >> entry;
            index[d] = entry;
        }
        
        // read real part of coefficients
        real_t coeffs_real[D];
        for (dim_t d = 0; d < D; d++) {
            real_t coeff_real;
            csv >> coeff_real;
            coeffs_real[d] = coeff_real;
        }
        
        // read imaginary part of coefficients
        for (dim_t d = 0; d < D; d++) {
            real_t coeff_imag;
            csv >> coeff_imag;
            ref(d,0) = complex_t(coeffs_real[d],coeff_imag);
        }
        
        if (csv.good()) {
            ++lines;
            
            int islice = std::accumulate(index.begin(), index.end(), int(0));
            
            if (enumeration->contains(index)) {
                Eigen::Matrix<complex_t,D,1> coeff;
                for (int d = 0; d < D; d++) {
                    coeff(d,0) = gradient[d].coefficients()[ enumeration->slice(islice).offset() + enumeration->slice(islice).find(index) ];
                }
                
                real_t error = (coeff - ref).norm()/ref.norm();
                
                if ( error > tol) {
                    std::cout << "   [FAILURE] at line " << lines << ": coefficient at node " << 
                            index << " does not match entry in csv file. error = " << error << std::endl;
                    std::cout << "     computed value: " << coeff.transpose() << std::endl;
                    std::cout << "     reference value: " << ref.transpose() << std::endl;
                }
            }
            else {
                if (ref.norm() > tol) {
                    std::cout << "   [FAILURE] at line " << lines << ": reference file contains unexpected node: " << index << std::endl;
                    continue;
                }
            }
        }
    }
    std::cout << "   [INFO] lines processed: " << lines << std::endl;
    std::cout << "}" << std::endl;
}

}

#endif
