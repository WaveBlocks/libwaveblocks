#ifndef SAMPLE_WAVEFUNC_HPP
#define SAMPLE_WAVEFUNC_HPP

#include "waveblocks/basic_types.hpp"

namespace waveblocks {
    
    template<dim_t D>
    void createSampleParameters(HagedornParameterSet<D> &params)
    {
        params.eps = 0.9;
        for (dim_t i = 0; i < D; i++) {
            params.q(i,0) = std::cos(i+1);
            params.p(i,0) = std::sin(i+1);
            
            for (dim_t j = 0; j < D; j++) {
                params.Q(i,j) += complex_t( 0.3*std::sin(i+j+1), 0.3*std::cos(i+j+1) );
                params.P(i,j) += complex_t( 0.3*std::cos(i+j+1), 0.3*std::sin(i+j+1) );
            }
        }
    }
    
    template<dim_t D, class S>
    void createSampleCoefficients(const SlicedShapeEnumeration<D,S> &enumeration,
                                  std::vector<complex_t> &coeffs)
    {
        coeffs.resize(enumeration.size());
        
        std::size_t ordinal = 0;
        for (auto index : enumeration) {
            real_t falloff = 0.1;
            
            real_t x = 0.0;
            real_t y = 0.0;
            
            int sum = 0;
            for (dim_t d = 0; d < D; d++) {
                x += std::sin(index[d] + 1.0*(d+1)/real_t(D));
                y += std::cos(index[d] + 1.5*(d+1)/real_t(D));
                sum += index[d];
            }
            
            x *= std::exp(-falloff*sum);
            y *= std::exp(-falloff*sum);
            
            coeffs[ordinal] = complex_t(x,y);
            
            //std::cout << index << ": " << coeffs[ordinal] << std::endl;
            
            ordinal++;
        }
    }
    
}

#endif