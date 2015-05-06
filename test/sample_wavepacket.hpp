#ifndef SAMPLE_WAVEFUNC_HPP
#define SAMPLE_WAVEFUNC_HPP

#include <memory>
#include <valarray>

#include "waveblocks/basic_types.hpp"
#include "waveblocks/sliced_shape_enumeration.hpp"
#include "waveblocks/hagedorn_parameter_set.hpp"

namespace waveblocks {
template<dim_t D>
std::shared_ptr< HagedornParameterSet<D> > createSampleParameters()
{
    auto params = std::make_shared< HagedornParameterSet<D> >();
    
    for (dim_t i = 0; i < D; i++) {
        params->q(i,0) = std::cos(i+1);
        params->p(i,0) = std::sin(i+1);

        for (dim_t j = 0; j < D; j++) {
            params->Q(i,j) += complex_t( 0.3*std::sin(i+j+1), 0.3*std::cos(i+j+1) );
            params->P(i,j) += complex_t( 0.3*std::cos(i+j+1), 0.3*std::sin(i+j+1) );
        }
    }

    return params;
}

template<dim_t D>
std::shared_ptr< std::valarray<complex_t> > createSampleCoefficients(const std::shared_ptr< const ShapeEnumeration<D> > &enumeration)
{
    auto coeffs = std::make_shared< std::valarray<complex_t> >(enumeration->size());
    
    real_t norm = 0.0;
    
    std::size_t ordinal = 0;
    for (auto index : *enumeration) {
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
        
        (*coeffs)[ordinal] = complex_t(x,y);
        
        //std::cout << index << ": " << coeffs[ordinal] << std::endl;
        
        norm += x*x + y*y;
        
        ordinal++;
    }
    
    //normalize wavepacket
    for (auto & c : *coeffs)
        c /= norm;

    return coeffs;
}

}

#endif
