#ifndef SAMPLE_WAVEFUNC_HPP
#define SAMPLE_WAVEFUNC_HPP

#include <vector>
#include <complex>

#include "waveblocks/basic_types.hpp"
#include "waveblocks/shape_enum.hpp"
#include "waveblocks/hagedorn_parameter_set.hpp"

namespace waveblocks {
template<dim_t D>
HagedornParameterSet<D> createSampleParameters()
{
    Eigen::Matrix<real_t,D,1> q,p;
    Eigen::Matrix<complex_t,D,D> Q,P;
    
    for (dim_t i = 0; i < D; i++) {
        q(i,0) = std::cos(i+1);
        p(i,0) = std::sin(i+1);
        
        for (dim_t j = 0; j < D; j++) {
            Q(i,j) = complex_t( 0.3*std::sin(i+j+1), 0.3*std::cos(i+j+1) );
            P(i,j) = complex_t( 0.3*std::cos(i+j+1), 0.3*std::sin(i+j+1) );
        }
        
        Q(i,i) += complex_t(1.0,0.0);
        P(i,i) += complex_t(0.0,1.0);
    }
    
    return {q,p,Q,P};
}

template<dim_t D, class MultiIndex>
std::vector<complex_t> createSampleCoefficients(const ShapeEnum<D,MultiIndex>& enumeration)
{
    std::vector<complex_t> coeffs(enumeration.n_entries());
    
    real_t norm = 0.0;
    
    std::size_t ordinal = 0;
    for (int islice = 0; islice < enumeration.n_slices(); islice++) {
        for (auto index : enumeration.slice(islice)) {
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

            norm += x*x + y*y;

            ordinal++;
        }
    }
    
    //normalize wavepacket
    for (auto & c : coeffs)
        c /= norm;
    
    return coeffs;
}

}

#endif
