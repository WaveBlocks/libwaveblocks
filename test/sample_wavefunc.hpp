#ifndef SAMPLE_WAVEFUNC_HPP
#define SAMPLE_WAVEFUNC_HPP

#include <memory>

#include "waveblocks/basic_types.hpp"
#include "waveblocks/hagedorn_parameter_set.hpp"
#include "waveblocks/hagedorn_coefficient_vector.hpp"

namespace waveblocks {
template<dim_t D>
std::shared_ptr< HagedornParameterSet<D> > createSampleParameters()
{
    auto params = std::make_shared< HagedornParameterSet<D> >();

    params->eps = 0.9;
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

template<dim_t D, class S>
std::shared_ptr< CoefficientVector<D,S> > createSampleCoefficients(const std::shared_ptr< const SlicedShapeEnumeration<D,S> > &enumeration)
{
    auto coeffs = std::make_shared< CoefficientVector<D,S> >(enumeration);

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

        ordinal++;
    }

    return coeffs;
}

}

#endif
