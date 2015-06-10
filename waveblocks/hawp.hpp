#ifndef WAVEBLOCKS_HAGEDORN_WAVEPACKET
#define WAVEBLOCKS_HAGEDORN_WAVEPACKET

#include <Eigen/Core>
#include <Eigen/Dense>

#include <vector>
#include <array>
#include <memory>
#include <initializer_list>

#include "basic_types.hpp"
#include "math_util.hpp"
#include "hagedorn_parameter_set.hpp"
#include "hawp_evaluator.hpp"
#include "shape_enum.hpp"

#include "kahan_sum.hpp"

namespace waveblocks {

template<dim_t D, class MultiIndex>
class HaWpBasis
{
public:
    double eps;
    const HagedornParameterSet<D>& parameters;
    const ShapeEnum<D,MultiIndex>& enumeration;
    
    HaWpBasis(double eps, 
              const HagedornParameterSet<D>& parameters,
              const ShapeEnum<D,MultiIndex>& enumeration)
        : eps(eps)
        , parameters(parameters)
        , enumeration(enumeration)
    { }
    
    template<int N>
    HaWpEvaluator<D, MultiIndex,N> at(const Eigen::Matrix<complex_t,D,N> &x) const
    {
        return {eps, parameters, enumeration, x};
    }
    
    template<int N>
    HaWpEvaluator<D, MultiIndex,N> at(const Eigen::Matrix<real_t,D,N> &x) const
    {
        Eigen::Matrix<complex_t,D,N> x_complex = x.template cast<complex_t>();
        return {eps, parameters, enumeration, x_complex};
    }
};

namespace hawp {

template<dim_t D, class MultiIndex>
HaWpBasis<D,MultiIndex> basis(double eps, 
                              const HagedornParameterSet<D>& parameters,
                              const ShapeEnum<D,MultiIndex>& enumeration)
{
    return {eps, parameters, enumeration};
}

template<dim_t D>
complex_t prefactor(const HagedornParameterSet<D>& paramaters)
{
    return real_t(1)/paramaters.sqrt_detQ();;
}

}

}

#endif