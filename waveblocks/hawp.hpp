#ifndef WAVEBLOCKS_HAGEDORN_WAVEPACKET
#define WAVEBLOCKS_HAGEDORN_WAVEPACKET

#include <Eigen/Core>
#include <Eigen/Dense>

#include <vector>
#include <array>
#include <memory>
#include <initializer_list>
#include <memory>

#include "basic_types.hpp"
#include "math_util.hpp"
#include "hawp_paramset.hpp"
#include "shape_enum.hpp"

#include "hawp_evaluator.hpp"
#include "hawp_gradient_evaluator.hpp"

namespace waveblocks {

/**
 * \brief represents a hagedorn wavepacket basis
 * 
 */
template<dim_t D, class MultiIndex>
struct HaWpBasis
{
public:
    double eps;
    const HaWpParamSet<D>* parameters;
    const ShapeEnum<D,MultiIndex>* enumeration;
    
    HaWpBasis(double eps, 
              const HaWpParamSet<D>* parameters,
              const ShapeEnum<D,MultiIndex>* enumeration)
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
        Eigen::Matrix<complex_t,D,N> xtmp = x.template cast<complex_t>();
        return {eps, parameters, enumeration, xtmp};
    }
};

/**
 * \brief represents a hagedorn wavepacket
 * 
 */
template<dim_t D, class MultiIndex>
struct HaWp
{
public:
    HaWpBasis<D, MultiIndex> basis;
    const std::vector<complex_t>* coefficients;
    
    HaWp(double eps, 
         const HaWpParamSet<D>* parameters,
         const ShapeEnum<D,MultiIndex>* enumeration,
         const std::vector<complex_t>* coefficients)
        : basis(eps, parameters, enumeration)
        , coefficients(coefficients)
    { }
    
    HaWp(const HaWpBasis<D,MultiIndex>& basis,
         const std::vector<complex_t>* coefficients)
        : basis(basis)
        , coefficients(coefficients)
    { }
};

namespace hawp {

template<dim_t D, class MultiIndex>
HaWpBasis<D,MultiIndex> basis(double eps, 
                              const HaWpParamSet<D>* parameters,
                              const ShapeEnum<D,MultiIndex>* enumeration)
{
    return {eps, parameters, enumeration};
}

template<dim_t D>
complex_t prefactor(const HaWpParamSet<D>& paramaters)
{
    return real_t(1)/paramaters.sqrt_detQ();;
}

template<dim_t D, class MultiIndex>
HaWpGradientEvaluator<D,MultiIndex> nabla(double eps, 
                                     const HaWpParamSet<D>* parameters,
                                     const ShapeEnum<D,MultiIndex>* base_enum,
                                     const ShapeEnum<D,MultiIndex>* grad_enum)
{
    return {eps, parameters, base_enum, grad_enum};
}

template<dim_t D, class MultiIndex>
HaWpGradientEvaluator<D,MultiIndex> nabla(const HaWpBasis<D,MultiIndex>& basis,
                                     const ShapeEnum<D,MultiIndex>* grad_enum)
{
    return {basis.eps, basis.parameters, basis.enumeration, grad_enum};
}

} // namespace hawp

} // namespace waveblocks

#endif