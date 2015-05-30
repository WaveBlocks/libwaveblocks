#ifndef WAVEBLOCKS_HAGEDORN_WAVEPACKET
#define WAVEBLOCKS_HAGEDORN_WAVEPACKET

#include <Eigen/Core>
#include <Eigen/Dense>

#include <vector>
#include <array>
#include <valarray>
#include <memory>
#include <initializer_list>

#include "basic_types.hpp"
#include "math_util.hpp"
#include "hagedorn_parameter_set.hpp"
#include "hagedorn_basis_evaluator.hpp"
#include "shape_enumeration_base.hpp"

#include "kahan_sum.hpp"

namespace waveblocks {

/**
 * \tparam D degrees of freedom
 * \tparam C number of components if this wavepacket is vector-valued
 * \parblock
 * Components of vector-valued wavepackets share the same paramaters and enumeration,
 * but use different coefficients.
 * \endparblock
 */
template<dim_t D, dim_t C = 1>
class HagedornWavepacket
{
public:
    typedef Eigen::Matrix<complex_t,D,1> CMatrixD1;
    typedef Eigen::Matrix<real_t,D,1> RMatrixD1;
    
    /**
     * complex square matrix
     * size: (#dimensions)
     */
    typedef Eigen::Matrix<complex_t,D,D> CMatrixDD;
    
    /**
     * complex column vector
     * size: (#dimensions)
     */
    typedef Eigen::Matrix<complex_t,C,1> CMatrixC1;
    
private:
    real_t eps_;
    std::shared_ptr< HagedornParameterSet<D> > parameters_;
    std::shared_ptr< ShapeEnumeration<D> > enumeration_;
    std::array< std::shared_ptr< std::valarray<complex_t> >, C> coefficients_;
    
public:
    HagedornWavepacket(real_t eps,
                       std::shared_ptr< HagedornParameterSet<D> > parameters, 
                       std::shared_ptr< ShapeEnumeration<D> > enumeration,
                       std::initializer_list< std::shared_ptr< std::valarray<complex_t> > > coefficients)
        : eps_(eps)
        , parameters_(parameters)
        , enumeration_(enumeration)
        , coefficients_()
    { 
        dim_t c = 0;
        for (auto & item : coefficients)
            coefficients_[c++] = item;
        assert (c == C);
    }
    
    HagedornWavepacket(real_t eps,
                       std::shared_ptr< HagedornParameterSet<D> > parameters, 
                       std::shared_ptr< ShapeEnumeration<D> > enumeration,
                       std::array< std::shared_ptr< std::valarray<complex_t> >, C> coefficients)
        : eps_(eps)
        , parameters_(parameters)
        , enumeration_(enumeration)
        , coefficients_(coefficients)
    { }
    
    HagedornWavepacket(real_t eps,
                       std::shared_ptr< HagedornParameterSet<D> > parameters, 
                       std::shared_ptr< ShapeEnumeration<D> > enumeration)
        : eps_(eps)
        , parameters_(parameters)
        , enumeration_(enumeration)
        , coefficients_()
    {
        for (dim_t c = 0; c < C; c++)
            coefficients_[c] = std::make_shared< std::valarray<complex_t> >(enumeration->size());
    }
    
    HagedornWavepacket(const HagedornWavepacket<D,C> &other)
        : eps_(other.eps_)
        , parameters_(other.parameters_)
        , enumeration_(other.enumeration_)
        , coefficients_(other.coefficients_)
    { }
    
    HagedornWavepacket<D,C> &operator=(const HagedornWavepacket<D,C> &other)
    {
        eps_ = other.eps_;
        parameters_ = other.parameters_;
        enumeration_ = other.enumeration_;
        coefficients_ = other.coefficients_;
        
        return *this;
    }
    
    HagedornWavepacket<D,1> component(dim_t component) const
    {
        HagedornWavepacket<D,1> slice(parameters_, enumeration_, coefficients_[component]);
        
        return slice;
    }
    
    real_t eps() const
    {
        return eps_;
    }
    
    std::shared_ptr< HagedornParameterSet<D> > parameters() const
    {
        return parameters_;
    }
    
    std::shared_ptr< ShapeEnumeration<D> > enumeration() const
    {
        return enumeration_;
    }
    
    std::shared_ptr< std::valarray<complex_t> > coefficients(dim_t component) const
    {
        return coefficients_[component];
    }
    
    std::array< std::shared_ptr< std::valarray<complex_t> >, C> coefficients() const
    {
        return coefficients_;
    }
    
    CMatrixC1 coefficient(std::size_t ordinal) const
    {
        Eigen::Matrix<complex_t,C,1> coeff;
        for (dim_t c = 0; c < C; c++)
            coeff(c,0) = coefficients_[c]->operator[](ordinal);
        
        return coeff;
    }
    
    template<int N>
    Eigen::Matrix<complex_t,C,N> evaluate_and_reduce_(const Eigen::Matrix<complex_t,D,N> &x) const
    {
        //Determine number of quadrature points. We cannot use template parameter N for
        //that as N could be 'Eigen::Dynamic'.
        int npts = x.cols();

        //Use Kahan's algorithm to accumulate bases with O(1) numerical error instead of O(Sqrt(N))
        KahanSum< Eigen::Matrix<complex_t,C,N> > psi( Eigen::Matrix<complex_t,C,N>::Zero(C,npts) );

        Evaluator<D,C,N> evaluator(eps_, parameters_, enumeration_, x);

        Eigen::Array<complex_t,Eigen::Dynamic,N> prev_slice_values(0,npts);
        Eigen::Array<complex_t,Eigen::Dynamic,N> curr_slice_values(0,npts);
        Eigen::Array<complex_t,Eigen::Dynamic,N> next_slice_values(1,npts);

        evaluator.do_recursion(0, prev_slice_values, curr_slice_values, next_slice_values);

        psi += coefficient(0)*next_slice_values.row(0).matrix();

        for (std::size_t islice = 1; islice < enumeration_->count_slices(); islice++) {
            //exchange slices
            std::swap(prev_slice_values, curr_slice_values);
            std::swap(curr_slice_values, next_slice_values);
            next_slice_values = Eigen::Array<complex_t,Eigen::Dynamic,N>( enumeration_->slice(islice).size(), npts );

            evaluator.do_recursion(islice, prev_slice_values, curr_slice_values, next_slice_values);

            for (long j = 0; j < next_slice_values.rows(); j++) {
                CMatrixC1 cj = coefficient(enumeration_->slice(islice).offset() + j);
                psi += cj*next_slice_values.row(j).matrix();
            }
        }

        return psi();
    }
    
public:
    /**
     * 
     */
    complex_t prefactor() const
    {
        return real_t(1)/parameters_->sqrt_detQ();
    }
    
    /**
     * 
     */
    template<int N>
    Eigen::Matrix<complex_t,C,N> operator()(const Eigen::Matrix<real_t,D,N> &x) const
    {
        Eigen::Matrix<complex_t,D,N> xtmp = x.template cast<complex_t>();
        return evaluate_and_reduce_<N>(xtmp);
    }
    
    /**
     * 
     */
    template<int N>
    Eigen::Matrix<complex_t,C,N> operator()(const Eigen::Matrix<complex_t,D,N> &x) const
    {
        return evaluate_and_reduce_<N>(x);
    }
};

}

#endif