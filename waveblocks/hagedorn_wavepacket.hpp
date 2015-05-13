#ifndef WAVEBLOCKS_HAGEDORN_WAVEPACKET
#define WAVEBLOCKS_HAGEDORN_WAVEPACKET

#include <Eigen/Core>
#include <Eigen/Dense>

#include <vector>
#include <valarray>
#include <memory>
#include <initializer_list>

#include "basic_types.hpp"
#include "math_util.hpp"
#include "hagedorn_parameter_set.hpp"
#include "sliced_shape_enumeration.hpp"

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
    
    Eigen::Matrix<complex_t,C,1> coefficient(std::size_t ordinal) const
    {
        Eigen::Matrix<complex_t,C,1> coeff;
        for (dim_t c = 0; c < C; c++)
            coeff(c,0) = coefficients_[c]->operator[](ordinal);
        
        return coeff;
    }
    
    
private:
    class EvaluatorBase_
    {
    private:
        real_t eps_;
        std::shared_ptr< HagedornParameterSet<D> > parameters_;
        std::shared_ptr< ShapeEnumeration<D> > enumeration_;
        
        /**
         * precomputed expression: x - q
         */
        Eigen::Matrix<complex_t,D,1> dx;
        
        /**
         * precomputed expression: Q^{-1}
         */
        Eigen::Matrix<complex_t,D,D> Qinv;
        
        /**
         * precomputed expression: Q^H * Q^{-T}
         */
        Eigen::Matrix<complex_t,D,D> QhQinvt;
        
        /**
         * precomputed expression: Q^{-1} * (x - q)
         */
        Eigen::Matrix<complex_t,D,1> Qinv_dx;
        
        int islice;
        
        complex_t evaluateGroundState() const
        {
            complex_t pr1 = dx.transpose()*parameters_->P*parameters_->Q.inverse()*dx;
            complex_t pr2 = parameters_->p.transpose()*dx;
            
            complex_t exponent = complex_t(0.0, 1.0)/(eps_*eps_) * (0.5*pr1 + pr2);
            
            return 1.0/std::pow(pi<real_t>()*eps_*eps_, D/4.0) * std::exp(exponent);
        }
        
        Eigen::Matrix<complex_t,C,1> evaluateBasis(dim_t axis, 
                                                   MultiIndex<D> k, 
                                                   const Eigen::Matrix<complex_t,C,1> &curr_basis, 
                                                   const Eigen::Matrix<complex_t,C,D> &prev_bases) const
        {
            //compute {sqrt(k[i])*phi[k-e[i]]}
            //  e[i]: unit vector aligned to i-th axis
            Eigen::Matrix<complex_t,D,C> prev_bases_scaled = prev_bases.transpose();
            for (dim_t d = 0; d < D; d++)
                prev_bases_scaled.row(d) *= std::sqrt( real_t(k[d]) );
            
            Eigen::Matrix<complex_t,C,1> pr1 = curr_basis.transpose() * std::sqrt(2.0)/eps_ * Qinv_dx(axis,0);
            Eigen::Matrix<complex_t,C,1> pr2 = QhQinvt.row(axis)*prev_bases_scaled;
            
            return (pr1 - pr2) / std::sqrt( real_t(k[axis])+1.0);
        }
        
    public:
        EvaluatorBase_(const HagedornWavepacket<D,C> &wavepacket, const Eigen::Matrix<complex_t,D,1> &x)
            : eps_(wavepacket.eps_)
            , parameters_(wavepacket.parameters_)
            , enumeration_(wavepacket.enumeration_)
            , islice(-1)
        {
            // precompute ...
            dx = x - parameters_->q.template cast<complex_t>();
            Qinv = parameters_->Q.inverse();
            QhQinvt = parameters_->Q.adjoint()*Qinv.transpose();
            Qinv_dx = Qinv*dx;
        }
        
        /**
         * \param[in] islice index of next slice
         * \param[in] prev_slice_values
         * \param[in] curr_slice_values
         * \param[out] next_slice_values
         */
        void do_recursion(std::size_t islice,
                          typename std::vector< Eigen::Matrix<complex_t,C,1> >::const_iterator prev_slice_values,
                          typename std::vector< Eigen::Matrix<complex_t,C,1> >::const_iterator curr_slice_values,
                          typename std::vector< Eigen::Matrix<complex_t,C,1> >::iterator next_slice_values)
        {
            if (islice == 0) {
                next_slice_values[0] = Eigen::Matrix<complex_t,C,1>::Constant( evaluateGroundState() );
            } else {
                //loop over all multi-indices within next slice [j = position of multi-index within next slice]
                for (std::size_t j = 0; j < enumeration_->slice(islice).size(); j++) {
                    MultiIndex<D> next_index = enumeration_->slice(islice)[j];
                    //find valid precursor: find first non-zero entry
                    dim_t axis = D;
                    for (dim_t d = 0; d < D; d++) {
                        if (next_index[d] != 0) {
                            axis = d;
                            break;
                        }
                    }
                    
                    assert(axis != D); //assert that multi-index contains some non-zero entries
                    
                    //retrieve the basis value within current slice
                    MultiIndex<D> curr_index = next_index; curr_index[axis] -= 1; //get backward neighbour
                    std::size_t curr_ordinal = enumeration_->slice(islice-1).find(curr_index);
                    
                    assert(curr_ordinal < enumeration_->slice(islice-1).size()); //assert that multi-index has been found within current slice
                    Eigen::Matrix<complex_t,C,1> curr_basis = curr_slice_values[curr_ordinal];
                    
                    //retrieve the basis values within previous slice
                    Eigen::Matrix<complex_t,C,D> prev_bases;
                    for (dim_t d = 0; d < D; d++) {
                        if (curr_index[d] == 0) {
                            //precursor is out of shape therefore set this precursor value to zero
                            prev_bases.col(d) = Eigen::Matrix<complex_t,C,1>::Zero();
                        }
                        else {
                            MultiIndex<D> prev_index = curr_index; prev_index[d] -= 1; //get backward neighbour
                            std::size_t prev_ordinal = enumeration_->slice(islice-2).find(prev_index);
                            assert (prev_ordinal < enumeration_->slice(islice-2).size()); //assert that multi-index has been found within previous slice
                            prev_bases.col(d) = prev_slice_values[prev_ordinal];
                        }
                    }
                    
                    //compute basis value within next slice
                    Eigen::Matrix<complex_t,C,1> next_basis = evaluateBasis(axis, curr_index, curr_basis, prev_bases);
                    next_slice_values[j] = next_basis;
                }
            }
        }
    };
    
public:
    /**
     * This class evaluates one slice by one.
     * \tparam T type of quadrature/evaluation points (either real_t or complex_t)
     */
    template<class T>
    class Evaluator : public EvaluatorBase_
    {
    public:
        Evaluator(const HagedornWavepacket<D,C> &wavepacket, const Eigen::Matrix<T,D,1> &x)
            : EvaluatorBase_(wavepacket, x.template cast<complex_t>())
        { }
    };
    
private:
    template<class T>
    Eigen::Matrix<complex_t,C,1> evaluate_and_reduce_(const Eigen::Matrix<T,D,1> &x) const
    {
        //Use Kahan's algorithm to accumulate bases with O(1) numerical error instead of O(Sqrt(N))
        KahanSum< Eigen::Matrix<complex_t,C,1> > psi;
        
        Evaluator<T> evaluator(*this, x);
        
        std::vector< Eigen::Matrix<complex_t,C,1> > prev_slice_values(0);
        std::vector< Eigen::Matrix<complex_t,C,1> > curr_slice_values(0);
        std::vector< Eigen::Matrix<complex_t,C,1> > next_slice_values(1);
        
        evaluator.do_recursion(0, prev_slice_values.begin(), curr_slice_values.begin(), next_slice_values.begin());
        
        psi += coefficient(0).cwiseProduct( next_slice_values[0] );
        
        for (std::size_t islice = 1; islice < enumeration_->count_slices(); islice++) {
            //exchange slices
            std::swap(prev_slice_values, curr_slice_values);
            std::swap(curr_slice_values, next_slice_values);
            next_slice_values.resize( enumeration_->slice(islice).size() );
            
            evaluator.do_recursion(islice, prev_slice_values.begin(), curr_slice_values.begin(), next_slice_values.begin());
            
            for (std::size_t j = 0; j < next_slice_values.size(); j++) {
                psi += coefficient(enumeration_->slice(islice).offset() + j).cwiseProduct( next_slice_values[j] );
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
    Eigen::Matrix<complex_t,C,1> operator()(const Eigen::Matrix<real_t,D,1> &x) const
    {
        return evaluate_and_reduce_<real_t>(x);
    }
    
    /**
     * 
     */
    Eigen::Matrix<complex_t,C,1> operator()(const Eigen::Matrix<complex_t,D,1> &x) const
    {
        return evaluate_and_reduce_<complex_t>(x);
    }
};

}

#endif