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
    
    typedef Eigen::Matrix<complex_t,D,1> CMatrixD1;
    typedef Eigen::Matrix<real_t,D,1> RMatrixD1;
    
    /**
     * matrix size: (#dimensions, #dimensions)
     */
    typedef Eigen::Matrix<complex_t,D,D> CMatrixDD;
    
    typedef Eigen::Matrix<complex_t,C,1> CMatrixC1;
    
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
    
    /**
     * \tparam N number of quadrature points (if unknown: use Eigen::Dynamic)
     */
    template<int N>
    class Evaluator
    {
    public:
        /**
         * array size: (1, #quadrature points)
         */
        typedef Eigen::Array<complex_t,1,N> CArray1N;
        
        /**
         * matrix size: (#dimensions, #quadrature points)
         */
        typedef Eigen::Matrix<complex_t,D,N> CMatrixDN;
        typedef Eigen::Matrix<real_t,D,N> RMatrixDN;
        
        /**
         * matrix size: (/unspecified/, #quadrature points)
         */
        typedef Eigen::Array<complex_t,Eigen::Dynamic,N> CArrayXN;
        
        /**
         * matrix size: (#components, #quadrature points)
         */
        typedef Eigen::Matrix<complex_t,C,N> CMatrixCN;
        
    private:
        real_t eps_;
        std::shared_ptr< HagedornParameterSet<D> > parameters_;
        std::shared_ptr< ShapeEnumeration<D> > enumeration_;
        
        /**
         * number of quadrature points
         */
        int npts_;
        
        /**
         * precomputed expression: x - q
         */
        CMatrixDN dx_;
        
        /**
         * precomputed expression: Q^{-1}
         */
        CMatrixDD Qinv_;
        
        /**
         * precomputed expression: Q^H * Q^{-T}
         */
        CMatrixDD Qh_Qinvt_;
        
        /**
         * precomputed expression: Q^{-1} * (x - q)
         */
        CMatrixDN Qinv_dx_;
        
        /**
         * lookup-table for sqrt
         */
        std::vector<real_t> sqrt_;
        
        /**
         * Complexity: theta(N*D)
         */
        CArray1N evaluateGroundState() const
        {
            auto & P = parameters_->P;
            auto & Q = parameters_->Q;
            auto & p = parameters_->p;
            
            CMatrixDN P_Qinv_dx = P*Qinv_dx_;
            
            CArray1N pr1 = ( dx_.array() * P_Qinv_dx.array() ).colwise().sum();
            CArray1N pr2 = ( p.transpose()*dx_ ).array();
            
            CArray1N e = complex_t(0.0, 1.0)/(eps_*eps_) * (0.5*pr1 + pr2);
            
            return e.exp() / std::pow(pi<real_t>()*eps_*eps_, D/4.0);
        }
        
        /**
         * Complexity: theta(N*D*C)
         */
        CArray1N evaluateBasis(dim_t axis, 
                                MultiIndex<D> k, 
                                const CArray1N &curr_basis, 
                                const std::array< CArray1N, D > &prev_bases) const
        {
            CArray1N pr1 = curr_basis * Qinv_dx_.row(axis).array() * std::sqrt(2.0)/eps_ ;
            
            CArray1N pr2 = CArray1N::Zero(1,npts_);
            for (dim_t d = 0; d < D; d++) {
                pr2 += prev_bases[d] * Qh_Qinvt_(axis,d) * sqrt_[ k[d] ];
            }
            
            return (pr1 - pr2) / sqrt_[ 1+k[axis] ];
        }
        
    public:
        Evaluator(const HagedornWavepacket<D,C> &wavepacket, const CMatrixDN &x)
            : eps_(wavepacket.eps_)
            , parameters_(wavepacket.parameters_)
            , enumeration_(wavepacket.enumeration_)
            , npts_(x.cols())
            , sqrt_()
        {
            auto & q = parameters_->q;
            auto & Q = parameters_->Q;
            
            // precompute ...
            dx_ = x.colwise() - q.template cast<complex_t>();
            Qinv_ = Q.inverse();
            Qh_Qinvt_ = Q.adjoint()*Qinv_.transpose();
            Qinv_dx_ = Qinv_*dx_;
            
            //precompute sqrt lookup table
            {
                MultiIndex<D> limits = enumeration_->limits();
                int limit = 0;
                for (dim_t d = 0; d < D; d++)
                    limit = std::max(limit, (int)limits[d]);
                
                for (int i = 0; i <= limit+1; i++)
                    sqrt_.push_back( std::sqrt( real_t(i) ) );
            }
        }
        
        /**
         * Complexity: theta(S*N*D*C)
         * 
         * \param[in] islice index of next slice
         * \param[in] prev_slice_values
         * \param[in] curr_slice_values
         * \param[out] next_slice_values
         */
        void do_recursion(std::size_t islice,
                          const CArrayXN &prev_slice_values,
                          const CArrayXN &curr_slice_values,
                          CArrayXN &next_slice_values)
        {
            assert(prev_slice_values.cols() == npts_);
            assert(curr_slice_values.cols() == npts_);
            assert(next_slice_values.cols() == npts_);
            
            if (islice == 0) {
                next_slice_values.row(0) = evaluateGroundState();
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
                    CArray1N curr_basis = curr_slice_values.row(curr_ordinal);
                    
                    //retrieve the basis values within previous slice
                    std::array< CArray1N, D > prev_bases;
                    for (dim_t d = 0; d < D; d++) {
                        if (curr_index[d] == 0) {
                            //precursor is out of shape therefore set this precursor value to zero
                            prev_bases[d] = CArray1N::Zero(1,npts_);
                        }
                        else {
                            MultiIndex<D> prev_index = curr_index; prev_index[d] -= 1; //get backward neighbour
                            std::size_t prev_ordinal = enumeration_->slice(islice-2).find(prev_index);
                            assert (prev_ordinal < enumeration_->slice(islice-2).size()); //assert that multi-index has been found within previous slice
                            prev_bases[d] = prev_slice_values.row(prev_ordinal);
                        }
                    }
                    
                    //compute basis value within next slice
                    next_slice_values.row(j) = evaluateBasis(axis, curr_index, curr_basis, prev_bases);
                }
            }
        }
        
        static CMatrixCN evaluate_and_reduce_(const HagedornWavepacket &wavepacket, const CMatrixDN &x)
        {
            //Determine number of quadrature points. We cannot use template parameter N for
            //that as N could be 'Eigen::Dynamic'.
            int npts = x.cols();
            
            //Use Kahan's algorithm to accumulate bases with O(1) numerical error instead of O(Sqrt(N))
            KahanSum< CMatrixCN > psi( CMatrixCN::Zero(C,npts) );
            
            Evaluator<N> evaluator(wavepacket, x);
            
            CArrayXN prev_slice_values(0,npts);
            CArrayXN curr_slice_values(0,npts);
            CArrayXN next_slice_values(1,npts);
            
            evaluator.do_recursion(0, prev_slice_values, curr_slice_values, next_slice_values);
            
            psi += wavepacket.coefficient(0)*next_slice_values.row(0).matrix();
            
            for (std::size_t islice = 1; islice < wavepacket.enumeration_->count_slices(); islice++) {
                //exchange slices
                std::swap(prev_slice_values, curr_slice_values);
                std::swap(curr_slice_values, next_slice_values);
                next_slice_values = CArrayXN( wavepacket.enumeration_->slice(islice).size(), npts );
                
                evaluator.do_recursion(islice, prev_slice_values, curr_slice_values, next_slice_values);
                
                for (long j = 0; j < next_slice_values.rows(); j++) {
                    CMatrixC1 cj = wavepacket.coefficient(wavepacket.enumeration_->slice(islice).offset() + j);
                    psi += cj*next_slice_values.row(j).matrix();
                }
            }
            
            return psi();
        }
    };
    
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
        return Evaluator<N>::evaluate_and_reduce_(*this, xtmp);
    }
    
    /**
     * 
     */
    template<int N>
    Eigen::Matrix<complex_t,C,N> operator()(const Eigen::Matrix<complex_t,D,N> &x) const
    {
        return Evaluator<N>::evaluate_and_reduce_(*this, x);
    }
};

}

#endif