#ifndef WAVEBLOCKS_HAGEDORN_BASIS_EVALUATOR
#define WAVEBLOCKS_HAGEDORN_BASIS_EVALUATOR

#include <functional>

#include <Eigen/Core>

#include "hawp_paramset.hpp"

#include "shape_enum.hpp"
#include "kahan_sum.hpp"

namespace waveblocks {

/**
 * 
 * 
 * \tparam N number of quadrature points (if unknown: use Eigen::Dynamic)
 */
template<dim_t D, class MultiIndex, int N>
class HaWpEvaluator
{
public:
    /**
     * complex square matrix
     * size: (#dimensions)
     */
    typedef Eigen::Matrix<complex_t,D,D> CMatrixDD;

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

private:
    real_t eps_;
    const HaWpParamSet<D>* parameters_;
    const ShapeEnum<D,MultiIndex>* enumeration_;
    
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
    
public:
    HaWpEvaluator(real_t eps, 
              const HaWpParamSet<D>* parameters, 
              const ShapeEnum<D,MultiIndex>* enumeration,
              const CMatrixDN &x)
        : eps_(eps)
        , parameters_(parameters)
        , enumeration_(enumeration)
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
            int limit = 0;
            for (dim_t d = 0; d < D; d++)
                limit = std::max(limit, enumeration_->limit(d) );
            
            sqrt_.resize(limit+2);
            
            for (int i = 0; i <= limit+1; i++)
                sqrt_[i] = std::sqrt( real_t(i) );
        }
    }
    
    /**
     * \brief Evaluates value of the basis function with multi-index (0,0,...,0)
     * 
     * \return (1,N)-Matrix
     */
    CArray1N seed() const
    {
        auto & P = parameters_->P;
        auto & p = parameters_->p;

        CMatrixDN P_Qinv_dx = P*Qinv_dx_;

        CArray1N pr1 = ( dx_.array() * P_Qinv_dx.array() ).colwise().sum();
        CArray1N pr2 = ( p.transpose()*dx_ ).array();

        CArray1N e = complex_t(0.0, 1.0)/(eps_*eps_) * (0.5*pr1 + pr2);

        return e.exp() / std::pow(pi<real_t>()*eps_*eps_, D/4.0);
    }
    
    /**
     * \param[in] islice ordinal of current slice
     * \param[in] prev_basis basis values of previous slice
     * \param[in] curr_basis basis values of current slice
     * \return computed basis values of next slice
     */
    CArrayXN step(std::size_t islice,
                  const CArrayXN& prev_basis,
                  const CArrayXN& curr_basis) const
    {
        auto & prev_enum = enumeration_->slice(islice-1);
        auto & curr_enum = enumeration_->slice(islice);
        auto & next_enum = enumeration_->slice(islice+1);
        
        CArrayXN next_basis = CArrayXN::Zero(next_enum.size(), npts_);
        
        //loop over all multi-indices within next slice [j = position of multi-index within next slice]
        for (std::size_t j = 0; j < next_enum.size(); j++) {
            std::array<int,D> next_index = next_enum[j];
            //find valid precursor: find first non-zero entry
            dim_t axis = D;
            for (dim_t d = 0; d < D; d++) {
                if (next_index[d] != 0) {
                    axis = d;
                    break;
                }
            }
            
            assert(axis != D); //assert that multi-index contains some non-zero entries
            
            
            // compute contribution of current slice
            std::array<int,D> curr_index = next_index;
            curr_index[axis] -= 1; //get backward neighbour
            std::size_t curr_ordinal = curr_enum.find(curr_index);
            
            assert(curr_ordinal < curr_enum.size()); //assert that multi-index has been found within current slice
            
            CArray1N pr1 = curr_basis.row(curr_ordinal) * Qinv_dx_.row(axis).array() * std::sqrt(2.0)/eps_ ;
            
            
            // compute contribution of previous slice
            std::array< std::size_t,D > prev_ordinals = prev_enum.find_backward_neighbours(curr_index);
            
            CArray1N pr2{1,npts_};
            
            for (dim_t d = 0; d < D; d++) {
                if (curr_index[d] != 0) {
                    pr2 += prev_basis.row(prev_ordinals[d]) * Qh_Qinvt_(axis,d) * sqrt_[ curr_index[d] ];
                }
            }
            
            
            // compute basis value within next slice
            next_basis.row(j) = (pr1 - pr2) / sqrt_[ 1+curr_index[axis] ];
        }
        
        return next_basis;
    }
    
    /**
     * \brief evaluates wavepacket in a memory efficient manner
     */
    Eigen::Array<complex_t,1,N> reduce(const std::vector<complex_t>& coefficients)
    {
        // use Kahan's algorithm to accumulate bases with O(1) numerical error instead of O(Sqrt(N))
        KahanSum< Eigen::Array<complex_t,1,N> > psi( Eigen::Matrix<complex_t,1,N>::Zero(1,npts_) );
        
        Eigen::Array<complex_t,Eigen::Dynamic,N> prev_basis(0,npts_);
        Eigen::Array<complex_t,Eigen::Dynamic,N> curr_basis(0,npts_);
        Eigen::Array<complex_t,Eigen::Dynamic,N> next_basis(1,npts_);
        
        next_basis = seed();
        
        psi += coefficients[0]*next_basis.row(0).matrix();
        
        for (int islice = 0; islice < enumeration_->n_slices(); islice++) {
            prev_basis = std::move(curr_basis);
            curr_basis = std::move(next_basis);
            
            next_basis = step(islice, prev_basis, curr_basis);
            
            std::size_t offset = enumeration_->slice(islice+1).offset();
            
            for (long j = 0; j < next_basis.rows(); j++) {
                complex_t cj = coefficients[offset + j];
                
                //prints: multi-index -> basis -> coefficient
                //std::cout << enumeration->slice(islice+1)[j] << " -> " << next_basis.row(j).matrix() << " * " << cj << std::endl;
                
                psi += cj*next_basis.row(j).matrix();
            }
        }
        
        return psi();
    }
};

}

#endif