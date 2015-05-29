#ifndef WAVEBLOCKS_HAGEDORN_BASIS_EVALUATOR
#define WAVEBLOCKS_HAGEDORN_BASIS_EVALUATOR

#include <Eigen/Core>

#include "hagedorn_parameter_set.hpp"

#include "sliced_shape_enumeration.hpp"

namespace waveblocks {

/**
 * \tparam N number of quadrature points (if unknown: use Eigen::Dynamic)
 */
template<dim_t D, dim_t C, int N>
class Evaluator
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
        auto & p = parameters_->p;

        CMatrixDN P_Qinv_dx = P*Qinv_dx_;

        CArray1N pr1 = ( dx_.array() * P_Qinv_dx.array() ).colwise().sum();
        CArray1N pr2 = ( p.transpose()*dx_ ).array();

        CArray1N e = complex_t(0.0, 1.0)/(eps_*eps_) * (0.5*pr1 + pr2);

        return e.exp() / std::pow(pi<real_t>()*eps_*eps_, D/4.0);
    }

    /**
     * Complexity: theta(N*D)
     */
    CArray1N evaluateBasis(dim_t axis,
                           const std::array<int,D> &k,
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
    Evaluator(real_t eps, 
              std::shared_ptr< HagedornParameterSet<D> > parameters, 
              std::shared_ptr< ShapeEnumeration<D> > enumeration, 
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
                limit = std::max(limit, enumeration_->bbox(d));

            for (int i = 0; i <= limit+1; i++)
                sqrt_.push_back( std::sqrt( real_t(i) ) );
        }
    }

    /**
     * Complexity: theta(S*N*D)
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
                std::array<int,D> next_index = enumeration_->slice(islice)[j];
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
                std::array<int,D> curr_index = next_index;
                curr_index[axis] -= 1; //get backward neighbour
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
                        std::array<int,D> prev_index = curr_index;
                        prev_index[d] -= 1; //get backward neighbour
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
};

}

#endif